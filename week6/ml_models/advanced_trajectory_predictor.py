# ml_models/advanced_trajectory_predictor.py
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import numpy as np
from typing import Tuple, List, Dict, Optional
import h5py

class MultiModalMissileDataset(Dataset):
    """
    Dataset for loading multi-modal missile signature data
    """
    
    def __init__(self, csv_path: str, h5_path: str, sequence_length: int = 50):
        self.csv_data = pd.read_csv(csv_path)
        self.h5_file = h5py.File(h5_path, 'r')
        self.sequence_length = sequence_length
        
        # Group data by missile track
        self.tracks = self._group_by_track()
        
    def _group_by_track(self) -> List[List[int]]:
        """Group measurements by continuous tracks"""
        tracks = []
        current_track = []
        
        for i, row in self.csv_data.iterrows():
            if i == 0 or (row['t'] - self.csv_data.iloc[i-1]['t']) > 5.0:
                # New track (time gap > 5 seconds)
                if current_track:
                    tracks.append(current_track)
                current_track = [i]
            else:
                current_track.append(i)
        
        if current_track:
            tracks.append(current_track)
            
        return tracks
    
    def __len__(self) -> int:
        # Number of valid sequences across all tracks
        return sum(max(0, len(track) - self.sequence_length) for track in self.tracks)
    
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        # Find which track and position within track
        track_idx, seq_start = self._get_track_position(idx)
        track = self.tracks[track_idx]
        
        # Extract sequence
        sequence_indices = track[seq_start:seq_start + self.sequence_length]
        
        # Time series features
        ts_features = []
        hrr_profiles = []
        micro_doppler = []
        
        for i in sequence_indices:
            row = self.csv_data.iloc[i]
            
            # Time series: position, velocity, RCS, SNR, quality
            ts_features.append([
                row['px'], row['py'], row['pz'],
                row['vx'], row['vy'], row['vz'],
                row['rcs_dbsm'], row['snr_db'], row['quality_score'],
                row['velocity_mach'], row['altitude_m']
            ])
            
            # Load signatures from HDF5
            sig_id = row['signature_id']
            hrr_profiles.append(self.h5_file['hrr_profiles'][sig_id][:])
            micro_doppler.append(self.h5_file['micro_doppler'][sig_id][:])
        
        # Convert to tensors
        ts_tensor = torch.tensor(ts_features, dtype=torch.float32)
        hrr_tensor = torch.tensor(np.stack(hrr_profiles), dtype=torch.float32)
        md_tensor = torch.tensor(np.stack(micro_doppler), dtype=torch.float32)
        
        # Target: next position/velocity (for trajectory prediction)
        if seq_start + self.sequence_length < len(track):
            next_row = self.csv_data.iloc[track[seq_start + self.sequence_length]]
            target = torch.tensor([
                next_row['px'], next_row['py'], next_row['pz'],
                next_row['vx'], next_row['vy'], next_row['vz']
            ], dtype=torch.float32)
        else:
            # Last available state (for terminal prediction)
            target = ts_tensor[-1, :6]
            
        return ts_tensor, hrr_tensor, md_tensor, target

class AdvancedTrajectoryPredictor(nn.Module):
    """
    Multi-modal transformer-based trajectory predictor with real-time adaptation
    """
    
    def __init__(self, config: Dict):
        super().__init__()
        self.config = config
        
        # Time series encoder
        self.ts_encoder = nn.LSTM(
            input_size=11,  # px,py,pz,vx,vy,vz,rcs,snr,quality,mach,alt
            hidden_size=config['ts_hidden'],
            num_layers=config['ts_layers'],
            batch_first=True,
            dropout=config.get('dropout', 0.1)
        )
        
        # HRR profile encoder (1D CNN)
        self.hrr_encoder = nn.Sequential(
            nn.Conv1d(config['sequence_length'], 64, kernel_size=7, padding=3),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(64, 128, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(128, 256, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.AdaptiveAvgPool1d(32),
            nn.Flatten()
        )
        
        # Micro-doppler encoder (2D CNN)
        self.md_encoder = nn.Sequential(
            nn.Conv3d(1, 32, kernel_size=(3, 3, 3), padding=(1, 1, 1)),
            nn.ReLU(),
            nn.MaxPool3d((2, 2, 2)),
            nn.Conv3d(32, 64, kernel_size=(3, 3, 3), padding=(1, 1, 1)),
            nn.ReLU(),
            nn.MaxPool3d((2, 2, 2)),
            nn.AdaptiveAvgPool3d((4, 4, 4)),
            nn.Flatten()
        )
        
        # Attention mechanism for feature fusion
        self.attention = nn.MultiheadAttention(
            embed_dim=config['fusion_dim'],
            num_heads=config['num_heads'],
            dropout=config.get('dropout', 0.1)
        )
        
        # Feature projection layers
        self.ts_proj = nn.Linear(config['ts_hidden'], config['fusion_dim'])
        self.hrr_proj = nn.Linear(256 * 32, config['fusion_dim'])
        self.md_proj = nn.Linear(64 * 64, config['fusion_dim'])
        
        # Multi-task prediction heads
        self.trajectory_head = nn.Sequential(
            nn.Linear(config['fusion_dim'], 256),
            nn.ReLU(),
            nn.Dropout(config.get('dropout', 0.1)),
            nn.Linear(256, 6)  # px,py,pz,vx,vy,vz
        )
        
        self.maneuver_head = nn.Sequential(
            nn.Linear(config['fusion_dim'], 128),
            nn.ReLU(),
            nn.Linear(128, 4)  # cruise, maneuver, terminal, evasive
        )
        
        self.uncertainty_head = nn.Sequential(
            nn.Linear(config['fusion_dim'], 64),
            nn.ReLU(),
            nn.Linear(64, 6)  # uncertainty for each trajectory component
        )
        
        # Real-time adaptation parameters
        self.adaptation_rate = nn.Parameter(torch.tensor(0.01))
        self.context_memory = nn.Parameter(torch.zeros(config['fusion_dim']))
        
    def forward(self, ts_seq: torch.Tensor, hrr_seq: torch.Tensor, 
                md_seq: torch.Tensor) -> Dict[str, torch.Tensor]:
        
        batch_size = ts_seq.size(0)
        
        # Encode time series
        ts_out, (h_n, c_n) = self.ts_encoder(ts_seq)
        ts_features = self.ts_proj(h_n[-1])  # Use last hidden state
        
        # Encode HRR profiles
        hrr_flat = hrr_seq.view(batch_size, self.config['sequence_length'], -1)
        hrr_out = self.hrr_encoder(hrr_flat)
        hrr_features = self.hrr_proj(hrr_out)
        
        # Encode micro-doppler (add sequence dimension as channels)
        md_input = md_seq.unsqueeze(1)  # Add channel dimension
        md_out = self.md_encoder(md_input)
        md_features = self.md_proj(md_out)
        
        # Stack features for attention
        features = torch.stack([ts_features, hrr_features, md_features], dim=0)
        
        # Self-attention fusion
        fused_features, attention_weights = self.attention(features, features, features)
        fused_features = fused_features.mean(dim=0)  # Average across modalities
        
        # Real-time adaptation: blend with context memory
        adapted_features = (1 - self.adaptation_rate) * fused_features + \
                          self.adaptation_rate * self.context_memory.expand_as(fused_features)
        
        # Multi-task predictions
        trajectory_pred = self.trajectory_head(adapted_features)
        maneuver_pred = F.softmax(self.maneuver_head(adapted_features), dim=-1)
        uncertainty_pred = F.softplus(self.uncertainty_head(adapted_features))
        
        return {
            'trajectory': trajectory_pred,
            'maneuver_probabilities': maneuver_pred,
            'uncertainty': uncertainty_pred,
            'attention_weights': attention_weights,
            'adapted_features': adapted_features
        }
    
    def update_context_memory(self, new_context: torch.Tensor, learning_rate: float = 0.1):
        """Real-time adaptation of context memory"""
        with torch.no_grad():
            self.context_memory.data = (1 - learning_rate) * self.context_memory.data + \
                                     learning_rate * new_context.mean(dim=0)

class RealTimeAdaptationManager:
    """
    Manages real-time model adaptation based on tracking performance
    """
    
    def __init__(self, model: AdvancedTrajectoryPredictor, adaptation_config: Dict):
        self.model = model
        self.config = adaptation_config
        
        # Performance tracking
        self.recent_errors = []
        self.adaptation_threshold = adaptation_config['error_threshold']
        self.window_size = adaptation_config['window_size']
        
        # Online learning components
        self.optimizer = torch.optim.Adam(model.parameters(), lr=adaptation_config['learning_rate'])
        self.adaptation_loss = nn.MSELoss()
        
    def should_adapt(self, prediction_error: float) -> bool:
        """Determine if model should adapt based on recent performance"""
        self.recent_errors.append(prediction_error)
        
        if len(self.recent_errors) > self.window_size:
            self.recent_errors.pop(0)
        
        if len(self.recent_errors) >= self.window_size:
            mean_error = np.mean(self.recent_errors)
            return mean_error > self.adaptation_threshold
        
        return False
    
    def adapt_model(self, recent_data: torch.Tensor, true_outcomes: torch.Tensor):
        """Perform online learning adaptation"""
        self.model.train()
        
        # Forward pass
        predictions = self.model(recent_data['ts'], recent_data['hrr'], recent_data['md'])
        
        # Compute adaptation loss
        traj_loss = self.adaptation_loss(predictions['trajectory'], true_outcomes)
        
        # Backward pass
        self.optimizer.zero_grad()
        traj_loss.backward()
        
        # Gradient clipping for stability
        torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
        
        self.optimizer.step()
        
        # Update context memory
        self.model.update_context_memory(predictions['adapted_features'])
        
        self.model.eval()
        
        return traj_loss.item()

# Training data generation
def generate_synthetic_training_data(n_scenarios: int = 10000, 
                                   output_dir: str = "training_data/") -> str:
    """
    Generate synthetic training scenarios covering various missile types and scenarios
    """
    from trajectory_generators.supermaneuverable_generator import SupermaneuverableTrajectoryGenerator
    
    scenarios = []
    
    missile_types = ['ballistic', 'cruise', 'hgv', 'supermaneuverable']
    flight_phases = ['boost', 'cruise', 'terminal', 'evasive']
    
    for scenario_id in range(n_scenarios):
        # Random scenario parameters
        missile_type = np.random.choice(missile_types)
        has_threats = np.random.random() < 0.3  # 30% scenarios have threats
        
        # Generate trajectory
        start_lat = np.random.uniform(35, 55)  # Approximate threat region
        start_lon = np.random.uniform(-120, 50)
        target_lat = np.random.uniform(35, 55)
        target_lon = np.random.uniform(-120, 50)
        
        threat_locations = []
        if has_threats:
            n_threats = np.random.randint(1, 4)
            for _ in range(n_threats):
                threat_locations.append((
                    np.random.uniform(35, 55),
                    np.random.uniform(-120, 50)
                ))
        
        # Generate trajectory points
        trajectory_points = generate_dynamic_trajectory_blender(
            start_lat, start_lon, target_lat, target_lon, threat_locations, n_points=300
        )
        
        # Simulate measurements along trajectory
        dt = 1.0  # 1 second intervals
        scenario_data = []
        
        for i, point in enumerate(trajectory_points[:-1]):
            t = i * dt
            next_point = trajectory_points[i + 1]
            
            # Estimate velocity
            velocity = (next_point - point) / dt
            
            # Determine flight phase based on time and position
            total_time = len(trajectory_points) * dt
            if t < total_time * 0.1:
                phase = 'boost'
            elif t < total_time * 0.8:
                phase = 'cruise'
            else:
                phase = 'terminal'
            
            # Add evasive phase if threats present
            if has_threats and phase == 'terminal':
                if np.random.random() < 0.3:  # 30% chance of evasive maneuver
                    phase = 'evasive'
            
            missile_context = {
                'missile_type': missile_type,
                'flight_state': phase,
                'current_g_load': np.random.uniform(0, 12),
                'spin_hz': np.random.uniform(0, 2)
            }
            
            # Simulate radar measurements
            measurements = simulate_advanced_radar_measurements(
                point, velocity, t, {}, np.random.default_rng(), missile_context
            )
            
            scenario_data.extend(measurements)
        
        scenarios.append({
            'scenario_id': scenario_id,
            'missile_type': missile_type,
            'has_threats': has_threats,
            'measurements': scenario_data
        })
    
    # Save to files
    import pickle
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    with open(output_path / 'training_scenarios.pkl', 'wb') as f:
        pickle.dump(scenarios, f)
    
    print(f"Generated {n_scenarios} training scenarios")
    return str(output_path)
