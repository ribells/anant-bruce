# experiments/analyze_blender_trajectory.py
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))  # add project root

import numpy as np
import matplotlib.pyplot as plt
from filters.ekf import EKF
from sim.measurements import simulate_radar_measurements  
from radars import RADARS_ECEF

def load_blender_trajectory():
    """Extract trajectory from your Blender script logic"""
    # Copy the trajectory computation from your Blender script
    # ... (ll_to_ecef, great_circle_points, altitude_profile functions)
    
    # Recompute the same trajectory
    LAUNCH_LAT, LAUNCH_LON = 40.730610, -73.935242
    TARGET_LAT, TARGET_LON = 38.907192, -77.036871
    # ... rest of trajectory computation
    
    return world_pts_m, dt_simulation

def run_ekf_analysis():
    world_pts_m, dt_sim = load_blender_trajectory()
    
    # Same EKF setup as Option 1...
    # ... EKF initialization and tracking loop
    
    # Analysis and plotting
    true_positions = np.array(world_pts_m)
    estimated_positions = np.array([est[:3] for est in ekf_estimates])
    
    # Compute RMSE
    position_errors = np.linalg.norm(true_positions - estimated_positions, axis=1)
    rmse = np.sqrt(np.mean(position_errors**2))
    
    # Plots
    fig, ((ax1,ax2), (ax3,ax4)) = plt.subplots(2,2, figsize=(12,8))
    
    # Position error over time
    t_vec = np.arange(len(position_errors)) * dt_sim
    ax1.plot(t_vec, position_errors)
    ax1.set_ylabel('Position Error (m)')
    ax1.set_title(f'RMSE: {rmse:.1f} m')
    
    # NIS consistency
    ax2.plot(t_vec, nis_values)
    ax2.axhline(y=4, color='r', linestyle='--', label='Expected (4 DOF)')
    ax2.set_ylabel('NIS')
    ax2.set_title('Filter Consistency')
    ax2.legend()
    
    # 3D trajectory comparison
    ax3 = fig.add_subplot(223, projection='3d')
    ax3.plot(*true_positions.T, 'b-', label='True', alpha=0.7)
    ax3.plot(*estimated_positions.T, 'r--', label='EKF', alpha=0.7)
    ax3.legend()
    ax3.set_title('Trajectory Comparison')
    
    # Error breakdown by axis
    errors_xyz = true_positions - estimated_positions
    ax4.plot(t_vec, errors_xyz)
    ax4.legend(['X error', 'Y error', 'Z error'])
    ax4.set_ylabel('Position Error (m)')
    ax4.set_xlabel('Time (s)')
    
    plt.tight_layout()
    plt.savefig('runs/ekf_blender_analysis.png', dpi=150)
    plt.show()
    
    print(f"Analysis complete:")
    print(f"  RMSE: {rmse:.1f} m")
    print(f"  Mean NIS: {np.mean(nis_values):.2f} (expect ~4.0)")
    print(f"  NIS in 95% bounds: {np.mean((nis_values > 0.48) & (nis_values < 11.14))*100:.1f}%")

if __name__ == "__main__":
    run_ekf_analysis()
