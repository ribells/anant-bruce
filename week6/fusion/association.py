# tracking/association.py
from __future__ import annotations
import numpy as np
from typing import Dict, List, Tuple, Optional
from scipy.optimize import linear_sum_assignment

def gating_mahalanobis(z: np.ndarray, z_pred: np.ndarray, S: np.ndarray, gamma: float) -> bool:
    """
    Ellipsoidal gating using Mahalanobis distance: (z - z_pred)^T S^{-1} (z - z_pred) <= gamma
    """
    dz = z - z_pred
    try:
        Sinv = np.linalg.inv(S)
    except np.linalg.LinAlgError:
        return False
    d2 = float(dz.T @ Sinv @ dz)
    return d2 <= gamma

def innovation_distance(z: np.ndarray, z_pred: np.ndarray, S: np.ndarray) -> float:
    dz = z - z_pred
    try:
        Sinv = np.linalg.inv(S)
    except np.linalg.LinAlgError:
        return np.inf
    return float(dz.T @ Sinv @ dz)

def rcs_cost(rcs_meas_db: float, rcs_track_db: float, sigma_db: float = 5.0) -> float:
    """Gaussian cost in dB for RCS consistency."""
    d = (rcs_meas_db - rcs_track_db) / sigma_db
    return d*d

def snr_weight(snr_db: float, min_db: float = -10.0, max_db: float = 30.0) -> float:
    """
    Map SNR(dB) to [0,1] to scale measurement reliability.
    """
    return float(np.clip((snr_db - min_db) / (max_db - min_db), 0.0, 1.0))

def feature_cost_hrr(hrr_meas: np.ndarray, hrr_ref: np.ndarray) -> float:
    """
    Simple cosine distance between HRR profiles.
    """
    a = hrr_meas.astype(np.float32); b = hrr_ref.astype(np.float32)
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    if na < 1e-6 or nb < 1e-6:
        return 1.0
    cos = float(np.dot(a, b) / (na * nb))
    return 1.0 - cos  # 0 best, 1 worst

def feature_cost_md(md_meas: np.ndarray, md_ref: np.ndarray) -> float:
    """
    L2 distance on micro-doppler spectrograms after normalization.
    """
    a = md_meas.astype(np.float32); b = md_ref.astype(np.float32)
    a /= (np.linalg.norm(a) + 1e-6)
    b /= (np.linalg.norm(b) + 1e-6)
    return float(np.linalg.norm(a - b))

def build_cost_matrix(
    tracks: List[Dict],
    measurements: List[Dict],
    kinematic_models: Dict[int, Dict],
    w_kin: float = 1.0,
    w_rcs: float = 0.2,
    w_hrr: float = 0.3,
    w_md: float  = 0.3,
) -> np.ndarray:
    """
    Build an association cost matrix combining kinematics with advanced features.
    Each track dict should hold:
      - 'z_pred': predicted measurement vector (4,)
      - 'S': innovation covariance (4x4)
      - 'rcs_ref_db': running estimate of RCS in dB
      - 'hrr_ref': reference HRR profile (1D)
      - 'md_ref': reference Micro-Doppler (2D)
    Each measurement dict may include:
      - 'z', 'rcs_dbsm', 'hrr_profile', 'micro_doppler', 'snr_db', 'quality'
    """
    T = len(tracks); M = len(measurements)
    C = np.full((T, M), fill_value=1e6, dtype=float)

    for i, trk in enumerate(tracks):
        z_pred = trk.get('z_pred', None)
        S = trk.get('S', None)
        if z_pred is None or S is None:
            continue

        for j, meas in enumerate(measurements):
            z = np.array(meas["z"], dtype=float)
            # Kinematic innovation distance
            kin = innovation_distance(z, z_pred, S)
            if not np.isfinite(kin):
                continue

            # Feature costs
            rcs_m = float(meas.get("rcs_dbsm", 0.0))
            rcs_t = float(trk.get("rcs_ref_db", rcs_m))
            c_rcs = rcs_cost(rcs_m, rcs_t)

            hrr_m = np.array(meas.get("hrr_profile", []), dtype=np.float32)
            hrr_t = np.array(trk.get("hrr_ref", hrr_m), dtype=np.float32)
            c_hrr = feature_cost_hrr(hrr_m, hrr_t) if hrr_m.size and hrr_t.size else 0.0

            md_m = np.array(meas.get("micro_doppler", []), dtype=np.float32)
            md_t = np.array(trk.get("md_ref", md_m), dtype=np.float32)
            c_md = feature_cost_md(md_m, md_t) if md_m.size and md_t.size else 0.0

            # SNR-based weighting (reliability)
            w_meas = snr_weight(float(meas.get("snr_db", 0.0)))

            # Combined cost (lower is better)
            C[i, j] = w_kin * kin + w_meas * (w_rcs * c_rcs + w_hrr * c_hrr + w_md * c_md)

    return C

def associate_tracks_measurements(
    tracks: List[Dict],
    measurements: List[Dict],
    kinematic_models: Dict[int, Dict],
    gating_threshold: float = 16.0,  # ~chi-square 95% for 4 DOF
    **cost_kwargs,
) -> Tuple[List[Tuple[int,int]], List[int], List[int]]:
    """
    Returns:
      - matches: list of (track_index, measurement_index)
      - unassigned_tracks: list of track indices
      - unassigned_measurements: list of measurement indices
    """
    if not tracks or not measurements:
        return [], list(range(len(tracks))), list(range(len(measurements)))

    # Pre-gate to mask impossible pairs
    gating_mask = np.zeros((len(tracks), len(measurements)), dtype=bool)
    for i, trk in enumerate(tracks):
        z_pred = trk.get('z_pred', None)
        S = trk.get('S', None)
        if z_pred is None or S is None:
            continue
        for j, meas in enumerate(measurements):
            z = np.array(meas["z"], dtype=float)
            gating_mask[i, j] = gating_mahalanobis(z, z_pred, S, gating_threshold)

    # Build cost matrix
    C = build_cost_matrix(tracks, measurements, kinematic_models, **cost_kwargs)

    # Apply gating mask: inf cost where gated out
    C_masked = C.copy()
    C_masked[~gating_mask] = 1e9

    # Hungarian assignment
    row_ind, col_ind = linear_sum_assignment(C_masked)
    matches = []
    assigned_tracks = set()
    assigned_meas = set()

    for r, c in zip(row_ind, col_ind):
        if C_masked[r, c] < 1e8:  # only accept non-gated pairs
            matches.append((r, c))
            assigned_tracks.add(r)
            assigned_meas.add(c)

    unassigned_tracks = [i for i in range(len(tracks)) if i not in assigned_tracks]
    unassigned_measurements = [j for j in range(len(measurements)) if j not in assigned_meas]

    return matches, unassigned_tracks, unassigned_measurements
