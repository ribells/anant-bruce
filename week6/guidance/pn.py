# tracking/pn.py
from __future__ import annotations
import numpy as np

def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-9 else v

def los_vector(target_pos: np.ndarray, interceptor_pos: np.ndarray) -> np.ndarray:
    """Line-of-sight unit vector from interceptor to target."""
    return _unit(target_pos - interceptor_pos)

def los_rate(target_pos: np.ndarray, target_vel: np.ndarray,
             interceptor_pos: np.ndarray, interceptor_vel: np.ndarray) -> np.ndarray:
    """
    LOS rate vector (rad/s) approximation in 3D:
    ω_LOS ≈ (r × r_dot) / ||r||^2, where r = target_pos - interceptor_pos.
    """
    r = target_pos - interceptor_pos
    rdot = target_vel - interceptor_vel
    rnorm2 = float(np.dot(r, r))
    if rnorm2 < 1e-6:
        return np.zeros(3)
    return np.cross(r, rdot) / rnorm2

def pn_command(target_pos: np.ndarray, target_vel: np.ndarray,
               interceptor_pos: np.ndarray, interceptor_vel: np.ndarray,
               N: float = 3.0) -> np.ndarray:
    """
    Classical proportional navigation guidance acceleration command in m/s^2:
    a_cmd = N * Vc * (ω_LOS × v_hat)
    where Vc is closing speed, ω_LOS is LOS-rate, and v_hat is interceptor velocity unit vector.
    """
    v_hat = _unit(interceptor_vel)
    rdot = float(np.dot(target_vel - interceptor_vel, _unit(target_pos - interceptor_pos)))
    Vc = -rdot  # closing speed (positive if closing)
    omega_los = los_rate(target_pos, target_vel, interceptor_pos, interceptor_vel)
    a_cmd = N * Vc * np.cross(omega_los, v_hat)
    return a_cmd

def limited_accel(a_cmd: np.ndarray, a_max: float) -> np.ndarray:
    """Clip commanded acceleration to interceptor physical limits."""
    n = float(np.linalg.norm(a_cmd))
    if n <= a_max:
        return a_cmd
    return a_cmd * (a_max / n)

def evasion_vector(target_pos: np.ndarray, target_vel: np.ndarray,
                   interceptor_pos: np.ndarray) -> np.ndarray:
    """
    Evasion direction heuristic: perpendicular to LOS and target velocity, in the plane of max escape.
    """
    r = target_pos - interceptor_pos
    u_los = _unit(r)
    v_perp = target_vel - np.dot(target_vel, u_los) * u_los
    if np.linalg.norm(v_perp) < 1e-6:
        # pick any perpendicular direction
        basis = np.array([1.0, 0.0, 0.0])
        v_perp = np.cross(u_los, basis)
        if np.linalg.norm(v_perp) < 1e-6:
            basis = np.array([0.0, 1.0, 0.0])
            v_perp = np.cross(u_los, basis)
    return _unit(v_perp)

def evasive_accel(target_pos: np.ndarray, target_vel: np.ndarray,
                  interceptor_pos: np.ndarray, a_max: float,
                  scale: float = 1.0) -> np.ndarray:
    """
    Evasive acceleration at up to a_max, directed to maximize miss distance based on LOS geometry.
    """
    dir_evade = evasion_vector(target_pos, target_vel, interceptor_pos)
    return dir_evade * (a_max * scale)
