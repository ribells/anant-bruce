# sim/measurements.py
from __future__ import annotations
import os
import sys
import math
import importlib.util
from typing import Dict, List, Tuple
import numpy as np

# ------------- CONFIG: set your main sim path here -------------
# Example: MAIN_PATH = "/absolute/path/to/main_sim.py"
MAIN_PATH = os.environ.get("SIM_MAIN_PATH", "").strip()
# ---------------------------------------------------------------

# Helper to dynamically import a module from a file path
# Ref: importlib util recipe (Python docs)
def import_from_path(module_name: str, file_path: str):
    spec = importlib.util.spec_from_file_location(module_name, file_path)  # [2]
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {module_name} from {file_path}")
    mod = importlib.util.module_from_spec(spec)                            # [2]
    sys.modules[module_name] = mod
    spec.loader.exec_module(mod)                                           # [2]
    return mod

# Import radars (must be importable via PYTHONPATH or package)
from radars import RADARS_ECEF, should_radar_emit  # [9]

PI = math.pi

def wrap_pi(a: float) -> float:
    return (a + PI) % (2.0 * PI) - PI

def clamp_el(e: float) -> float:
    half = 0.5 * PI
    if e > half: return half
    if e < -half: return -half
    return e

def geom_to_measurement(
    p_ecef: Tuple[float, float, float],
    v_ecef: Tuple[float, float, float],
    r_ecef: Tuple[float, float, float],
) -> Tuple[float, float, float, float]:
    """
    Compute ideal [range, az, el, doppler] from ECEF geometry.
    """
    px, py, pz = p_ecef
    rx, ry, rz = r_ecef
    dx, dy, dz = px - rx, py - ry, pz - rz
    R = math.sqrt(dx*dx + dy*dy + dz*dz)
    if R <= 1e-6:
        # Degenerate: return zeros (caller can skip)
        return 0.0, 0.0, 0.0, 0.0
    az = math.atan2(dy, dx)
    el = math.asin(dz / R)
    # Doppler: radial component of velocity along line-of-sight
    vx, vy, vz = v_ecef
    losx, losy, losz = dx / R, dy / R, dz / R
    dop = losx * vx + losy * vy + losz * vz
    return R, az, el, dop

def simulate_radar_measurements(
    p_ecef: Tuple[float, float, float],
    v_ecef: Tuple[float, float, float],
    t_now: float,
    last_times: Dict[str, float],
    rng: np.random.Generator,
) -> List[Dict]:
    """
    Generate noisy radar measurements for all radars that emit at t_now.

    Inputs:
      - p_ecef, v_ecef: target position/velocity in ECEF (m, m/s)
      - t_now: current sim time (s)
      - last_times: dict radar_id -> last emission time (s), updated in-place
      - rng: numpy random Generator for reproducibility

    Output list items:
      { "id": radar_id, "z": np.array([range, az, el, dop]),
        "R_diag": [σ_r², σ_az², σ_el², σ_dop²],
        "true": (R_true, az_true, el_true, dop_true) }
    """
    outputs: List[Dict] = []
    for radar in RADARS_ECEF:  # [9]
        rid = radar["id"]
        rate = radar["update_rate_hz"]
        if not should_radar_emit(last_times.get(rid, -1e12), t_now, rate):  # [9]
            continue

        R_true, az_true, el_true, dop_true = geom_to_measurement(p_ecef, v_ecef, radar["ecef_m"])

        # Skip if degenerate geometry (very rare)
        if R_true <= 1e-6:
            continue

        # Draw noise using per-radar sigmas
        sr = radar["sigma_range_m"]
        saz = radar["sigma_az_rad"]
        sel = radar["sigma_el_rad"]
        sd = radar["sigma_doppler_mps"]

        Rn  = rng.normal(0.0, sr)
        azn = rng.normal(0.0, saz)
        eln = rng.normal(0.0, sel)
        dn  = rng.normal(0.0, sd)

        z_range = R_true + Rn
        z_az    = wrap_pi(az_true + azn)
        z_el    = clamp_el(el_true + eln)
        z_dop   = dop_true + dn

        outputs.append({
            "id": rid,
            "z": np.array([z_range, z_az, z_el, z_dop], dtype=float),
            "R_diag": radar["R_diag"],
            "true": (R_true, az_true, el_true, dop_true),
        })

        # Update emission time
        last_times[rid] = t_now

    return outputs

def main_entry():
    """
    Example entry which loads your main sim module dynamically using MAIN_PATH
    and runs a single measurement step (for demonstration).
    Set SIM_MAIN_PATH env var or edit MAIN_PATH above.
    """
    if not MAIN_PATH:
        raise RuntimeError("Set SIM_MAIN_PATH environment variable to the path of your main sim file, or edit MAIN_PATH in measurements.py")

    # Import the main module dynamically (no refactor needed)  [2][3]
    main_mod = import_from_path("sim_main_module", MAIN_PATH)

    # Expect main module to expose current truth or a getter for (p_ecef, v_ecef, t_now)
    # Adjust names to your actual main API.
    if hasattr(main_mod, "get_truth_state_ecef"):
        p_ecef, v_ecef, t_now = main_mod.get_truth_state_ecef()
    else:
        # Example fallback: look for variables
        p_ecef = getattr(main_mod, "P_ECEF", (0.0, 0.0, 0.0))
        v_ecef = getattr(main_mod, "V_ECEF", (0.0, 0.0, 0.0))
        t_now  = getattr(main_mod, "T_NOW", 0.0)

    # Seed RNG for reproducibility in experiments
    rng = np.random.default_rng(12345)

    # Track last emission times
    last_times: Dict[str, float] = {}

    # Generate measurements for this time
    meas = simulate_radar_measurements(p_ecef, v_ecef, t_now, last_times, rng)

    # Example: print a quick summary
    for m in meas:
        rid = m["id"]
        z = m["z"]
        print(f"{rid}: range={z:.1f} m, az={z:.4f} rad, el={z:.4f} rad, dop={z:.2f} m/s")

if __name__ == "__main__":
    main_entry()
