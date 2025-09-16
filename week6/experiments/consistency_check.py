# experiments/consistency_check.py
import csv
import os
import numpy as np

from radars import RADARS_ECEF
from sim.truth_provider import DynamicTruth
from sim.measurements import simulate_radar_measurements
from tracking.ekf import EKF

# Import the Blender-backed state function; try common locations to keep this script portable.
# Ensure one of these exists in the project; otherwise raise a clear error.
try:
    # Preferred: a dedicated bridge/module that exposes blender_state_fn(t)
    from sim.blender_bridge import blender_state_fn  # noqa: F401
except Exception:
    try:
        # Alternate: main simulation script exposes blender_state_fn(t)
        from missile_simulation_main import blender_state_fn  # noqa: F401
    except Exception:
        try:
            # Fallback: a generic main module
            from main import blender_state_fn  # noqa: F401
        except Exception as e:
            raise ImportError(
                "Could not import blender_state_fn from sim.blender_bridge, "
                "missile_simulation_main, or main. Ensure the Blender-side "
                "function blender_state_fn(t)->(p_ecef,v_ecef) is importable."
            ) from e

# Wire the truth provider to call Blender each time step without changing measurement/EKF code.
provider = DynamicTruth(blender_state_fn)


def _get_radar_ecef_by_id(rid):
    """
    Robustly fetch a radar ECEF by radar id from RADARS_ECEF, whether it is a list of dicts
    or a dict keyed by id -> { 'ecef_m': ... }.
    """
    # Case 1: list of radar dicts: [{'id': 'X', 'ecef_m': np.array([...])}, ...]
    if isinstance(RADARS_ECEF, (list, tuple)):
        for r in RADARS_ECEF:
            if isinstance(r, dict) and r.get("id") == rid:
                return r["ecef_m"]
        raise KeyError(f"Radar id {rid} not found in RADARS_ECEF list")
    # Case 2: dict keyed by id
    if isinstance(RADARS_ECEF, dict):
        if rid in RADARS_ECEF and isinstance(RADARS_ECEF[rid], dict) and "ecef_m" in RADARS_ECEF[rid]:
            return RADARS_ECEF[rid]["ecef_m"]
        # Degenerate single-anchor dict support (legacy)
        if "ecef_m" in RADARS_ECEF and not isinstance(RADARS_ECEF["ecef_m"], dict):
            return RADARS_ECEF["ecef_m"]
    raise TypeError("RADARS_ECEF must be a list/tuple of dicts or a dict keyed by radar id")


def run(duration_s=60.0, dt=0.5, qa=0.1,
        sigmas=(50.0, 0.0010, 0.0010, 3.0),
        out_csv="runs/run_consistency.csv"):
    """
    Consistency run that:
    - pulls ground-truth state from Blender via blender_state_fn(t),
    - simulates radar measurements,
    - runs the EKF predict/update,
    - writes a single CSV with truth, NIS, per-radar measurements, and R diag entries.
    """
    # Randomness for measurement noise, gating, and scheduling
    rng = np.random.default_rng(42)

    # Keep when a radar last measured for scheduling/cadence inside simulate_radar_measurements
    last_times = {}

    truth = provider  # DynamicTruth(blender_state_fn)

    # EKF init near first truth (t=0)
    p0, v0 = truth.step(0.0, dt)
    x0 = np.array([*p0, *v0], dtype=float)
    P0 = np.diag([1e4, 1e4, 1e4, 1e2, 1e2, 1e2])
    ekf = EKF(x0, P0, qa)

    # Measurement covariance diagonal from sigmas (range, az, el, dop)
    # Fix: compute elementwise square from tuple -> array
    R_diag = np.array(sigmas, dtype=float) ** 2

    # Ensure output directory exists
    out_dir = os.path.dirname(out_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        header = [
            "t",
            "px", "py", "pz",
            "vx", "vy", "vz",
            "nis",
            "radar_id",
            "range", "az", "el", "dop",
            "R11", "R22", "R33", "R44"
        ]
        w.writerow(header)

        t = 0.0
        # Small epsilon to ensure inclusion of final sample due to float addition
        while t <= duration_s + 1e-12:
            # Truth state from Blender-backed provider
            p, v = truth.step(t, dt)

            # EKF time update
            ekf.predict(dt)

            # Generate all radar measurements available at time t
            meas_list = simulate_radar_measurements(p, v, t, last_times, rng)

            # Process each measurement independently (same R for all; if per-radar R differs, adjust here)
            R = np.diag(R_diag)
            for meas in meas_list:
                rid = meas["id"]
                z = np.array(meas["z"], dtype=float)

                radar_ecef = _get_radar_ecef_by_id(rid)

                # EKF measurement update for this radar
                nis = ekf.update_radar(z, R, radar_ecef)

                # Persist a row per processed measurement
                w.writerow([t, *p, *v, nis, rid, *z, *R_diag])

            t += dt


if __name__ == "__main__":
    # Sensible default run for a quick end-to-end check
    run()
