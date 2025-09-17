# sim/measurements.py
from __future__ import annotations
import os
import sys
import math
import importlib.util
from typing import Dict, List, Tuple, Optional
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

def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-9 else v

def _aspect_angle_rad(p_ecef: np.ndarray, v_ecef: np.ndarray, r_ecef: np.ndarray) -> float:
    """
    Aspect angle between radar line-of-sight (radar->target) and target body axis proxy (velocity).
    """
    los = np.array(p_ecef) - np.array(r_ecef)
    u_los = _unit(los)
    u_vel = _unit(np.array(v_ecef))
    c = float(np.clip(np.dot(u_los, u_vel), -1.0, 1.0))
    return math.acos(c)

# === Advanced measurement surrogates (fast, deterministic, tunable) ===

def _simulate_rcs_dbsm(missile_type: str, aspect_angle_rad: float,
                       altitude_m: float, flight_state: str, rng: Optional[np.random.Generator]) -> float:
    """
    Very lightweight RCS surrogate (dBsm) keyed by missile type and aspect.
    Tuned for a supersonic cruise missile baseline.
    """
    # Baselines (Brahmos-like defaults)
    base = {
        "supersonic_cruise": -3.0,
        "ballistic": -12.0,
        "hgv": -15.0,
        "booster": 8.0,
        "warhead": -22.0
    }.get(missile_type, -8.0)

    # Aspect dependence: nose-on low, broadside higher
    aspect_gain = 8.0 * abs(math.sin(aspect_angle_rad))

    # Altitude dependence (multipath at very low altitude)
    if altitude_m < 100:
        alt_adj = 2.0
    elif altitude_m > 15000:
        alt_adj = -1.0
    else:
        alt_adj = 0.0

    # Flight state modifier
    state_adj = {
        "boost": 4.0,
        "cruise": 0.0,
        "terminal": 1.5,
        "maneuver": 3.0,
        "evasive": 3.0
    }.get(flight_state, 0.0)

    jitter = (rng.normal(0, 0.3) if rng is not None else 0.0)  # small fluctuation
    return base + aspect_gain + alt_adj + state_adj + jitter

def _simulate_snr_db(range_m: float, rcs_dbsm: float,
                     radar_power_db: float = 100.0,
                     noise_floor_db: float = -110.0,
                     processing_gain_db: float = 20.0,
                     path_loss_coeff: float = 40.0) -> float:
    """
    Simplified SNR in dB: SNR ≈ P + RCS - PL(range) - NF + PG
    """
    rk = max(range_m / 1000.0, 0.1)
    path_loss_db = path_loss_coeff * math.log10(rk)
    snr_db = radar_power_db + rcs_dbsm - path_loss_db - noise_floor_db + processing_gain_db
    return max(snr_db, -20.0)

def _simulate_hrr_profile(missile_type: str, aspect_angle_rad: float,
                          velocity_mach: float, flight_state: str,
                          n_bins: int = 128, rng: Optional[np.random.Generator] = None) -> np.ndarray:
    rng = rng or np.random.default_rng()
    prof = np.zeros(n_bins, dtype=np.float32)

    if missile_type == "booster":
        a = slice(n_bins//3, n_bins//3+10)
        b = slice(2*n_bins//3-5, 2*n_bins//3+5)
        prof[a] += 1.0
        prof[b] += 0.6
    elif missile_type == "warhead":
        c = slice(n_bins//2-2, n_bins//2+2)
        prof[c] += 1.0
    elif missile_type == "supersonic_cruise":
        # Major scattering centers: nose, body junction, intake
        nose = int(0.15*n_bins); prof[nose-2:nose+3] = 1.0
        body = int(0.4*n_bins);  prof[body-3:body+4] = 0.7
        intake = int(0.75*n_bins)
        prof[intake-2:intake+3] = 0.6 * abs(math.cos(aspect_angle_rad))
        if flight_state in ("maneuver", "evasive", "terminal"):
            # More spread under dynamic states
            span = 6
            for center in [int(0.25*n_bins), int(0.55*n_bins), int(0.8*n_bins)]:
                prof[max(center-span,0):min(center+span,n_bins)] += 0.2
    else:
        # Generic elongated body
        span = 14
        s = n_bins//2 - span//2
        prof[s:s+span] += np.hanning(span).astype(np.float32)

    noise = 0.05 if flight_state in ("maneuver","evasive","terminal") else 0.02
    prof += noise * (rng.standard_normal(n_bins).astype(np.float32))
    prof = np.clip(prof, 0.0, None)
    return prof

def _simulate_micro_doppler(missile_type: str, flight_state: str,
                            spin_hz: float = 0.0, H: int = 64, W: int = 64,
                            rng: Optional[np.random.Generator] = None) -> np.ndarray:
    rng = rng or np.random.default_rng()
    md = np.zeros((H, W), dtype=np.float32)
    mid = H // 2
    if missile_type == "warhead":
        amp = max(2, int(H*0.08))
        for x in range(W):
            y = int(mid + amp*math.sin(2*math.pi*spin_hz*x/W))
            if 0 <= y < H: md[y, x] = 1.0
        md += 0.02*rng.standard_normal((H, W)).astype(np.float32)
    elif missile_type == "booster":
        md += 0.2*rng.random((H, W)).astype(np.float32)
    elif missile_type == "supersonic_cruise":
        if flight_state in ("maneuver","evasive","terminal"):
            # intermittent streaks indicating control activity
            for x in range(0, W, 8):
                md[:, x:x+2] = 0.5
            md += 0.05*rng.standard_normal((H, W)).astype(np.float32)
        else:
            # steady cruise: faint vibrations
            for x in range(W):
                y = int(mid + 2*math.sin(2*math.pi*50*x/W))
                if 0 <= y < H: md[y, x] = 0.3
            md += 0.02*rng.standard_normal((H, W)).astype(np.float32)
    else:
        md += 0.03*rng.standard_normal((H, W)).astype(np.float32)
    return np.clip(md, 0.0, 1.0)

# === Existing geometric measurement ===

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

# === Main measurement API (backwards compatible) ===

def simulate_radar_measurements(
    p_ecef: Tuple[float, float, float],
    v_ecef: Tuple[float, float, float],
    t_now: float,
    last_times: Dict[str, float],
    rng: np.random.Generator,
    meas_context: Optional[Dict] = None,
) -> List[Dict]:
    """
    Generate noisy radar measurements for all radars that emit at t_now.

    Inputs:
      - p_ecef, v_ecef: target position/velocity in ECEF (m, m/s)
      - t_now: current sim time (s)
      - last_times: dict radar_id -> last emission time (s), updated in-place
      - rng: numpy random Generator for reproducibility
      - meas_context: optional dict to drive advanced fields, e.g.:
          {
            "missile_type": "supersonic_cruise" | "ballistic" | "hgv" | "booster" | "warhead",
            "flight_state": "boost" | "cruise" | "maneuver" | "terminal" | "evasive",
            "spin_hz": 0.0,
            "radar_power_db": 100.0
          }

    Output list items (extends previous schema with advanced fields):
      {
        "id": radar_id,
        "z": np.array([range, az, el, dop]),
        "R_diag": [σ_r², σ_az², σ_el², σ_dop²],
        "true": (R_true, az_true, el_true, dop_true),
        "rcs_dbsm": float,
        "snr_db": float,
        "quality": float in [0,1],
        "hrr_profile": List[float],   # 1D array serialized
        "micro_doppler": List[List[float]]  # 2D array serialized
      }
    """
    meas_context = meas_context or {}
    missile_type = meas_context.get("missile_type", "supersonic_cruise")
    flight_state = meas_context.get("flight_state", "cruise")
    spin_hz = float(meas_context.get("spin_hz", 0.0))
    default_power_db = float(meas_context.get("radar_power_db", 100.0))

    outputs: List[Dict] = []
    p_ecef = np.array(p_ecef, dtype=float)
    v_ecef = np.array(v_ecef, dtype=float)
    alt_m = float(np.linalg.norm(p_ecef) - 6371000.0)
    mach = float(np.linalg.norm(v_ecef) / 343.0)

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

        # Advanced fields
        aspect = _aspect_angle_rad(p_ecef, v_ecef, np.array(radar["ecef_m"]))
        rcs_db = _simulate_rcs_dbsm(missile_type, aspect, alt_m, flight_state, rng)
        snr_db = _simulate_snr_db(R_true, rcs_db,
                                  radar_power_db=float(radar.get("power_db", default_power_db)))
        # Only synthesize heavier signatures when SNR sufficient
        if snr_db > 10.0:
            hrr = _simulate_hrr_profile(missile_type, aspect, mach, flight_state, n_bins=128, rng=rng)
            md  = _simulate_micro_doppler(missile_type, flight_state, spin_hz, H=64, W=64, rng=rng)
        else:
            hrr = (0.1 * rng.standard_normal(128)).astype(np.float32)
            md  = (0.1 * rng.standard_normal((64,64))).astype(np.float32)
            hrr = np.clip(hrr, 0.0, None)
            md  = np.clip(md,  0.0, 1.0)

        # Quality score from SNR (0..1)
        quality = float(min(max((snr_db - 5.0) / 20.0, 0.0), 1.0))

        outputs.append({
            "id": rid,
            "z": np.array([z_range, z_az, z_el, z_dop], dtype=float),
            "R_diag": radar["R_diag"],
            "true": (R_true, az_true, el_true, dop_true),
            "rcs_dbsm": float(rcs_db),
            "snr_db": float(snr_db),
            "quality": quality,
            "hrr_profile": hrr.tolist(),
            "micro_doppler": md.tolist(),
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
    # Generate measurements for this time (with a default supersonic context)
    meas = simulate_radar_measurements(
        p_ecef, v_ecef, t_now, last_times, rng,
        meas_context={"missile_type":"supersonic_cruise","flight_state":"cruise"}
    )
    # Example: print a quick summary
    for m in meas:
        rid = m["id"]
        z = m["z"]
        print(f"{rid}: range={z[0]:.1f} m, az={z[1]:.4f} rad, el={z[2]:.4f} rad, dop={z[3]:.2f} m/s | "
              f"RCS={m['rcs_dbsm']:.1f} dBsm, SNR={m['snr_db']:.1f} dB, Q={m['quality']:.2f}")

if __name__ == "__main__":
    main_entry()
