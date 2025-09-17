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
    c = float(np
