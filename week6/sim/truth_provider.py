# sim/truth_provider.py
from typing import Tuple, Protocol
Vec3 = Tuple[float, float, float]

class TruthProvider(Protocol):
    def reset(self) -> None: ...
    def step(self, t: float, dt: float) -> Tuple[Vec3, Vec3]: ...
    # returns (position_ecef_m, velocity_ecef_mps)

class ScriptedTruth:
    def __init__(self, anchor_ecef: Vec3):
        self.ax, self.ay, self.az = anchor_ecef
    def reset(self) -> None:
        pass
    def step(self, t: float, dt: float):
        px = self.ax + 50_000.0 + 300.0 * t
        py = self.ay + 10_000.0 + 50.0  * t
        pz = self.az + 1_000.0
        vx, vy, vz = 300.0, 50.0, 0.0
        return (px, py, pz), (vx, vy, vz)

class DynamicTruth:
    def __init__(self, state_fn):
        # state_fn(t) -> (position_ecef, velocity_ecef)
        self.state_fn = state_fn
    def reset(self) -> None:
        pass
    def step(self, t: float, dt: float):
        return self.state_fn(t)
