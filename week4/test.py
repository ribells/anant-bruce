# https://www.w3schools.com/python/python_oop.asp
# https://www.youtube.com/watch?v=nmJqIaSZlRs

# ---------------------------------------------------------------
#  NYC  →  DC   3-D  VISUALISER  (Blender 3.x)
#  globe-scale fix • 60-s flight • animated dotted predictor
# ---------------------------------------------------------------
import math

# -------------------- 1. PARAMETERS ----------------------------
FPS          = 24                 # scene frame-rate
FLIGHT_SECS  = 60                 # make it 1 minute
TOTAL_FRAMES = FPS * FLIGHT_SECS  # 1440 at 24 fps
SCALE        = 1/100              # 1 BU = 100 km  →  R⊕ ≈ 64 BU
DASH_EVERY   = 10                 # put a dot every N simulated points
DOT_RADIUS   = 0.6                # size of each “dash”

# -------------------- 2. CONSTANTS & HELPERS -------------------
R_EARTH = 6371000                 # metres
G       = 9.80665
NYC     = (40.730610, -73.935242)
DC      = (38.907192, -77.036871)
LAUNCH_ANGLE_DEG = 45.0
N_STEPS  = 600
DT       = FLIGHT_SECS / N_STEPS  # keep same duration whatever FPS

def ll_to_ecef(lat, lon, alt=0.0):
    lat, lon = math.radians(lat), math.radians(lon)
    x = (R_EARTH+alt)*math.cos(lat)*math.cos(lon)
    y = (R_EARTH+alt)*math.cos(lat)*math.sin(lon)
    z = (R_EARTH+alt)*math.sin(lat)
    return np.array([x, y, z])

# -------------------- 6. MISSILE OBJECT ------------------------
class Missile:
    def __init__(self, radius, depth, location):
        self.location = location
        self.vector = [0, 0, 0]
        self.acceleration = [0, 0, 0]

    def updatePosition(self, dt):
        # Update position of missile
        self.location[0] = self.location[0] + self.vector[0] * dt;
        self.location[1] = self.location[1] + self.vector[1] * dt;
        self.location[2] = self.location[2] + self.vector[2] * dt;

    def updateVelocity(self, dt):
        self.vector[0] += self.acceleration[0] * dt;
        self.vector[1] += self.acceleration[1] * dt;
        self.vector[2] += self.acceleration[2] * dt;

    def updateAcceleration(self, dt):
        ax = 0;
        ay = 0;
        az = 0;

        dx = self.location[0] * -1
        dy = self.location[1] * -1
        dz = self.location[2] * -1

        distSq = dx * dx + dy * dy + dz * dz;

        f = 0
        g = 0
        if (distSq != 0):
            g = G / (distSq * math.sqrt(distSq));
            f = 0.99

        ax += dx * f
        ay += dy * f
        az += dz * g

        self.acceleration[0] = ax
        self.acceleration[1] = ay
        self.acceleration[2] = az

    def shoot(self, duration, vector, velocity):
        import time
        self.vector = vector
        self.velocity = velocity
        start_time = time.time()

        print("Looping...")
        i = 1
        while ((time.time() - start_time < duration) and self.location[2] >= 0):
            # Your code to be executed in the loop
            self.updateAcceleration(.01)
            self.updateVelocity(.01)
            self.updatePosition(.01)

            print(i, self.location)
            i = i+1

        if(self.location[2] < 0):
            self.location[2] == 0

missile_1 = Missile(0.8, 5, [0,0,0])
missile_1.shoot(.1, [1, 0, 14], 100)  # shoot a missile and simulate for 100 timesteps

import numpy as np

def get_force_components(direction_vector, total_force):
    """
    direction_vector: tuple or list (dx, dy, dz)
    total_force: scalar magnitude of force
    
    Returns: (Fx, Fy, Fz)
    """
    d = np.array(direction_vector, dtype=float)
    magnitude = np.linalg.norm(d)
    
    if magnitude == 0:
        raise ValueError("Direction vector magnitude cannot be zero.")
    
    unit_vector = d / magnitude
    components = total_force * unit_vector
    return tuple(components)

# Example get_force_components usage
direction = (1, 0, 14)
F_total = 100  # Newtons

Fx, Fy, Fz = force_components(direction, F_total)
print(f"[{Fx:.2f}, {Fy:.2f}, {Fz:.2f}]")

def get_force_from_angles(F_total, azimuth_deg, elevation_deg):
    """
    F_total: total force magnitude
    azimuth_deg: angle in degrees from x-axis in the x-y plane (0° = along +x)
    elevation_deg: angle in degrees above the x-y plane (0° = flat, +90° = straight up)
    
    Returns: (Fx, Fy, Fz)
    """
    # Convert to radians
    azimuth = math.radians(azimuth_deg)
    elevation = math.radians(elevation_deg)
    
    Fx = F_total * math.cos(elevation) * math.cos(azimuth)
    Fy = F_total * math.cos(elevation) * math.sin(azimuth)
    Fz = F_total * math.sin(elevation)
    
    return Fx, Fy, Fz

# Example get_force_from_angles usage:
theta = 45     # azimuth in degrees
phi = 30       # elevation in degrees
F_total = 100  # Newtons

Fx, Fy, Fz = get_force_from_angles(F_total, theta, phi)
print(f"[{Fx:.2f}, {Fy:.2f}, {Fz:.2f}]")
