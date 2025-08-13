# ---------------------------------------------------------------
#  NYC  →  DC   3-D  VISUALISER  (Blender 3.x)
#  globe-scale fix • 60-s flight • animated dotted predictor
#  Refactored to use a Ballistic class for the trajectory
# ---------------------------------------------------------------
import bpy
import math
import numpy as np

# -------------------- 1. PARAMETERS ----------------------------
FPS          = 24                 # scene frame-rate
FLIGHT_SECS  = 60                 # make it 1 minute
TOTAL_FRAMES = FPS * FLIGHT_SECS  # 1440 at 24 fps
SCALE        = 1/100              # 1 BU = 100 km  →  R⊕ ≈ 64 BU
DASH_EVERY   = 10                 # put a dot every N simulated points
DOT_RADIUS   = 0.6                # size of each “dash”

# -------------------- 2. CONSTANTS & HELPERS -------------------
R_EARTH = 6371e3                  # metres
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
    return np.array([x, y, z], dtype=float)

# -------------------- 3. CLEAN SCENE & BUILD GLOBE -------------
bpy.ops.object.select_all(action="SELECT")
bpy.ops.object.delete()

bpy.context.scene.frame_start = 1
bpy.context.scene.frame_end   = TOTAL_FRAMES

# Globe (scaled)
bpy.ops.mesh.primitive_uv_sphere_add(radius=R_EARTH*SCALE, segments=96, ring_count=48)
earth = bpy.context.object
earth.name = "Earth"
mat = bpy.data.materials.new("EarthMat")
mat.diffuse_color = (0.05,0.25,0.55,1)
earth.data.materials.append(mat)

# -------------------- 4. BALLISTIC CLASS -----------------------
class Ballistic:
    """
    Builds a great-circle chord plane from start→end on a sphere, then
    integrates a 2D ballistic path in that plane (x along the chord; z up
    along the start radial).
    """
    def __init__(self, start_ll, end_ll,
                 launch_angle_deg=45.0,
                 g=G, n_steps=N_STEPS, dt=DT, scale=SCALE,
                 name="Trajectory"):
        self.start_ll = start_ll
        self.end_ll   = end_ll
        self.theta    = math.radians(launch_angle_deg)
        self.g        = float(g)
        self.n_steps  = int(n_steps)
        self.dt       = float(dt)
        self.scale    = float(scale)
        self.name     = name

        # ECEF anchors
        self.v_start = ll_to_ecef(*self.start_ll)
        self.v_end   = ll_to_ecef(*self.end_ll)

        # Local basis: ex (along chord), ez (up/radial at start)
        chord   = self.v_end - self.v_start
        self.range_m = float(np.linalg.norm(chord))
        self.ex = chord / self.range_m
        self.ez = self.v_start / np.linalg.norm(self.v_start)  # "up" at start

        # Initial speed to span 'range_m' at angle theta (planar formula)
        self.v0 = math.sqrt(self.range_m*self.g / math.sin(2*self.theta))
        self.vx0, self.vz0 = self.v0*math.cos(self.theta), self.v0*math.sin(self.theta)

        self.states = None     # list of [x,z,vx,vz]
        self.pts    = None     # list of local points (scaled), origin at v_start

    def _rk4_step(self, state):
        # state = [x, z, vx, vz]; a_x = 0, a_z = -g (as in original)
        x, z, vx, vz = state
        k1 = np.array([vx, vz, 0.0, -self.g], dtype=float)
        k2 = np.array([vx, vz - 0.5*self.g*self.dt, 0.0, -self.g], dtype=float)
        k3 = k2.copy()
        k4 = np.array([vx, vz - self.g*self.dt, 0.0, -self.g], dtype=float)
        return state + self.dt*(k1 + 2*k2 + 2*k3 + k4)/6.0

    def simulate(self):
        """Integrate and build scaled local points (relative to v_start)."""
        s = np.array([0.0, 0.0, self.vx0, self.vz0], dtype=float)
        states = []
        for _ in range(self.n_steps):
            states.append(s.copy())
            s = self._rk4_step(s)
        self.states = states

        pts = []
        for s in states:
            x, z, _, _ = s
            world = self.v_start + self.ex*x + self.ez*z
            local = (world - self.v_start) * self.scale   # local coordinates, scaled
            pts.append(tuple(local))
        self.pts = pts
        return pts

    def make_curve(self, bevel_depth=0.4, material=None):
        """Create a 3D poly curve from simulated points."""
        if self.pts is None:
            self.simulate()
        curve_data = bpy.data.curves.new(self.name, type='CURVE')
        curve_data.dimensions = '3D'
        curve_data.bevel_depth = bevel_depth

        spline = curve_data.splines.new('POLY')
        spline.points.add(len(self.pts) - 1)
        for p, q in zip(spline.points, self.pts):
            p.co = (*q, 1.0)

        curve_obj = bpy.data.objects.new(self.name, curve_data)
        bpy.context.collection.objects.link(curve_obj)
        curve_obj.location = earth.location
        if material:
            curve_obj.data.materials.append(material)
        return curve_obj

# -------------------- 5. BUILD TRAJECTORY VIA CLASS ------------
B = Ballistic(NYC, DC, launch_angle_deg=LAUNCH_ANGLE_DEG, g=G,
              n_steps=N_STEPS, dt=DT, scale=SCALE, name="Traj")
pts = B.simulate()  # list of local (scaled) points
curve_obj = B.make_curve(bevel_depth=0.4, material=mat)

start_frame = 1

# -------------------- 6. MISSILE OBJECT ------------------------
class Missile:
    def __init__(self, radius, depth, location, material=None):
        # body
        bpy.ops.mesh.primitive_cylinder_add(radius=radius, depth=depth, location=location)
        missile = bpy.context.object
        if material:
            missile.data.materials.append(material)
        # simple nose cone
        nose_z = missile.location.z + depth/2 + 1.0
        bpy.ops.mesh.primitive_cone_add(radius1=radius*0.8, depth=2.0, location=(location[0], location[1], nose_z))
        cone = bpy.context.object
        if material:
            cone.data.materials.append(material)
        # parent cone to body
        cone.parent = missile

        missile.keyframe_insert("location", frame=start_frame)
        self.missile = missile

    def shoot(self, duration_seconds=0):
        """No-op placeholder; follow-path will drive motion."""
        pass

missile_1 = Missile(0.8, 5.0, pts[0], material=mat)

# Follow Path constraint + path animation
con = missile_1.missile.constraints.new('FOLLOW_PATH')
con.target = curve_obj
con.use_curve_follow = True
curve_data = curve_obj.data
curve_data.use_path = True
curve_data.path_duration = TOTAL_FRAMES
missile_1.missile.location = (0, 0, 0)
curve_data.eval_time = 0
curve_data.keyframe_insert(data_path="eval_time", frame=1)
curve_data.eval_time = TOTAL_FRAMES
curve_data.keyframe_insert(data_path="eval_time", frame=TOTAL_FRAMES)

# -------------------- 7. DOTTED PREDICTION LINE ----------------
for idx, p in enumerate(pts[::DASH_EVERY]):
    bpy.ops.mesh.primitive_uv_sphere_add(radius=DOT_RADIUS, location=p)
    dot = bpy.context.object
    dot.data.materials.append(mat)
    f_on  = int(idx*TOTAL_FRAMES/len(pts))
    dot.hide_viewport = dot.hide_render = True
    dot.keyframe_insert("hide_viewport", frame=f_on-1)
    dot.keyframe_insert("hide_render",   frame=f_on-1)
    dot.hide_viewport = dot.hide_render = False
    dot.keyframe_insert("hide_viewport", frame=f_on)
    dot.keyframe_insert("hide_render",   frame=f_on)

# -------------------- 8. CAMERA -------------------------------
cam_dist = R_EARTH*SCALE*2
bpy.ops.object.camera_add(location=(0, -cam_dist*2.2, cam_dist), rotation=(math.radians(60),0,0))
bpy.context.scene.camera = bpy.context.object

print(f"Initial speed ≈ {B.v0/1000:.2f} km/s  •  Flight time {FLIGHT_SECS} s")