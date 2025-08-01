# https://www.w3schools.com/python/python_oop.asp
# https://www.youtube.com/watch?v=nmJqIaSZlRs

# ---------------------------------------------------------------
#  NYC  →  DC   3-D  VISUALISER  (Blender 3.x)
#  globe-scale fix • 60-s flight • animated dotted predictor
# ---------------------------------------------------------------
import bpy, math, numpy as np

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

# -------------------- 4. BALLISTIC TRAJECTORY ------------------
v_start = ll_to_ecef(*NYC)
v_end   = ll_to_ecef(*DC)
chord   = v_end - v_start
range_m = np.linalg.norm(chord)

theta = math.radians(LAUNCH_ANGLE_DEG)
v0    = math.sqrt(range_m*G / math.sin(2*theta))
vx0, vz0 = v0*math.cos(theta), v0*math.sin(theta)

ex = chord / range_m
ez = v_start / np.linalg.norm(v_start)
def rk4(state):
    x,z,vx,vz = state
    k1 = np.array([vx, vz, 0, -G])
    k2 = np.array([vx, vz - 0.5*G*DT, 0, -G])
    k3 = k2
    k4 = np.array([vx, vz - G*DT, 0, -G])
    return state + DT*(k1 + 2*k2 + 2*k3 + k4)/6

states=[]
s = np.array([0,0,vx0,vz0])
for _ in range(N_STEPS):
    states.append(s.copy())
    s = rk4(s)

# world points (scaled)
pts=[]
for s in states:
    x,z,_,_ = s
    pos = v_start + ex*x + ez*z
    pts.append((pos - v_start)*SCALE/1000.0)  # metres→km→BU

# -------------------- 5. TRAJECTORY CURVE ----------------------
curve_data = bpy.data.curves.new("Traj", type='CURVE')
curve_data.dimensions='3D'; curve_data.bevel_depth=0.4
spline = curve_data.splines.new('POLY')
spline.points.add(len(pts)-1)
for p,q in zip(spline.points, pts):
    p.co = (*q,1)

curve_obj = bpy.data.objects.new("TrajectoryCurve", curve_data)
bpy.context.collection.objects.link(curve_obj)
curve_obj.location = earth.location

start_frame = 1
# ------------------------- 6. MISSILES --------------------------
class Missile:
  def __init__(self, radius, depth, location):
    self.location = location
    
    #create the visual missile
    bpy.ops.mesh.primitive_cylinder_add(radius=radius, depth=depth)
    missile = bpy.context.object
    missile.location = location
    nose = missile.data.vertices[-1].co.z + 2.5
    bpy.ops.object.editmode_toggle() 
    bpy.ops.mesh.primitive_cone_add(radius1=0.8, depth=2, location=(0,0,nose)) 
    bpy.ops.object.editmode_toggle()
    missile.data.materials.append(mat)
    
    missile.keyframe_insert("location", frame=start_frame)
    
    self.missile = missile
    
  def shoot(self, duration):
	import time
	duration = duration  # Loop for the duration
	start_time = time.time()
	keyframe_time = 10
	
	print("Looping...")
	while time.time() - start_time < duration:
	  # Your code to be executed in the loop
	  self.location[0] += 1.0
	  self.missile.location = self.location
	  self.missile.keyframe_insert("location", frame=keyframe_time)
	  keyframe_time = keyframe_time + 10
	  time.sleep(1)  # Sleep one second while missile moves

missile_1 = Missile(0.8, 5, pts[0])
missile_1.shoot(10)

# follow path constraint
#con = missile.constraints.new('FOLLOW_PATH')
#con.target = curve_obj
#con.use_curve_follow = True
#curve_data.use_path=True
#curve_data.path_duration = TOTAL_FRAMES
#curve_data.eval_time = 0
#curve_data.keyframe_insert(data_path="eval_time", frame=1)
#curve_data.eval_time = TOTAL_FRAMES
#curve_data.keyframe_insert(data_path="eval_time", frame=TOTAL_FRAMES)

# -------------------- 7. DOTTED PREDICTION LINE ----------------
for idx,p in enumerate(pts[::DASH_EVERY]):
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

print(f"Initial speed ≈ {v0/1000:.2f} km/s  •  Flight time {FLIGHT_SECS} s")
