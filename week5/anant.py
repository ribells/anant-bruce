# Ready-to-run: textured Earth + missile animated along curved trajectory
# Paste into Blender Scripting editor and Run

import bpy, os, math, urllib.request, pathlib
import numpy as np
from mathutils import Vector, Matrix

# ---------------- USER PARAMETERS ----------------
LAUNCH_LAT, LAUNCH_LON = 40.730610, -73.935242   # NYC
TARGET_LAT, TARGET_LON = 38.907192, -77.036871   # DC

# visual / scale
BU_PER_KM = 1/100.0           # 1 BU = 100 km
BU_PER_M = BU_PER_KM / 1000.0 # 1e-5 BU per meter

# animation
FPS = 24
ANIM_SECONDS = 60
TOTAL_FRAMES = FPS * ANIM_SECONDS

# altitude profile (meters)
LAUNCH_ALT_M = 50.0
PEAK_ALT_M   = 40000.0   # 40 km peak (tweak up/down)
DESCEND_TO_SURFACE = True

# missile visual
MISSILE_SCALE_BU = 0.02   # Blender units (visual size)
DOT_EVERY = 8
DOT_RADIUS_BU = 0.02

# texture (download once)
TEX_URL = "https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57730/land_ocean_ice_2048.jpg"
TEX_DIR = pathlib.Path.home() / "blender_textures"
TEX_DIR.mkdir(parents=True, exist_ok=True)
TEX_PATH = TEX_DIR / "earth_texture.jpg"

# earth constants
R_EARTH_M = 6371e3
EARTH_RADIUS_BU = R_EARTH_M * BU_PER_M

# ---------------- helpers ----------------
def ll_to_ecef(lat_deg, lon_deg, alt_m=0.0):
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    r = R_EARTH_M + alt_m
    x = r * math.cos(lat) * math.cos(lon)
    y = r * math.cos(lat) * math.sin(lon)
    z = r * math.sin(lat)
    return np.array([x,y,z], dtype=float)

def great_circle_points(a, b, n_points):
    # a, b = ECEF vectors (meters)
    ra = np.linalg.norm(a); rb = np.linalg.norm(b)
    # normalize to sphere radius
    a_n = a / ra
    b_n = b / rb
    cosg = np.dot(a_n, b_n)
    cosg = max(-1.0, min(1.0, cosg))
    gamma = math.acos(cosg)
    pts = []
    if gamma == 0:
        return [a.copy()] * n_points
    for i in range(n_points):
        t = i/(n_points-1)
        # Slerp on sphere
        part1 = math.sin((1-t)*gamma) / math.sin(gamma)
        part2 = math.sin(t*gamma) / math.sin(gamma)
        p = part1 * a_n + part2 * b_n
        p = p * R_EARTH_M   # back to meters at surface
        pts.append(p)
    return pts

def altitude_profile(n_points, launch_alt, peak_alt, descent_to_surface=True):
    # simple smooth rise and fall: use a raised-cosine shape
    profile = []
    for i in range(n_points):
        t = i/(n_points-1)
        # ramp up to 0.5 then ramp down: use sin shaped envelope
        h = peak_alt * math.sin(math.pi * t)  # 0→peak→0
        # scale such that near edges it's at launch_alt
        # Blend launch altitude in
        h = launch_alt*(1 - 0.5*(1 - math.cos(math.pi*t))) + h*0.5
        if not descent_to_surface:
            # keep above launch altitude on descent
            h = max(h, launch_alt)
        profile.append(h)
    return profile

def world_to_bu(p_m):
    return tuple((p_m * BU_PER_M).tolist())

# ---------------- scene setup ----------------
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete()

scene = bpy.context.scene
scene.frame_start = 1
scene.frame_end = TOTAL_FRAMES
scene.render.fps = FPS

# ---------------- texture download (once) ----------------
if not TEX_PATH.exists():
    try:
        print("Downloading Earth texture to", TEX_PATH)
        urllib.request.urlretrieve(TEX_URL, str(TEX_PATH))
        print("Downloaded texture.")
    except Exception as e:
        print("Texture download failed; continuing without texture:", e)

# ---------------- create earth ----------------
bpy.ops.mesh.primitive_uv_sphere_add(radius=EARTH_RADIUS_BU, segments=128, ring_count=64)
earth = bpy.context.object
earth.name = "Earth"
bpy.ops.object.shade_smooth()

# Add a Sun light (directional light)
# The location defines the origin point, but the direction is determined by rotation.
bpy.ops.object.light_add(type='SUN', align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))

# Get a reference to the newly created light object
sun_light = bpy.context.object

# Set the rotation to define the light direction (e.g., pointing down and slightly forward)
# Rotation is in Euler angles (X, Y, Z) in radians
sun_light.rotation_euler = (0.785, -0.785, 0) # Approximately 45 degrees on X and -45 degrees on Y

# Adjust the sun light's properties (optional)
sun_light.data.energy = 2.0  # Set the strength of the light
sun_light.data.color = (1.0, 0.9, 0.8) # Set the color (RGB values from 0.0 to 1.0)

# material with texture
mat = bpy.data.materials.new("EarthMat")
mat.diffuse_color = (0.05,0.25,0.55,1)
mat.emissive_color = (0.05,0.25,0.55,1)
#turn on intensity here
mat.use_nodes = True
nodes = mat.node_tree.nodes
links = mat.node_tree.links
nodes.clear()
output = nodes.new(type='ShaderNodeOutputMaterial')
bsdf = nodes.new(type='ShaderNodeBsdfPrincipled')
bsdf.inputs['Roughness'].default_value = 0.8
links.new(bsdf.outputs['BSDF'], output.inputs['Surface'])

# image texture node
if TEX_PATH.exists():
    img_node = nodes.new(type='ShaderNodeTexImage')
    try:
        img_node.image = bpy.data.images.load(str(TEX_PATH))
        # Use UV coords from the sphere's default UV map
        tex_coord = nodes.new(type='ShaderNodeTexCoord')
        links.new(tex_coord.outputs['UV'], img_node.inputs['Vector'])
        links.new(img_node.outputs['Color'], bsdf.inputs['Base Color'])
    except Exception as e:
        print("Failed to load image texture:", e)
        bsdf.inputs['Base Color'].default_value = (0.05,0.25,0.55,1.0)
else:
    bsdf.inputs['Base Color'].default_value = (0.05,0.25,0.55,1.0)

earth.data.materials.append(mat)

# ensure UVs exist (UV sphere has UVs by default)
# but to be safe, set active UV map if none
if not earth.data.uv_layers:
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.uv.smart_project(angle_limit=66)
    bpy.ops.object.mode_set(mode='OBJECT')

# ---------------- compute trajectory (curved + elevated) ----------------
n_samples = max(64, TOTAL_FRAMES)  # number of path samples
p_launch = ll_to_ecef(LAUNCH_LAT, LAUNCH_LON, 0.0)
p_target = ll_to_ecef(TARGET_LAT, TARGET_LON, 0.0)

# great-circle points at surface (meters)
gc_pts = great_circle_points(p_launch, p_target, n_samples)

# altitude profile (meters)
alt_profile = altitude_profile(n_samples, LAUNCH_ALT_M, PEAK_ALT_M, DESCEND_TO_SURFACE)

# build world points (ECEF) with altitude added along radial direction
world_pts_m = []
for p_surf, h in zip(gc_pts, alt_profile):
    radial = np.array(p_surf) / np.linalg.norm(p_surf)  # unit radial
    world = p_surf + radial * h
    world_pts_m.append(world)

# ---------------- create curve and dotted predictor ----------------
pts_bu = [world_to_bu(p) for p in world_pts_m]

# curve
curve_data = bpy.data.curves.new("Trajectory", type='CURVE')
curve_data.dimensions = '3D'
curve_data.bevel_depth = 0.02
spline = curve_data.splines.new('POLY')
spline.points.add(len(pts_bu)-1)
for i, p in enumerate(pts_bu):
    x,y,z = p
    spline.points[i].co = (x,y,z,1.0)
curve_obj = bpy.data.objects.new("TrajectoryCurve", curve_data)
bpy.context.collection.objects.link(curve_obj)

# dotted predictor
dot_objs = []
for i,p in enumerate(pts_bu):
    if i % DOT_EVERY == 0:
        bpy.ops.mesh.primitive_uv_sphere_add(radius=DOT_RADIUS_BU, location=p)
        dot = bpy.context.object
        dot.data.materials.append(mat)
        dot_objs.append(dot)

# ---------------- missile mesh ----------------
# Create a small missile model and keep it local (we will keyframe it)
# Cylinder body aligned so +Y points forward
bpy.ops.mesh.primitive_cylinder_add(radius=0.02, depth=0.08, location=(0,0,0))
body = bpy.context.object
bpy.ops.mesh.primitive_cone_add(radius1=0.02, depth=0.04, location=(0, 0.06, 0))
nose = bpy.context.object
# join nose -> body
body.select_set(True)
nose.select_set(True)
bpy.context.view_layer.objects.active = body
bpy.ops.object.join()
missile = bpy.context.object
missile.name = "Missile"
missile.data.materials.append(mat)
bpy.ops.object.shade_smooth()

# ---------------- animate missile along sampled frames ----------------
frame_start = 1
frame_end = TOTAL_FRAMES
scene.frame_end = frame_end

# Sample trajectory to frames uniformly
for frame_idx in range(frame_start, frame_end+1):
    t = (frame_idx - frame_start) / (frame_end - frame_start)  # 0..1
    sample_idx = int(t * (len(world_pts_m)-1))
    pos_m = world_pts_m[sample_idx]
    pos_bu = Vector(pos_m * BU_PER_M)
    missile.location = pos_bu
    missile.keyframe_insert(data_path="location", frame=frame_idx)
    # orient missile to face along next tangent (use next point for direction)
    next_idx = min(sample_idx+1, len(world_pts_m)-1)
    dir_vec = Vector(world_pts_m[next_idx] - world_pts_m[sample_idx])
    # tangent direction projected to local tangent plane at point: remove radial component
    radial = Vector(world_pts_m[sample_idx]).normalized()
    dir_tangent = (dir_vec - radial * dir_vec.dot(radial))
    if dir_tangent.length > 1e-6:
        # align +Y to dir_tangent
        forward = dir_tangent.normalized()
        # build rotation matrix: +Y -> forward, pick +Z from radial
        right = forward.cross(radial).normalized()
        up = right.cross(forward).normalized()
        R = Matrix(((right.x, forward.x, up.x),
                    (right.y, forward.y, up.y),
                    (right.z, forward.z, up.z)))
        missile.rotation_mode = 'XYZ'
        missile.rotation_euler = R.to_euler('XYZ')
    missile.keyframe_insert(data_path="rotation_euler", frame=frame_idx)

# ---------------- camera ----------------
cam_dist = EARTH_RADIUS_BU * 2.3
bpy.ops.object.camera_add(location=(0, -cam_dist, cam_dist*0.5), rotation=(math.radians(60), 0, 0))
scene.camera = bpy.context.object

print("Done. Texture cached at:", TEX_PATH)
print("If you don't see the texture, switch viewport to Material Preview or Rendered shading.")
