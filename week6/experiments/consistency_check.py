# experiments/consistency_check.py
import csv, numpy as np
from radars import RADARS_ECEF
from sim.truth_provider import ScriptedTruth
from sim.measurements import simulate_radar_measurements
from tracking.ekf import EKF

def run(duration_s=60.0, dt=0.5, qa=0.1,
        sigmas=(50.0, 0.0010, 0.0010, 3.0),
        out_csv="runs/run_consistency.csv"):
    rng = np.random.default_rng(42)
    last_times = {}
    anchor = RADARS_ECEF["ecef_m"]
    truth = ScriptedTruth(anchor)

    # EKF init near first truth
    p0, v0 = truth.step(0.0, dt)
    x0 = np.array([*p0, *v0], dtype=float)
    P0 = np.diag([1e4,1e4,1e4, 1e2,1e2,1e2])
    ekf = EKF(x0, P0, qa)

    R_diag = np.array([sigmas**2, sigmas**2, sigmas**2, sigmas**2])

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        header = ["t","px","py","pz","vx","vy","vz","nis","radar_id",
                  "range","az","el","dop","R11","R22","R33","R44"]
        w.writerow(header)

        t = 0.0
        while t <= duration_s:
            p, v = truth.step(t, dt)
            ekf.predict(dt)

            meas_list = simulate_radar_measurements(p, v, t, last_times, rng)
            for meas in meas_list:
                rid = meas["id"]
                z = np.array(meas["z"], dtype=float)
                R = np.diag(R_diag)
                radar_ecef = next(r["ecef_m"] for r in RADARS_ECEF if r["id"]==rid)
                nis = ekf.update_radar(z, R, radar_ecef)
                w.writerow([t, *p, *v, nis, rid, *z, *R_diag])

            t += dt
