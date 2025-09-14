# tracking/ekf.py
import numpy as np

def F_cv(dt: float):
    F = np.eye(6)
    F[0,3] = dt; F[1,4] = dt; F[2,5] = dt
    return F

def Q_cv(dt: float, qa: float):
    # qa: accel PSD (m^2/s^3) tuning param
    q = qa
    dt2, dt3, dt4 = dt*dt, dt*dt*dt, dt*dt*dt*dt
    Q = np.zeros((6,6))
    block = np.array([[dt4/4, dt3/2],
                      [dt3/2, dt2   ]]) * q
    Q[np.ix_([0,3],[0,3])] = block
    Q[np.ix_([1,4],[1,4])] = block
    Q[np.ix_([2,5],[2,5])] = block
    return Q

def wrap_angle(a):
    return (a + np.pi) % (2*np.pi) - np.pi

def h_and_H_radar(x, radar_ecef):
    # x = [px,py,pz,vx,vy,vz]
    px,py,pz,vx,vy,vz = x
    rx,ry,rz = radar_ecef
    dx,dy,dz = px-rx, py-ry, pz-rz
    r_xy2 = dx*dx + dy*dy
    r2 = r_xy2 + dz*dz
    r = np.sqrt(r2)
    az = np.arctan2(dy, dx)
    el = np.arctan2(dz, np.sqrt(r_xy2))
    dhat = np.array([dx,dy,dz]) / r
    dop = vx*dhat + vy*dhat + vz*dhat
    z = np.array([r, az, el, dop])

    # Jacobian H (6x4 transposed -> return 4x6)
    H = np.zeros((4,6))
    # range partials
    H[0,0] = dx/r; H[0,1] = dy/r; H[0,2] = dz/r
    # azimuth partials
    H[1,0] = -dy/r_xy2; H[1,1] = dx/r_xy2
    # elevation partials
    r_xy = np.sqrt(r_xy2) + 1e-12
    denom = r2
    H[2,0] = -(dx*dz)/(r_xy*denom)
    H[2,1] = -(dy*dz)/(r_xy*denom)
    H[2,2] = r_xy/denom
    # doppler partials
    H[3,3] = dhat; H[3,4] = dhat; H[3,5] = dhat
    return z, H

class EKF:
    def __init__(self, x0, P0, qa):
        self.x = x0.copy()
        self.P = P0.copy()
        self.qa = qa

    def predict(self, dt):
        F = F_cv(dt)
        Q = Q_cv(dt, self.qa)
        self.x = F @ self.x
        self.P = F @ self.P @ F.T + Q

    def update_radar(self, z_meas, R, radar_ecef):
        z_pred, H = h_and_H_radar(self.x, radar_ecef)
        y = z_meas - z_pred
        y = wrap_angle(y)
        y = wrap_angle(y)
        S = H @ self.P @ H.T + R
        K = self.P @ H.T @ np.linalg.inv(S)
        self.x = self.x + K @ y
        I = np.eye(6)
        self.P = (I - K @ H) @ self.P @ (I - K @ H).T + K @ R @ K.T
        nis = float(y.T @ np.linalg.inv(S) @ y)
        return nis
