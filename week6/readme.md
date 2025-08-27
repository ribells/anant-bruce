Blender-Based Multi-Radar Missile Tracking, Prediction, and Interception
This project is an open, reproducible Blender/Python simulation that generates noisy multi-radar measurements, fuses them with an Extended Kalman Filter (EKF), predicts short-horizon target motion with a lightweight RNN, and guides an interceptor using proportional navigation (PN).

Goals
Simulate a network of ground radars producing range/azimuth/elevation/Doppler with controllable noise, update rates, and dropouts to emulate realistic tracking conditions.

Fuse asynchronous radar measurements into a continuous 3D state estimate using EKF with a nonlinear measurement model.

Learn a 2–5 s look-ahead predictor from EKF state histories using a small LSTM/GRU to improve near-term trajectory prediction under mild maneuvers.

Drive an interceptor with PN guidance informed by predicted target states; quantify miss distance improvements over PN on EKF-only predictions.

Provide scoped demos for MIRV split-up (multi-target association with JPDA-lite) and a maneuvering HGV track using adaptive filtering (e.g., UKF/IMM).

Why this matters
EKF remains a standard for radar-based tracking with nonlinear measurement models; combining multi-radar fusion improves convergence and robustness.

Short-horizon data-driven predictors can complement filters by capturing residual nonlinearities in near-term motion, especially during brief maneuvers.

PN is a widely used, interpretable guidance law that benefits from more accurate look-ahead state estimates for reduced miss distance.

JPDA-style association stabilizes multi-target tracking after MIRV separation, while adaptive filters (UKF/IMM) help maintain tracks on maneuvering HGV-like trajectories.

Features
Blender Earth scene with truth missile trajectory and radar nodes placed in ECEF; per-frame measurement synthesizer for each radar (range/az/el/Doppler).

EKF tracking in Cartesian state with radar measurement Jacobians and configurable process/measurement covariances.

Optional LSTM/GRU predictor that consumes the last K EKF states and outputs Δt-ahead state for 2–5 s look-ahead.

PN interceptor using LOS-rate and closing speed; batch Monte Carlo scripts to compute miss distance distributions.

MIRV MVP: bus separation into multiple RVs plus JPDA-lite association and simple interceptor-to-target assignment.

HGV MVP: glide-and-turn maneuver model with adaptive EKF/UKF/IMM comparison for track continuity under lateral maneuvers.

Repo structure
blender/: .blend scenes, assets, rendering helpers.

sim/: truth models, radar geometry, measurement generator, CSV logging.

filters/: ekf.py (EKF), ukf_or_imm.py (adaptive options).

fusion/: association.py (nearest-neighbor → JPDA-lite).

guidance/: pn.py (PN guidance and intercept logic).

data/: generated CSVs and Monte Carlo outputs.

paper/: outline, figures, and scripts to reproduce plots.

Getting started
Dependencies: Python 3.10+, Blender 4.x, numpy/scipy/matplotlib/pandas, and torch or tensorflow for the small RNN.

Run the Blender script to spawn radars and generate measurement CSVs; measurements include range/az/el/Doppler with Gaussian noise and optional dropouts.

Launch the EKF loop to fuse multi-radar data and log estimated states; verify RMSE convergence on baseline trajectories.

Train the short-horizon LSTM/GRU on EKF histories to predict 2–5 s ahead; integrate real-time inference to draw “ghost” future points.

Enable PN interceptor using predicted target states; run Monte Carlo scripts and plot miss-distance distributions.

Switch on MIRV mode for bus separation and JPDA-lite association; then try HGV mode with an adaptive filter for maneuvering targets.

Evaluation
Tracking: position/velocity RMSE vs. number of radars, SNR, and dropout rate; time-to-converge metrics.

Prediction: 2–5 s look-ahead MAE/RMSE for EKF-only vs. EKF+RNN on maneuver/noise sweeps.

Interception: PN miss distance/time-to-go using EKF-only vs. EKF+RNN predictions; distributions over 50–200 runs.

Multi-target: MIRV track continuity and identity maintenance metrics; assignment effectiveness vs. naive pairing.

Maneuvering: HGV track continuity and RMSE under lateral turns; adaptive filter vs. EKF-only comparison.

