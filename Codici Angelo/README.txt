README

Files:
1. mainHW4b_shooting.m - Solves optimal orbital transfer using indirect shooting method with Newton-Raphson
2. mainHW4b_fmincon.m - Solves same problem using MATLAB's built-in fmincon

Usage:
1. Run either file in MATLAB
2. Set save_pic = true to store automatically the results figures (default: false)

Both scripts:
- Compute time-optimal transfer from elliptical to circular orbit
- Generate plots of trajectory, states, and attitude

Key outputs::
- Transfer trajectory plots
- State variable histories and errors
- Quaternion and Bryant angle plots
- Thrust angle profiles
