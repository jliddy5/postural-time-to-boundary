# Postural Time-to-Boundary (TtB)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17460159.svg)](https://doi.org/10.5281/zenodo.17460159)

A MATLAB toolbox for calculating time-to-boundary, a spatiotemporal measure of postural control that quantifies how close the center of pressure (CoP) or center of mass (CoM) trajectory is to the boundaries of the base of support.

## Overview

Time-to-Boundary (TtB) estimates the time until the extrapolated CoP/CoM trajectory would cross the boundaries of the base of support based on instantaneous position, velocity, and acceleration. This measure provides insight into postural stability and fall risk.

## Installation

1. Download or clone this repository
2. Add the folder to your MATLAB path:
   ```matlab
   addpath('path/to/postural-time-to-boundary')
   ```
3. Verify installation:
   ```matlab
   which ttb
   ```

## Quick Start

```matlab
% Load your data
load('example_cop.mat');  % Contains cop_x, cop_y, bound_pts

% Define parameters
fs = 120;                 % Sampling rate (Hz)
extrap_method = 2;      % TtB method (1 = Riccio, 2 = Slobounov)

% Compute time-to-boundary
[ttb, ttb_bound, bound_crossed, bound_percent] = ttb(cop_x, cop_y, 1/fs, bound_pts, extrap_method);

% Compute boundary-specific statistics
bound_labels = {'F', 'FR', 'BR', 'B', 'BL', 'FL'};
[mean_ttb, med_ttb, min_ttb] = ttbBoundary(ttb_bound, 50, 1, bound_labels);
```

## Methods

### Method 1: Riccio (Linear)
Uses position and velocity to extrapolate CoP/CoM trajectory:

```
A·τ + B = 0
```

where:
- `A = [v_y(t) - m·v_x(t)]`
- `B = [(p_y(t) - y_b) - m·(p_x(t) - x_b)]`
- `m` = boundary slope

**Reference**: Riccio, G. E. (1993). Information in movement variability about the qualitative dynamics of posture and orientation. In *Variability and Motor Control*.

### Method 2: Slobounov (Quadratic) [Default]
Incorporates acceleration for improved accuracy:

```
A·τ² + B·τ + C = 0
```

where:
- `A = [a_y(t) - m·a_x(t)] / 2`
- `B = [v_y(t) - m·v_x(t)]`
- `C = [(p_y(t) - y_b) - m·(p_x(t) - x_b)]`

**Reference**: Slobounov, S. M., Slobounova, E. S., & Newell, K. M. (1997). Virtual time-to-collision and human postural control. *Journal of Motor Behavior, 29*(3), 263-281.

## Functions

### `ttb(r_x, r_y, dt, bounds, extrap_method)`
Main function to compute time-to-boundary.

**Inputs:**
- `r_x` - ML position vector (numeric column vector)
- `r_y` - AP position vector (numeric column vector)
- `dt` - Time step (1/sampling_rate)
- `bounds` - Boundary coordinates matrix (n×2, clockwise order)
- `extrap_method` - Extrapolation method: 1 or 2 (default: 2)

**Outputs:**
- `ttb` - Minimum TtB time series (n_samples × 1)
- `ttb_bound` - TtB for each boundary (n_samples × n_boundaries)
- `bound_crossed` - Which boundary had minimum TtB at each time point
- `bound_percent` - Percentage of contacts to each boundary

### `ttbBoundary(ttb_bound, n_min, plot_flag, labels)`
Calculate mean, median, and minimum TtB for each boundary.

**Inputs:**
- `ttb_bound` - TtB matrix for each boundary (n_samples × n_boundaries)
- `n_min` - Number of minima to average (default: 10% of samples)
- `plot_flag` - Generate plots: 0 or 1 (default: 0)
- `labels` - Cell array of boundary labels (required if plotting)

**Outputs:**
- `mean_ttb` - Mean TtB for each boundary
- `med_ttb` - Median TtB for each boundary
- `min_ttb` - Minimum TtB estimate for each boundary

### `ttbBoundaryPercent(which_boundary, n_boundaries)`
Calculate percentage of virtual contacts to each boundary.

**Inputs:**
- `which_boundary` - Vector of crossed boundaries at each time point
- `n_boundaries` - Total number of boundaries

**Outputs:**
- `boundary_percent` - Percentage vector (n_boundaries × 1)

### `ttbMinN(ttb, n_min)`
Estimate minimum TtB by averaging the n_min lowest values (robust to noise).

**Inputs:**
- `ttb` - TtB time series vector
- `n_min` - Number of minimum values to average

**Outputs:**
- `min_ttb` - Average of n_min lowest values

## Example Analysis

See `test_ttb.m` for a complete example analyzing 30 seconds of quiet standing data:

```matlab
% The script demonstrates:
% - Loading CoP data
% - Defining base of support boundaries
% - Computing TtB with different methods
% - Generating plots and statistics
% - Fitting log-normal distributions
```

## Data Format

### Center of Pressure / Center of Mass
- **Column vectors** of ML and AP positions
- Units: typically millimeters
- Same sampling rate for both directions

### Base of Support Boundaries
- **n × 2 matrix**: [x_coordinates, y_coordinates]
- Points ordered **clockwise** (or counterclockwise, but consistently)
- Can start at any point (polygon is closed automatically)
- Minimum 3 points required

Example:
```matlab
bounds = [
    toe_R_x,    toe_R_y;     % Right toe
    m5_R_x,     m5_R_y;      % Right 5th metatarsal
    heel_R_x,   heel_R_y;    % Right heel
    heel_L_x,   heel_L_y;    % Left heel
    m5_L_x,     m5_L_y;      % Left 5th metatarsal
    toe_L_x,    toe_L_y      % Left toe
];
```

## Requirements

- **MATLAB R2019b or later** (requires `arguments` validation blocks)
- **Statistics and Machine Learning Toolbox** (for `lognfit` in example script)

## Citation

If you use this toolbox in your research, please cite:

```
Liddy, J. (2025). Postural Time-to-Boundary MATLAB Toolbox (Version 2.0.0) [Software]. 
https://doi.org/10.5281/zenodo.17460159
```

**BibTeX:**
```bibtex
@software{liddy2025ttb,
  author       = {Liddy, Josh},
  title        = {Postural Time-to-Boundary MATLAB Toolbox},
  year         = 2025,
  version      = {2.0.0},
  doi          = {10.5281/zenodo.17460159},
  url          = {https://github.com/jliddy5/postural-time-to-boundary}
}
```

## References

1. Riccio, G. E. (1993). Information in movement variability about the qualitative dynamics of posture and orientation. In *Variability and Motor Control* (pp. 317-357). Human Kinetics Publishers.

2. Slobounov, S. M., Slobounova, E. S., & Newell, K. M. (1997). Virtual time-to-collision and human postural control. *Journal of Motor Behavior, 29*(3), 263-281.

3. Hertel, J., & Olmsted-Kramer, L. C. (2007). Deficits in time-to-boundary measures of postural control with chronic ankle instability. *Gait & Posture, 25*(1), 33-39.

4. Van Wegen, E. E. H., Van Emmerik, R. E. A., Wagenaar, R. C., & Ellis, T. (2001). Stability boundaries and lateral postural control in parkinson's disease. *Motor Control, 5*(3), 254-269.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions, issues, and feature requests are welcome! Feel free to check the [issues page](https://github.com/jliddy5/postural-time-to-boundary/issues).

## Contact

Josh Liddy, jliddy@umass.edu

---

**Version**: 2.0.0  
**Last Updated**: October 2025
