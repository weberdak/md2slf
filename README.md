# md2slf

A workflow to predict and simulate 15N-1H OS-ssNMR SLF spectra from MD simulations.


## Overview

A method and assortment of scripts to generate SLF spectra in SPARKY format from an MD simulation trajectory. A brief overview of the process:

1. Simulate a membrane protein of interest with any simulation package.

2. Load trajectory into VMD and run the slf.tcl script to compute 15N chemical shifts and 15N-1H dipolar couplings.

3. Convert the output of slf.tcl to SPARKY peak list format using the slf2sparky.py script.

4. Simulate a 15N-1H SLF spectrum in SPARKY format from the peak list file using the simulate_slf.py script with desired line widths and field strength (requires NMRGlue).

The simulated spectrum can then be loaded into SPARKY with the corresponding peak list file. Overlays with experimental spectra can  be used to aid assignments.

## Quickguide

Load peak list using "rp" command.


## Development
 
 The following are on the to-do list:

* Support for predicting side chains peaks.
* Scripts to transform slf.tcl outputs with user-defined vaues of local order parameters (i.e., from relaxation data) and trialling alternative flip angles without having to repeat the entire analysis.

## Limitations

Users should be aware that the accuracy of the results are limited to:

#### 15N chemical shift tensor

Systematic errors exist in the calculation of 15N chemical shifts. The reason for this is not clear yet,  but could be related to the accuracy of the elements chosen for the principle axis system of the 15N chemical shift tensor.

#### Force field parameters
 
15N chemical shifts may also be affected by the force field parameters used for MD simulation. From my own experience, fitting experimental spectra to ideal helical models requires adjustment of the CA-N-H bond angle parameter to properly fit the dispersion of chemical shifts. This might be a reflection of the polarizability of the N-H bond with respect to the dielectric constant of the surrounding medium. I suspect that these limitations in accuracy might be worse for single-pass membrane proteins, which are more exposed to the low-dielectric alkane environment of the lipids and wouldn't be accounted for in common force fields.

#### Highly dynamic regions

Inaccuracy is very high for peaks associated with highly dynamic residues (low order parameters). This only be overcome by increasing the timescale of simulations or by using accelerated sampling methods.

#### Use of PDB structures

While the slf.tcl script does allow calculation directly from PDB structures (provided hydrogens are added), the results will generally be nonsense since very minor distortions in backbone geometry will produce large errors in peak positions. About 50 ns of simulation should be sufficient to average out these distortions for the most rigid regions.

#### Rigid body motions

A lipid bilayer will tend to ripple over the timescale of a simulation and cause protein to have slow-timescale fluctuations about the Z-axis. Since chemical shifts, dipolar coupling and order parameters are computer from vectors reference to the Z-axis, these motions will cause considerable errors in the measurement if not adequately samples (which they will generally not be). To constrain these motions, the following could be tried:

* Applying "tilt" collective variable restaints on helical residues (available in NAMD).
* Aligning protein coordinates to the first frame of the trajectory. This works well for larger membrane proteins the have been pre-oriented in lipid bilayers in the OPM Database.

However, these motions may have actual importance as a real mode of disorder - i.e., addtional degrees of freedom of single-pass mebrane proteins. If there is reason to suspect such disorder, the "-ord" option in the slf.tcl script can be adjusted to account for it. 



