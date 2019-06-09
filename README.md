# md2slf

A workflow to predict and simulate 15N-1H OS-ssNMR SLF spectra from MD trajectories.


## Overview

A method and assortment of scripts to generate SLF spectra in SPARKY format from an MD simulation trajectory. A brief overview of the process:

1. Simulate a membrane protein of interest with any simulation package.

2. Load trajectory into VMD and run the [slf.tcl](slf.tcl) script to compute 15N chemical shifts and 15N-1H dipolar couplings.

3. Convert the output of [slf.tcl](slf.tcl) to SPARKY peak list format using the [slf2sparky.py](slf2sparky.py) script.

4. Simulate a 15N-1H SLF spectrum in SPARKY format from the peak list file using the [simulate_slf.py](simulate_slf.py) script with desired line widths and field strength (requires [NMRGlue](https://www.nmrglue.com/)).

The simulated spectrum can then be loaded into [SPARKY](https://nmrfam.wisc.edu/nmrfam-sparky-distribution/) with the corresponding peak list file. Overlays with experimental spectra can  be used to aid assignments.

## Quickstart guide

1. Copy the [slf.tcl](slf.tcl) into working directory and load simulation trajectory into VMD. I.e.,  from UNIX terminal:

		vmd simulation.psf simulation.dcd

2. Open "Extensions -> Tk Console"

3. Select all protein residues. Proline and terminal residues will be filtered out by the script. In the Tk Console:

		source slf.tcl
		set s [atomselect top "all protein"]
		slf $s -start 2500 -out slf_50-200ns

	This command will measure from frame 2500 until the end of the 	trajectory (10000). In this example, each frame equates to 0.02 ns. The first 50 ns of non-equilibrated trajectory is discarded from the analysis. A file named slf_50-200ns.dat will be created and contains five columns: 1. residue three-letter name, 2. residue number,  3. 15N chemical shift, 4. dipolar couping and 5. NH order parameter with respect to Z-axis used to scale down values.

	A more advanced example is:

		slf $s -start 2500 -step 2 -stop 7500 -flip 90.0 -ord 0.8 -out slf_50-200ns_unflipped

	This will calculate values from frame 2500 to 7500, skipping every second frame. Value will be tranformed with membrane normal perpendicular to the Z-axis (i.e., an unflipped bicelle). A global order parameter of 0.8 is also applied, which could account for fluctuation in the experimental alignment axis and/or rigid body fluctuation about the helical tilt axis.
		
4. Convert the slf.tcl output to SPARKY peak list format using the [slf2sparky.py](slf2sparky.py) script:

		python slf2sparky.py -i slf_50-200ns.dat -o slf_50-200ns.list -n

	Note that the "-n" flag is used to convert dipolar couplings to negative values, which are more natural to visualize in SPARKY.
	
5. Simulate an SLF spectrum from the peak list file using the [simulate_slf.py](simulate_slf.py) script. This requires the [NMRGlue](https://www.nmrglue.com/) library to be installed into Python 3.X.

		python simulate_slf.py \
       		-i slf_50-200ns.dat.list \
       		-o slf_50-200ns.ucsf \
       		--sw_n 20000.0 \
       		--sw_nh 31250.0 \
       		--size_n 2048 \
       		--size_nh 1024 \
       		--lw_n 100.0 \
       		--lw_nh 100.0 \
       		--freq_n 60.7639142 \
       		--carr_n 5600.0

	These are the default settings of the script, which specify a 15N spectral width of 20 kHz with 2048 points; a 15N-1H dipolar coupling spectral width of 31.25 kHz with 1024 points; linewidths of 100 Hz in both dimensions; 15N frequency of 60.7639142 MHz (i.e., from a 600 MHz spectrometer) and the carrier offset frequency set at 5600 Hz.

6. Load the spectrum into SPARKY:

		sparky slf_50-200ns.ucsf

7. Load the peak list file (slf_50-200ns.list; step 4) using the "rp" command.

## Development
 
 The following are on the to-do list:

* Support for predicting side chain peaks.
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



