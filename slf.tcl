# SLF-VMD
# -------
# Calculates a 15N-1H Separated Local Field (SLF) spectrum from a molecular dynamics simulation or PDB file.
# Chemical shifts are computed by projecting the 15N shift tensor onto the peptide plane. The average orientation
# of the shift tensor components over the trajectory is then used to calculate the chemical shift. 15N-1H dipolar
# couplings are found the same way, but using the average orientation of the N-H bond vector. Both shifts and dipolar
# coupling are scaled by a backbone order parameter, a user-defined order parameter and a user-defined flip angle.
#
# Written by: Daniel K. Weber
# Affiliation: University of Minnesota, MN, USA
# Funding sources: NIH R01GM064742 and NIH R01HL144130 (Gianluigi Veglia). AHA 19POST34420009 (DW).
#
# Original date: Apr 21 2019
# Last revision: Apr 22 2020 (no output results without scaling by order parameters)
#
# Basic usage
# -----------
# Load trajectory and copy slf.tcl into working directory. In Tk Console:
# > source slf.tcl
# > set p [atomselect top "all protein"]
# > slf $p
# This will run the script with default settings (see below).
#
# Full usage
# ----------
# As above but using optional flags.
# > slf $p -mol 1 -start 1000 -stop 5000 -step 10 -out mysystem.dat -ord 0.9 -flip 90.0
# Will run for molecule 1. Starts at frame 1000 and ends frame 5000, for every 10th frame.
# Outputs file mysystem.dat. An additional order parameter of 0.9 will be applied.
# Chemical shifts and dipolar couplings will transformed as if the system were aligned
# perpendicular to the magnetic field (i.e., an unflipped bicelle).

proc slf { selection { args } } {
    # Main function.
    #
    # Parameters
    # ----------
    # selection: obj
    #    VMD selection. The code below will automatically extract relevant info.
    #
    # Optional flags
    # --------------
    # -mol <molecule number>
    #    <molecule number>: "top" or int
    #       Working mol number.  Default: "top".
    # -start <frame>
    #    <frame>: int
    #       Starting frame. Default: 0.
    # -stop <frame>
    #    <frame>: int
    #       End frame. Default: [molinfo $mol get numframes] (last frame).
    # -step <frames>
    #    <frames>: int
    #       Compute only for every <frames>th frame. Default: 1.
    # -out <filename>
    #    <filename>: str
    #       Output file prefix. Default: slf.
    # -draw <switch>.
    #    <switch>: 0 (off) or 1 (on)
    #       Draw 1H-15N bond vector and 15N chemical shift PAS. Not recommended for trajectory. Default: 0.
    # -flip <angle>
    #    <angle>: float
    #       Flip angle in degrees. I.e., use 90.0 for unflipped bicelle. Default: 0.0.
    # -ord <order parameter>
    #    <order parameter>: float
    #       Addtional order parameter to apply on chemical shifts and dipolar couplings. Default: 1.0.
    # -pas <dxx> <dyy> <dzz>
    #    <dxx>: float
    #       Smallest component of 15N chemical shift PAS in PPM. Default: 57.3.
    #    <dyy>: float
    #       Middle component. Default: 81.2.
    #    <dzz>: float
    #       Largest component. Default: 228.1.
    # -dc <coupling constant>
    #    <coupling constant>: float
    #       Dipolar coupling constant in kHz. Default: 10.735.
    #
    # Outputs
    # -------
    # <prefix>: ASCII data file
    #    Residue name, residue number, 15N chemical shift, 15N-1H dipolar coupling and order parameter.

    # Handle args
    set mol [argparse $args "-mol" 1 top]
    set start [argparse $args "-start" 1 0]
    set stop [argparse $args "-stop" 1 [molinfo $mol get numframes]]
    set step [argparse $args "-step" 1 1]
    set out [argparse $args "-out" 1 "slf.dat"]
    set drawflag [argparse $args "-draw" 1 "null"]
    set flip [argparse $args "-flip" 1 0.0]
    set ord [argparse $args "-ord" 1 1.0]
    set dxx [argparse $args "-pas" 1 57.3]
    set dyy [argparse $args "-pas" 2 81.2]
    set dzz [argparse $args "-pas" 3 228.1]
    set beta [argparse $args "-beta" 1 17.0]
    set b [argparse $args "-dc" 1 10.735]

    set beta [expr $beta*-1]
    set B0 {0 0 1}
    set iso [expr {($dzz+$dyy+$dxx)/3}]
    set flip [expr {$flip*0.017453}]
    set outf [open $out w]

    # Filter selection for protein residues
    puts "Finding amides..."
    set resids [lsort -unique [$selection get resid]]
    set indices [$selection get index]
    # Filter out non-protein and N-terminal residues
    set presids {}
    foreach resid [lsort -integer $resids] {
	set temp1 [atomselect top "protein and index $indices and resid $resid and name N"]
	set nindex [$temp1 get index]
	set bonded [lindex [$temp1 getbonds] 0]
	set temp2 [atomselect top "index $nindex $bonded and name H N HN C CA"]
	set atoms($resid) [$temp2 get index]
	set num [llength $atoms($resid)]
	if { $num == 4 } {
	    lappend presids $resid
	}
	#puts "$resid: $atoms($resid)"
	$temp1 delete
	$temp2 delete
	unset temp1 temp2
    }
    puts "Found [llength $presids] residues!"
    puts "#resname\tresnum\tCS\tDC\tS\tCS_f\tDC_f"
    
    # Loop through each residue
    foreach resid $presids {
	set N [atomselect top "index $atoms($resid) and name N"]
 	set H [atomselect top "index $atoms($resid) and name HN H"]
	set C [atomselect top "index $atoms($resid) and name C"]
	set CA [atomselect top "index $atoms($resid) and name CA"]
	set resname [$N get resname]
	set mu_list {}
	set xx_list {}
	set yy_list {}
	set zz_list {}

	# Loop through every frame
	for {set i $start} {$i <= $stop} {incr i $step} {
	    $N frame $i
	    $H frame $i
	    $C frame $i
	    $CA frame $i

	    # Get coordinates of relavent atoms
	    set Ncoor [lindex [$N get "x y z"] 0]
	    set Hcoor [lindex [$H get "x y z"] 0]
	    set Ccoor [lindex [$C get "x y z"] 0]
	    set CAcoor [lindex [$CA get "x y z"] 0]

	    # Compute coodinates of 15N CS tensor components and dipolar coupling vector
	    set CACO [vecscale 0.5 [vecadd $Ccoor $CAcoor]]
	    set mu [vecsub $Hcoor $Ncoor]
	    set zz $mu
	    set yy [veccross [vecsub $CACO $Hcoor] [vecsub $CACO $CAcoor]]
	    set xx [veccross $zz $yy]

	    # Adjust PAS for beta rotation about syy
	    set rotation [ transabout $yy $beta deg ]
	    set xx [vectrans $rotation $xx]
	    set zz [vectrans $rotation $zz]

	    # Draw tensors (optional)
	    if { $drawflag != "null" } {
		draw_plane  $Ncoor [vecadd $Ncoor $yy] 0.05 0.75 gray
		draw_vector $Ncoor [vecadd $Ncoor $mu ] 1.5 white
		draw_vector $Ncoor [vecadd $Ncoor $yy] 1.5 green
		draw_vector $Ncoor [vecadd $Ncoor $xx] 1.5 red
		draw_vector $Ncoor [vecadd $Ncoor $zz] 1.5 blue	
	    }
	    
	    # Append lists for later averaging
	    lappend mu_list $mu
	    lappend xx_list $xx
	    lappend yy_list $yy
	    lappend zz_list $zz
	}
	
	# Compute normalised average coodinates and order parameter of 15N CS tensor components and dipolar coupling vectors for residue
	set xx_mean [vecnorm [vecaverage $xx_list]]
	set yy_mean [vecnorm [vecaverage $yy_list]]
	set zz_mean [vecnorm [vecaverage $zz_list]]
	set mu_mean [vecnorm [vecaverage $mu_list]]
	set mu_ord [vecorder $mu_list $mu_mean]
	
	# Compute chemical shift (apply flip and additional order parameter)
	set cszz [expr {$dzz*([vecdot $B0 $zz_mean])**2}]
	set csyy [expr {$dyy*([vecdot $B0 $yy_mean])**2}]
	set csxx [expr {$dxx*([vecdot $B0 $xx_mean])**2}]
	set cs1 [expr {$cszz+$csyy+$csxx}]
	set cs [format "%.2f" [expr {($cs1-$iso)*$ord*$mu_ord*0.5*(3*(cos($flip)**2)-1)+$iso}]]
	set csf [format "%.2f" [expr {($cs1-$iso)*$ord*0.5*(3*(cos($flip)**2)-1)+$iso}]]
	
	# Compute absolute dipolar coupling 
	set dc1 [expr {($b/2)*(3*(([vecdot $B0 $mu_mean])**2)-1)}]
	set dc [expr {$dc1*$ord*$mu_ord*0.5*(3*(cos($flip)**2)-1)}]
	set dc [format "%.2f" [expr abs($dc)]]
	set dcf [expr {$dc1*$ord*0.5*(3*(cos($flip)**2)-1)}]
	set dcf [format "%.2f" [expr abs($dcf)]]
	
	# Formatting
	set mu_ord [format "%.2f" $mu_ord]
	
	# Output to file
	puts "$resname\t$resid\t$cs\t$dc\t$mu_ord\t$csf\t$dcf"
	puts $outf "$resname\t$resid\t$cs\t$dc\t$mu_ord\t$csf\t$dcf"
	
	# Clean up
	$N delete
	$H delete
	$C delete
	$CA delete
	unset N H C CA
    }
    puts "CS and DC written to $out"
    close $outf
}


proc argparse { args flag input_index default } {
    # Crude argument parser that does the trick.
    #
    # Parameters
    # ----------
    # args: list of str and values
    #    List of arguments and values
    # flag: str
    #    <args> are searched for the presence of this flag
    # input_index: int
    #    Position that value occurs after <flag>. I.e., if "1", the
    #    value immediately following <flag> will be returned.
    # default: anything
    #    If <flag> isn't found in <args>, then this value is returned.
    #
    # Returns
    # -------
    # value: anything
    #    The value parsed in from <args> of default.
   
    set value $default
    if { [ lsearch $args $flag ] != -1 } { 
	set value [lindex $args [expr ([lsearch $args $flag ] + $input_index)]]
    } 
    return $value
}


proc draw_vector { coord1 coord2 length color } {
    # Draw arrow from <coord1> in the direction of <coord2> with <length>.
    #
    # Parameters
    # ----------
    # coord1: list of floats { x y z }
    #    Coordinate where arrow starts
    # coord2: list of floats { x y z }
    #    Coordinate in the direction arrow will point towards.
    # length: float or int
    #    Total length of arrow in angstroms.
    # color: str (must be a valid VMD color)
    #    Arrow color
    
    graphics top color $color
    graphics top material Opaque
    set scalar [expr $length / [veclength [vecsub $coord2 $coord1]]]
    set arrowtip [vecadd $coord1 [vecscale $scalar [vecsub $coord2 $coord1]]]
    set arrowbase [vecadd $coord1 [vecscale 0.8 [vecsub $arrowtip $coord1]]]
    graphics top cylinder $coord1 $arrowbase radius 0.1 resolution 100
    graphics top cone $arrowbase $arrowtip radius 0.2 resolution 100
}

proc draw_plane { coord1 coord2 width radius color } {
    graphics top color $color
    graphics top material Opaque
    set scalar [expr $width / [veclength [vecsub $coord2 $coord1]]]
    set point1 [vecadd $coord1 [vecscale $scalar [vecsub $coord2 $coord1]]]
    set point2 [vecadd $coord1 [vecscale $scalar [vecsub $coord1 $coord2]]]
    graphics top cylinder $point1 $point2 radius $radius resolution 100 filled yes
}


proc vecaverage { vectors } {
    # Compute average vector from list of vectors
    #
    # Parameters
    # ----------
    # vectors: list of lists of floats
    #    A list of vectors of any size
    #
    # Returns
    # -------
    # vecout: list of floats
    #    Average vector of list of vectors
    
    set numvectors [llength $vectors]
    set numelements [llength [lindex $vectors 0]]
    for { set i 0 } { $i < $numelements } { incr i } { set cumsum($i) 0.0 }

    # Iterate through each vector
    for { set i 0 } { $i < $numvectors } { incr i } {
	# Array cumsum element index
	for { set j 0 } { $j < $numelements } { incr j } {
	    set cumsum($j) [expr $cumsum($j) + [lindex [lindex $vectors $i] $j]]
	}
    }
    # Make outputlist
    for { set i 0 } { $i < $numelements } { incr i } { lappend vecout [expr $cumsum($i) / $numvectors] }
    return $vecout
}


proc vecorder { vectors vectors_average } {
    # Compute order parameter for list of vectors.
    #
    # Parameters
    # ----------
    # vectors: list of list of vectors
    #    A list of vectors. Three elements: {x y z}.
    # vectors_average: list
    #    The average of vectors list. Three elements {x y z}.
    #
    # Returns
    # -------
    # order: float
    #    Order parameter of vectors with respect to average vector (Z-axis).

    #set num [llength $vectors]

    # Get normalized average vector and align along X
    set vnorm_a [vecnorm $vectors_average]
    set vnorm_a [list [expr abs([lindex $vnorm_a 0])] [expr abs([lindex $vnorm_a 1])] [lindex $vnorm_a 2]]
    set vnorm_a_xy [vecnorm [list [lindex $vnorm_a 0] [lindex $vnorm_a 1]]]
    set rotation [transabout {0 0 1} [expr -acos([lindex $vnorm_a_xy 0])] rad ]
    set vnorm_a [vectrans $rotation $vnorm_a]
    
    set csum 0
    set num 0
    foreach vector $vectors {

	# Normalize vector and align along X
	# Skip bug that crashes analysis where vector is pefectly aligned with Z-axis. 
	if { [catch {
	    set vnorm_v [vecnorm $vector]
	    set vnorm_v [list [expr abs([lindex $vnorm_v 0])] [expr abs([lindex $vnorm_v 1])] [lindex $vnorm_v 2]]
	    set vnorm_v_xy [vecnorm [list [lindex $vnorm_v 0] [lindex $vnorm_v 1]]]
	    set rotation [transabout {0 0 1} [expr -acos([lindex $vnorm_v_xy 0])] rad ]
	    set vnorm_v [vectrans $rotation $vnorm_v]
	    set tmp [ expr {(3*([vecdot $vnorm_v $vnorm_a]**2)-1)/2}]
	    set csum [ expr {$csum + $tmp}]
	    incr num
	}
	     ]
	 } {
	    puts "WARNING! Vector $num will be skipped in order parameter calculation."
	}
    }
    return [expr {$csum/$num}]
}
