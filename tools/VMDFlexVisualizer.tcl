## VMD-compatible TCL script made to assist Robosample flex generation. 
## Example Usage inside VMD, after loading the desired system, with a molID of 0 (i.e):
##
## source Robosample_VMDFlexVisualizer.tcl ## This script
## colorFlex 0 system/all.cartesian.flex   ## Check a cartesian flex file
## colorFlex 0 system/all.pin.flex         ## Check a TD flex file
## colorFlex 0 system/rama.pin.flex
##

proc ::colorFlex {molID flexFileIn} {
    ## Load flex file
    set flexfile [open $flexFileIn "r"] 
    set lines [split [read $flexfile] "\n"]
    ## Make all bonds grey
    set numReps [molinfo $molID get numreps]
    for {set i 0} {$i < $numReps} {incr i 1} {
        mol delrep 0 $molID
    }
    mol selection "all"
    mol representation {Bonds 0.1}
    mol color {ColorID 2}
    mol addrep $molID 
    
    ## For each line in the flex file, generate a coloured representation
    mol representation {Bonds 0.3}
    set pinCounter 0
    set cartCounter 0
    set rootMobility "Weld"
    foreach line $lines {
        if {[lindex $line 2] == "Pin"} {
            mol color {ColorID 7}
            mol selection "index [lindex $line 0] [lindex $line 1]"
            mol addrep $molID
            incr pinCounter
        }
        if {[lindex $line 2] == "Cartesian"} {
            mol color {ColorID 11}
            mol selection " index [lindex $line 0] [lindex $line 1]"
            mol addrep $molID
            incr cartCounter
        }
        if {[lindex $line 0] == "-1"} {
	    set rootMobility [lindex $line 2]
	}
    }
    puts "Internal Pin joints: $pinCounter"
    puts "Internal Cartesian joints: $cartCounter"
    puts "Root mobility: $rootMobility"
}
