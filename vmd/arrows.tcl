# =====================================================================
# VMD VECTOR VISUALIZER  --  per-body force reporter
#
#   RED  = force   (straight arrow)  drawn at the body's representative atom
#   BLUE = torque  (curved arrow)    curling about the torque axis (right hand)
#
# The CSV's force/torque is whatever terms the reporter summed into them
# (net = OpenMM applied + reaction, applied-only, or reaction-only -- see
# enable_reaction_reporter's include_openmm / include_reaction flags). This
# script just draws the two vectors it finds; it does not know which terms.
#
# ---------------------------------------------------------------------
# SIZING (two axes)
#   1. SYSTEM SIZE: lengths are FRACTIONS of the molecule's bounding-box diagonal.
#        force_frac  -> largest force arrow  = force_frac  * system_diag
#        torque_frac -> largest torque arc R = torque_frac * system_diag
#   2. MAGNITUDE: size ~ |vector|, normalized to `ref` via mag_mode
#        (linear | sqrt | log). ref = per-file max (auto), or set *_ref > 0 to a
#        fixed magnitude for cross-run comparability.
#   thickness_frac -> shaft radius as a fraction of the arrow's length
#   min_frac       -> skip arrows below this fraction of full size
# =====================================================================


# ---------------------------------------------------------------------
# GLOBAL CONFIGURATION  (tune these)
# ---------------------------------------------------------------------
set force_frac      0.30   ;# largest force arrow  = 30% of the system diagonal
set torque_frac     0.18   ;# largest torque arc R = 18% of the system diagonal
set mag_mode        sqrt   ;# linear | sqrt | log   (how size tracks |vector|)
set force_ref       0      ;# 0 = auto (per-file max) ; >0 = fixed reference magnitude
set torque_ref      0      ;# 0 = auto (per-file max) ; >0 = fixed reference magnitude
set thickness_frac  0.05   ;# shaft radius = 5% of the arrow's length
set min_frac        0.03   ;# skip arrows smaller than 3% of full size
set col_force       red    ;# force  (straight)
set col_torque      blue   ;# torque (curved)

set csv_file "3SN6.lig.0.reactions.csv"


# =====================================================================
# |v| -> fraction in [0,1], per mag_mode
# =====================================================================
proc mag_to_frac {mag ref} {
    global mag_mode
    if {$ref <= 0.0} { return 0.0 }
    set r [expr {$mag / double($ref)}]
    if {$r > 1.0} { set r 1.0 }
    switch $mag_mode {
        linear  { return $r }
        log     { return [expr {log(1.0 + 9.0 * $r) / log(10.0)}] }
        default { return [expr {sqrt($r)}] }
    }
}

# =====================================================================
# System size = bounding-box diagonal of molid (Angstrom), measured once
# =====================================================================
proc system_diag {molid} {
    set sel [atomselect $molid all]
    set mm [measure minmax $sel]
    $sel delete
    return [veclength [vecsub [lindex $mm 1] [lindex $mm 0]]]
}


# =====================================================================
# DRAW FORCE ARROW  (straight).  Thickness is a fraction of its length.
# =====================================================================
proc draw_arrow {molid start end} {
    global thickness_frac
    set vec [vecsub $end $start]
    set amp [veclength $vec]
    if {$amp < 1e-6} return
    set dir [vecnorm $vec]

    set shaft_r  [expr {$thickness_frac * $amp}]
    set head_len [expr {min($shaft_r * 5.0, $amp * 0.4)}]
    set head_r   [expr {$shaft_r * 2.2}]
    set middle   [vecadd $end [vecscale -$head_len $dir]]

    graphics $molid cylinder $start $middle radius $shaft_r resolution 20
    graphics $molid cone     $middle $end   radius $head_r resolution 20
}

# =====================================================================
# DRAW TORQUE ARROW  (curved arc of radius `ar` about `axis`, right-handed).
# =====================================================================
proc draw_curved_arrow {molid center axis ar} {
    global thickness_frac
    if {$ar < 1e-6} return
    set shaft_r [expr {$thickness_frac * $ar}]

    set ax [lindex $axis 0]
    if {[expr {abs($ax)}] < 0.9} {
        set ref {1.0 0.0 0.0}
    } else {
        set ref {0.0 1.0 0.0}
    }
    set perp1 [vecnorm [vecsub $ref [vecscale [vecdot $ref $axis] $axis]]]
    set perp2 [veccross $axis $perp1]

    set n_seg 24
    set sweep [expr {270.0 * 3.14159265358979 / 180.0}]
    set prev_pt ""
    for {set i 0} {$i <= $n_seg} {incr i} {
        set t  [expr {$sweep * double($i) / double($n_seg)}]
        set pt [vecadd $center \
            [vecadd [vecscale [expr {$ar * cos($t)}] $perp1] \
                    [vecscale [expr {$ar * sin($t)}] $perp2]]]
        if {$prev_pt ne ""} {
            graphics $molid cylinder $prev_pt $pt radius $shaft_r resolution 12
        }
        set prev_pt $pt
    }

    set cos_e [expr {cos($sweep)}]
    set sin_e [expr {sin($sweep)}]
    set pt_end [vecadd $center \
        [vecadd [vecscale [expr {$ar * $cos_e}] $perp1] \
                [vecscale [expr {$ar * $sin_e}] $perp2]]]
    set tangent [vecnorm \
        [vecadd [vecscale [expr {-$sin_e}] $perp1] \
                [vecscale $cos_e $perp2]]]
    set pt_tip [vecadd $pt_end [vecscale [expr {$shaft_r * 5.0}] $tangent]]
    graphics $molid cone $pt_end $pt_tip radius [expr {$shaft_r * 2.2}] resolution 12
}

# =====================================================================
# One straight force arrow: length proportional to |v|, normalized to `ref`.
# =====================================================================
proc draw_force_vec {molid pos vx vy vz ref maxlen color} {
    global min_frac
    set mag [veclength [list $vx $vy $vz]]
    if {$mag <= 1e-12} return
    set frac [mag_to_frac $mag $ref]
    if {$frac < $min_frac} return
    set end [vecadd $pos [vecscale [expr {$frac * $maxlen}] [vecnorm [list $vx $vy $vz]]]]
    graphics $molid color $color
    graphics $molid material Opaque
    draw_arrow $molid $pos $end
}

# =====================================================================
# One curved torque arc: radius proportional to |v|, normalized to `ref`.
# =====================================================================
proc draw_torque_vec {molid pos vx vy vz ref maxlen color} {
    global min_frac
    set mag [veclength [list $vx $vy $vz]]
    if {$mag <= 1e-12} return
    set frac [mag_to_frac $mag $ref]
    if {$frac < $min_frac} return
    graphics $molid color $color
    graphics $molid material Opaque
    draw_curved_arrow $molid $pos [vecnorm [list $vx $vy $vz]] [expr {$frac * $maxlen}]
}


# =====================================================================
# LOAD VECTOR DATA  (records the per-file max |force| and max |torque|)
# =====================================================================
proc load_vector_data {filename} {
    global vector_data max_f max_t
    array unset vector_data
    set max_f 0.0
    set max_t 0.0

    set f [open $filename r]
    while {[gets $f line] >= 0} {
        set line [string trim $line]
        if {$line eq "" || [string match "#*" $line]} { continue }

        # CSV: frame,replica,body_idx,atom_idx,fx,fy,fz,tx,ty,tz
        set fields [split $line ,]
        if {[llength $fields] != 10} {
            puts "arrows.tcl: skipping line with [llength $fields] fields (expected 10:\
 frame,replica,body_idx,atom_idx,fx,fy,fz,tx,ty,tz). OLD-format CSV? Rebuild the reporter."
            continue
        }
        lassign $fields frame_idx replica_idx body_idx atom_idx fx fy fz tx ty tz

        set vector_data($frame_idx,$atom_idx) [list $fx $fy $fz $tx $ty $tz]

        set fm [veclength [list $fx $fy $fz]] ; if {$fm > $max_f} { set max_f $fm }
        set tm [veclength [list $tx $ty $tz]] ; if {$tm > $max_t} { set max_t $tm }
    }
    close $f
    puts "Loaded [array size vector_data] entries; max|F|=$max_f max|T|=$max_t"
}


# =====================================================================
# DRAW FRAME
# =====================================================================
proc draw_vectors_for_frame {frame_idx} {
    global vector_data max_f max_t
    global force_ref torque_ref force_len torque_len col_force col_torque

    set molid 0
    graphics $molid delete all

    set keys [array names vector_data "${frame_idx},*"]
    if {[llength $keys] == 0} {
        puts "No vector data for frame $frame_idx"
        return
    }

    set f_ref [expr {$force_ref  > 0 ? $force_ref  : $max_f}]
    set t_ref [expr {$torque_ref > 0 ? $torque_ref : $max_t}]

    foreach key $keys {
        set atom_idx [lindex [split $key ","] 1]
        lassign $vector_data($key) fx fy fz tx ty tz

        set sel_atom [atomselect $molid "index $atom_idx" frame $frame_idx]
        set pos_atom [lindex [$sel_atom get {x y z}] 0]
        $sel_atom delete

        draw_force_vec  $molid $pos_atom $fx $fy $fz $f_ref $force_len  $col_force
        draw_torque_vec $molid $pos_atom $tx $ty $tz $t_ref $torque_len $col_torque
    }
}


# =====================================================================
# FRAME CALLBACK
# =====================================================================
proc on_frame_change {args} {
    draw_vectors_for_frame [molinfo 0 get frame]
}


# =====================================================================
# INITIALIZE
# =====================================================================
load_vector_data $csv_file

set _diag [system_diag 0]
set force_len  [expr {$force_frac  * $_diag}]
set torque_len [expr {$torque_frac * $_diag}]
puts "system diagonal = $_diag A ; max force arrow = $force_len A ; max torque arc R = $torque_len A"

trace add variable ::vmd_frame(0) write on_frame_change
draw_vectors_for_frame [molinfo 0 get frame]
