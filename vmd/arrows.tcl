# ---------------------------------------------------------------------
# draw_arrow: straight arrow for forces
#   radius and length both scale with vector amplitude
# ---------------------------------------------------------------------
proc draw_arrow {molid start end base_radius} {
    set vec [vecsub $end $start]
    set amp [veclength $vec]
    if {$amp < 1e-6} return
    set r      [expr {$base_radius * $amp}]
    set middle [vecadd $start [vecscale 0.8 $vec]]
    graphics $molid cylinder $start $middle radius $r               resolution 20
    graphics $molid cone     $middle $end    radius [expr {$r*2.5}] resolution 20
}

# ---------------------------------------------------------------------
# draw_curved_arrow: 270-degree arc around bond axis for torques
#   center     : point the arc circles around (midpoint of bond)
#   axis       : unit vector of bond (inboard -> outboard)
#   signed_amp : FULL torque magnitude with sign from dot(T, axis);
#                positive = CCW viewed from outboard end (right-hand rule)
#   base_radius: shaft radius at amplitude = 1
#   arc_radius : radius of the circular arc in Å
# ---------------------------------------------------------------------
proc draw_curved_arrow {molid center axis signed_amp base_radius arc_radius} {
    if {[expr {abs($signed_amp)}] < 1e-6} return

    set sign [expr {$signed_amp >= 0 ? 1.0 : -1.0}]
    set amp  [expr {abs($signed_amp)}]
    set r    [expr {$base_radius * $amp}]
    set ar   [expr {$arc_radius  * $amp}]

    set ax [lindex $axis 0]
    if {[expr {abs($ax)}] < 0.9} { set ref {1.0 0.0 0.0} } \
    else                          { set ref {0.0 1.0 0.0} }
    set perp1 [vecnorm [vecsub $ref [vecscale [vecdot $ref $axis] $axis]]]
    set perp2 [vecscale $sign [veccross $axis $perp1]]

    set n_seg 24
    set sweep [expr {270.0 * 3.14159265358979 / 180.0}]

    set prev_pt ""
    for {set i 0} {$i <= $n_seg} {incr i} {
        set t   [expr {$sweep * double($i) / double($n_seg)}]
        set pt  [vecadd $center \
                    [vecadd [vecscale [expr {$ar*cos($t)}] $perp1] \
                            [vecscale [expr {$ar*sin($t)}] $perp2]]]
        if {$prev_pt ne ""} {
            graphics $molid cylinder $prev_pt $pt radius $r resolution 12
        }
        set prev_pt $pt
    }

    set cos_e  [expr {cos($sweep)}]
    set sin_e  [expr {sin($sweep)}]
    set pt_end [vecadd $center \
                    [vecadd [vecscale [expr {$ar*$cos_e}] $perp1] \
                            [vecscale [expr {$ar*$sin_e}] $perp2]]]
    set tangent [vecnorm \
                    [vecadd [vecscale [expr {-$sin_e}] $perp1] \
                            [vecscale   $cos_e          $perp2]]]
    set pt_tip  [vecadd $pt_end [vecscale [expr {$r * 5.0}] $tangent]]
    graphics $molid cone $pt_end $pt_tip radius [expr {$r * 2.5}] resolution 12
}

# ---------------------------------------------------------------------
# load_vector_data
#   Format: frame inboard_idx outboard_idx fx fy fz tx ty tz u uDot
# ---------------------------------------------------------------------
proc load_vector_data {filename} {
    global vector_data
    array unset vector_data
    set f [open $filename r]
    while {[gets $f line] >= 0} {
        set line [string trim $line]
        if {$line eq "" || [string match "#*" $line]} continue
        lassign $line frame_idx inboard_idx outboard_idx fx fy fz tx ty tz u uDot
        set vector_data($frame_idx,$inboard_idx) \
            [list $outboard_idx $fx $fy $fz $tx $ty $tz]
    }
    close $f
    puts "Loaded [array size vector_data] entries from $filename"
}

# ---------------------------------------------------------------------
# draw_vectors_for_frame
# ---------------------------------------------------------------------
proc draw_vectors_for_frame {frame_idx} {
    global vector_data force_scale torque_scale arrow_radius arc_radius
    set molid 0

    graphics $molid delete all

    set keys [array names vector_data "${frame_idx},*"]
    if {[llength $keys] == 0} {
        puts "Warning: no vector data for frame $frame_idx"
        return
    }

    foreach key $keys {
        set inboard_idx [lindex [split $key ","] 1]
        lassign $vector_data($key) outboard_idx fx fy fz tx ty tz

        # Positions of both atoms
        set sel_in  [atomselect $molid "index $inboard_idx"  frame $frame_idx]
        set sel_out [atomselect $molid "index $outboard_idx" frame $frame_idx]
        set pos_in  [lindex [$sel_in  get {x y z}] 0]
        set pos_out [lindex [$sel_out get {x y z}] 0]
        $sel_in  delete
        $sel_out delete

        # Bond axis — skip degenerate (same atom) entries
        set bond_vec [vecsub $pos_out $pos_in]
        if {[veclength $bond_vec] < 1e-6} {
            puts "Warning: inboard == outboard for index $inboard_idx, skipping torque"
            set bond_axis ""
        } else {
            set bond_axis [vecnorm $bond_vec]
        }
        set bond_mid [vecadd $pos_in [vecscale 0.5 $bond_vec]]

        # -- Force: straight arrow from inboard atom -------------------
        set f_raw [list $fx $fy $fz]
        set f_vec [vecscale $force_scale $f_raw]
        set f_end [vecadd $pos_in $f_vec]
        graphics $molid color red
        graphics $molid material Opaque
        draw_arrow $molid $pos_in $f_end $arrow_radius

        # -- Torque: curved arrow around bond axis ----------------------
        if {$bond_axis eq ""} continue
        set t_raw [list $tx $ty $tz]
        set t_mag [veclength $t_raw]                           ;# FULL magnitude for thickness

        if {$t_mag < 1e-6} continue
        # Sign from dot product: which way does T point along the bond?
        set t_sign [expr {[vecdot $t_raw $bond_axis] >= 0 ? 1.0 : -1.0}]
        # signed_amp carries direction; magnitude drives arc thickness and length
        set t_signed [expr {$t_sign * $t_mag * $torque_scale}]

        graphics $molid color blue
        graphics $molid material Opaque
        draw_curved_arrow $molid $bond_mid $bond_axis $t_signed $arrow_radius $arc_radius
    }
}

# ---------------------------------------------------------------------
# Frame-change callback
# ---------------------------------------------------------------------
proc on_frame_change {args} {
    draw_vectors_for_frame [molinfo 0 get frame]
}

# ---------------------------------------------------------------------
# Configuration
#   Your torque magnitudes are ~0.06, force magnitudes ~0.1–0.3.
#   Scale factors push those into visible Å-scale radii.
#   Tune arrow_radius first — it sets the maximum shaft thickness.
# ---------------------------------------------------------------------
set force_scale   5.0   ;# Å per unit force
set torque_scale  5.0   ;# multiplier on torque magnitude for shaft thickness
set arrow_radius  0.3   ;# max shaft radius in Å (at amplitude = 1 after scaling)
set arc_radius    1.0   ;# radius of the circular arc in Å

load_vector_data "vectors.dat"
trace add variable ::vmd_frame(0) write on_frame_change
draw_vectors_for_frame [molinfo 0 get frame]
