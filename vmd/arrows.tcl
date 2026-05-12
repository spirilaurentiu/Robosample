# =====================================================================
# VMD VECTOR VISUALIZER
#
#   RED  = forces  (straight arrows)
#   BLUE = torques (curved arrows)
#
# ---------------------------------------------------------------------
# GLOBAL KNOBS
#
# FORCE ARROWS:
#   force_scale       -> scales the displayed arrow LENGTH
#                        (multiplies the raw force vector)
#   force_size_scale  -> uniform multiplier on ALL red dimensions
#                        (shaft radius, head length, head radius)
#                        does NOT affect the length
#
# TORQUE ARROWS:
#   torque_scale         -> controls torque arc size
#   torque_radius_scale  -> controls torque shaft thickness
#
# COMMON:
#   arrow_radius  -> master radius multiplier (affects both colours)
#   arc_radius    -> base curvature radius for torque arcs
#
# =====================================================================


# =====================================================================
# DRAW FORCE ARROW
#
# Proportions are intentionally identical to draw_curved_arrow:
#   shaft_r  = arrow_radius * force_size_scale * amp^0.55
#   head_len = shaft_r * 5.0   (clamped to 0.5 * amp to prevent flip)
#   head_r   = shaft_r * 2.0
#
# force_size_scale is the single knob that uniformly enlarges or
# shrinks the red arrow independently of the force it represents.
# =====================================================================
proc draw_arrow {molid start end} {

    global arrow_radius
    global force_size_scale

    set vec [vecsub $end $start]
    set amp [veclength $vec]

    if {$amp < 1e-6} return

    set dir [vecnorm $vec]

    # ------------------------------------------------------------------
    # SHAFT RADIUS  — identical formula to draw_curved_arrow
    # ------------------------------------------------------------------
    set shaft_r \
        [expr {$arrow_radius * $force_size_scale * pow($amp, 0.55)}]

    # ------------------------------------------------------------------
    # HEAD — same shaft-relative ratios as draw_curved_arrow
    # Clamp head_len to half the arrow so it never flips on short arrows
    # ------------------------------------------------------------------
    set head_len [expr {min($shaft_r * 5.0, $amp * 0.5)}]
    set head_r   [expr {$shaft_r * 2.0}]

    # ------------------------------------------------------------------
    # GEOMETRY
    # ------------------------------------------------------------------
    set middle [vecadd $end [vecscale -$head_len $dir]]

    graphics $molid cylinder \
        $start $middle \
        radius $shaft_r \
        resolution 20

    graphics $molid cone \
        $middle $end \
        radius $head_r \
        resolution 20
}


# =====================================================================
# DRAW TORQUE ARROW  (unchanged)
# =====================================================================
proc draw_curved_arrow {molid center axis signed_amp} {

    global arrow_radius
    global torque_radius_scale
    global arc_radius

    if {[expr {abs($signed_amp)}] < 1e-6} return

    set sign [expr {$signed_amp >= 0 ? 1.0 : -1.0}]
    set amp  [expr {abs($signed_amp)}]

    # ------------------------------------------------------------------
    # THICKNESS
    # ------------------------------------------------------------------
    set shaft_r \
        [expr {$arrow_radius * $torque_radius_scale * pow($amp, 0.55)}]

    # ------------------------------------------------------------------
    # ARC SIZE
    # ------------------------------------------------------------------
    set ar \
        [expr {$arc_radius * (0.7 + 0.6*$amp)}]

    # ------------------------------------------------------------------
    # PERPENDICULAR BASIS
    # ------------------------------------------------------------------
    set ax [lindex $axis 0]

    if {[expr {abs($ax)}] < 0.9} {
        set ref {1.0 0.0 0.0}
    } else {
        set ref {0.0 1.0 0.0}
    }

    set perp1 [vecnorm \
        [vecsub $ref [vecscale [vecdot $ref $axis] $axis]]]

    set perp2 [vecscale $sign [veccross $axis $perp1]]

    # ------------------------------------------------------------------
    # ARC
    # ------------------------------------------------------------------
    set n_seg 24
    set sweep [expr {270.0 * 3.14159265358979 / 180.0}]
    set prev_pt ""

    for {set i 0} {$i <= $n_seg} {incr i} {

        set t [expr {$sweep * double($i) / double($n_seg)}]

        set pt [vecadd $center \
            [vecadd \
                [vecscale [expr {$ar * cos($t)}] $perp1] \
                [vecscale [expr {$ar * sin($t)}] $perp2]]]

        if {$prev_pt ne ""} {
            graphics $molid cylinder \
                $prev_pt $pt \
                radius $shaft_r \
                resolution 12
        }

        set prev_pt $pt
    }

    # ------------------------------------------------------------------
    # ARROWHEAD
    # ------------------------------------------------------------------
    set cos_e [expr {cos($sweep)}]
    set sin_e [expr {sin($sweep)}]

    set pt_end [vecadd $center \
        [vecadd \
            [vecscale [expr {$ar * $cos_e}] $perp1] \
            [vecscale [expr {$ar * $sin_e}] $perp2]]]

    set tangent [vecnorm \
        [vecadd \
            [vecscale [expr {-$sin_e}] $perp1] \
            [vecscale  $cos_e          $perp2]]]

    set pt_tip [vecadd $pt_end [vecscale [expr {$shaft_r * 5.0}] $tangent]]

    graphics $molid cone \
        $pt_end $pt_tip \
        radius [expr {$shaft_r * 2.0}] \
        resolution 12
}


# =====================================================================
# LOAD VECTOR DATA
# =====================================================================
proc load_vector_data {filename} {

    global vector_data

    array unset vector_data

    set f [open $filename r]

    while {[gets $f line] >= 0} {

        set line [string trim $line]

        if {$line eq "" || [string match "#*" $line]} {
            continue
        }

        lassign $line \
            frame_idx \
            inboard_idx \
            outboard_idx \
            fx fy fz \
            tx ty tz \
            u uDot

        set vector_data($frame_idx,$inboard_idx) \
            [list \
                $outboard_idx \
                $fx $fy $fz \
                $tx $ty $tz]
    }

    close $f

    puts "Loaded [array size vector_data] entries"
}


# =====================================================================
# DRAW FRAME
# =====================================================================
proc draw_vectors_for_frame {frame_idx} {

    global vector_data
    global force_scale
    global torque_scale

    set molid 0

    graphics $molid delete all

    set keys [array names vector_data "${frame_idx},*"]

    if {[llength $keys] == 0} {
        puts "No vector data for frame $frame_idx"
        return
    }

    foreach key $keys {

        set inboard_idx [lindex [split $key ","] 1]

        lassign $vector_data($key) \
            outboard_idx \
            fx fy fz \
            tx ty tz

        # --------------------------------------------------------------
        # POSITIONS
        # --------------------------------------------------------------
        set sel_in  [atomselect $molid "index $inboard_idx"  frame $frame_idx]
        set sel_out [atomselect $molid "index $outboard_idx" frame $frame_idx]

        set pos_in  [lindex [$sel_in  get {x y z}] 0]
        set pos_out [lindex [$sel_out get {x y z}] 0]

        $sel_in  delete
        $sel_out delete

        # --------------------------------------------------------------
        # BOND AXIS
        # --------------------------------------------------------------
        set bond_vec [vecsub $pos_out $pos_in]

        if {[veclength $bond_vec] < 1e-6} continue

        set bond_axis [vecnorm $bond_vec]
        set bond_mid  [vecadd $pos_in [vecscale 0.5 $bond_vec]]

        # ==============================================================
        # FORCE  (red)
        # ==============================================================
        set f_vec [vecscale $force_scale [list $fx $fy $fz]]
        set f_end [vecadd $pos_in $f_vec]

        graphics $molid color red
        graphics $molid material Opaque

        draw_arrow $molid $pos_in $f_end

        # ==============================================================
        # TORQUE  (blue)
        # ==============================================================
        set t_raw [list $tx $ty $tz]
        set t_mag [veclength $t_raw]

        if {$t_mag < 1e-6} continue

        set t_sign   [expr {[vecdot $t_raw $bond_axis] >= 0 ? 1.0 : -1.0}]
        set t_signed [expr {$t_sign * $t_mag * $torque_scale}]

        graphics $molid color blue
        graphics $molid material Opaque

        draw_curved_arrow $molid $bond_mid $bond_axis $t_signed
    }
}


# =====================================================================
# FRAME CALLBACK
# =====================================================================
proc on_frame_change {args} {
    draw_vectors_for_frame [molinfo 0 get frame]
}


# =====================================================================
# GLOBAL CONFIGURATION
# =====================================================================

# ----------------------------------------------------------------------
# FORCE ARROWS (red)
# ----------------------------------------------------------------------

# Scales the LENGTH of red arrows (multiplies the raw force vector)
set force_scale 10.0

# Uniform multiplier on ALL red arrow dimensions
# (shaft radius, head length, head radius — but NOT arrow length).
# Increase to make red arrows fatter/larger; decrease for slicker look.
#
# TRY:
#   0.2 = very slick
#   0.5 = medium
#   1.0 = same weight as torque arrows at torque_radius_scale 1.0
#
set force_size_scale 0.5


# ----------------------------------------------------------------------
# TORQUE ARROWS (blue)
# ----------------------------------------------------------------------

# Scales the arc SIZE of blue arrows
set torque_scale 5.0

# Controls blue shaft thickness (and therefore head size)
set torque_radius_scale 1.0


# ----------------------------------------------------------------------
# COMMON
# ----------------------------------------------------------------------

# Master radius multiplier (applied to both red and blue)
set arrow_radius 1.0

# Base curvature radius for torque arcs
set arc_radius 1.0


# =====================================================================
# INITIALIZE
# =====================================================================

load_vector_data "vectors.dat"

trace add variable ::vmd_frame(0) write on_frame_change

draw_vectors_for_frame [molinfo 0 get frame]
