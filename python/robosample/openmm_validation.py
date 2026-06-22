"""openmm_validation.py

Validate the C++ OpenMM build against a Python/OpenMM reference, and diagnose
common failure modes (infinite / NaN energy).

The C++ ``Context`` builds a single OpenMM ``System`` from the SoA
``system_topology`` (``Context.initialize_openmm``) and reports the potential
energy of the reference structure (``Context.calc_openmm_potential_energy``).
Here we build the *same* physical system directly with ParmEd + OpenMM from the
original ``prmtop``/``inpcrd`` and compare the two potential energies.

Because the potential energy is a sum over interactions, it is invariant to atom
ordering, so the BFS/compound ordering used on the C++ side and the prmtop
ordering used by ParmEd must agree to within floating-point/platform tolerance.

Per-force-group comparison
--------------------------
When the C++ side is initialized with separate force groups on
(``context.set_separate_force_groups(True)`` *before* ``initialize_openmm()``),
each OpenMM ``Force`` lands in its own force group and the C++ side can report
the potential energy of each group individually. :func:`compare_by_force_group`
asks for that breakdown, builds the matching per-group breakdown on the Python
reference side, and lines the two up by OpenMM force-class name. This makes it
obvious *which* term (nonbonded, bonds, angles, torsions, ...) is responsible
for a total-energy mismatch instead of only seeing the aggregate.

Notes
-----
* For an exact match, build both sides on the OpenMM ``Reference`` platform.
  Mixed-precision ``CUDA``/``OpenCL`` can differ by a small relative amount.
* The settings below (NoCutoff, no constraints, no implicit solvent, no CM
  motion removal) mirror the C++ defaults for a gas-phase single-point energy.
"""

from __future__ import annotations

import math
import os

import numpy as np
import openmm.unit as unit
import parmed as pmd

import openmm as mm
import robosample


def _build_reference_system(
    context: robosample.Context,
    prmtop_path: str | os.PathLike[str],
    inpcrd_path: str | os.PathLike[str],
) -> tuple[pmd.amber.AmberParm, mm.System]:
    """Build the ParmEd + OpenMM reference ``System`` (shared by the helpers).

    Kept in one place so the plain and per-group reference energies are
    guaranteed to come from an identically configured system. The nonbonded
    method, cutoff, implicit-solvent model, and (for Ewald/PME) the reciprocal-
    space tolerance are taken from ``context.system_topology`` so the reference
    mirrors whatever the C++ side built -- gas phase, GBSA, or explicit-solvent
    PME alike.
    """
    sys_top = context.system_topology
    implicit_solvent = mm.app.OBC2 if sys_top.use_gbsa_obc2 else None

    # Map the Robosample nonbonded method onto the OpenMM app constant. ParmEd's
    # createSystem reads the periodic box straight from `parm` (loaded with the
    # inpcrd), so PME/Ewald/CutoffPeriodic pick up the same box the C++ side uses.
    method_map = {
        robosample.NonbondedMethod.NoCutoff: mm.app.NoCutoff,
        robosample.NonbondedMethod.CutoffNonPeriodic: mm.app.CutoffNonPeriodic,
        robosample.NonbondedMethod.CutoffPeriodic: mm.app.CutoffPeriodic,
        robosample.NonbondedMethod.Ewald: mm.app.Ewald,
        robosample.NonbondedMethod.PME: mm.app.PME,
    }
    nb_method = method_map.get(sys_top.nonbonded_method, mm.app.NoCutoff)
    is_periodic = sys_top.nonbonded_method in (
        robosample.NonbondedMethod.CutoffPeriodic,
        robosample.NonbondedMethod.Ewald,
        robosample.NonbondedMethod.PME,
    )

    parm = pmd.load_file(str(prmtop_path), xyz=str(inpcrd_path))
    system = parm.createSystem(
        nonbondedMethod=nb_method,
        nonbondedCutoff=sys_top.nonbonded_cutoff,
        constraints=None,
        implicitSolvent=implicit_solvent,
        rigidWater=False,
        removeCMMotion=False,
    )

    # Match the C++ side's PME/Ewald reciprocal-space accuracy and dispersion
    # correction so the totals line up rather than differing by a tunable error.
    if is_periodic:
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setEwaldErrorTolerance(sys_top.ewald_error_tolerance)
                force.setUseDispersionCorrection(True)
    return parm, system


def reference_potential_energy(
    context: robosample.Context,
    prmtop_path: str | os.PathLike[str],
    inpcrd_path: str | os.PathLike[str],
    platform_name: str = "Reference",
) -> float:
    """Potential energy [kJ/mol] from a direct ParmEd + OpenMM build."""
    parm, system = _build_reference_system(context, prmtop_path, inpcrd_path)
    integrator = mm.VerletIntegrator(0.001)
    platform = mm.Platform.getPlatformByName(platform_name)
    mm_context = mm.Context(system, integrator, platform)
    mm_context.setPositions(parm.positions)
    mm_context.computeVirtualSites()  # place EPs from their frames (no-op if none)
    state = mm_context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def reference_potential_energy_by_group(
    context: robosample.Context,
    prmtop_path: str | os.PathLike[str],
    inpcrd_path: str | os.PathLike[str],
    platform_name: str = "Reference",
) -> tuple[float, dict[str, float]]:
    """
    Reference total + per-force-class potential energies [kJ/mol].

    Each OpenMM ``Force`` in the reference system is placed in its own force
    group, then queried individually. Energies are aggregated by force *class*
    name (``type(force).__name__``) so the result lines up with the C++ side,
    which tags its groups with the same class names. (Aggregating by class makes
    the mapping robust even when a class appears more than once or the two sides
    add forces in a different order.)

    Returns
    -------
    (total, by_class) : total PE and {force_class_name: summed PE}.
    """
    parm, system = _build_reference_system(context, prmtop_path, inpcrd_path)

    forces = list(system.getForces())
    # One group per force (0..n-1). OpenMM supports 32 groups (0-31).
    if len(forces) > 32:
        raise RuntimeError(
            f"reference system has {len(forces)} forces; only 32 force groups exist"
        )
    for i, force in enumerate(forces):
        force.setForceGroup(i)

    integrator = mm.VerletIntegrator(0.001)
    platform = mm.Platform.getPlatformByName(platform_name)
    mm_context = mm.Context(system, integrator, platform)
    mm_context.setPositions(parm.positions)
    mm_context.computeVirtualSites()  # place EPs from their frames (no-op if none)

    total = (
        mm_context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )

    by_class: dict[str, float] = {}
    for i, force in enumerate(forces):
        # getState's `groups` is a bitmask, matching the C++ side's `1 << group`.
        e = (
            mm_context.getState(getEnergy=True, groups=(1 << i))
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )
        name = type(force).__name__
        by_class[name] = by_class.get(name, 0.0) + e

    return total, by_class


def _normalize_cpp_groups(raw) -> dict[str, float]:
    """
    Aggregate the C++ per-group breakdown into {force_class_name: summed PE}.

    Accepts whatever ``calc_openmm_potential_energy_by_group`` hands back for the
    per-group list: each entry may be a ``(group, name, energy)`` tuple/list or
    an object exposing ``.name`` / ``.energy`` (the bound ``ForceGroupEnergy``).
    """
    by_class: dict[str, float] = {}
    for entry in raw:
        if hasattr(entry, "name") and hasattr(entry, "energy"):
            name = str(entry.name)
            energy = float(entry.energy)
        else:  # tuple/list: (group, name, energy)
            _group, name, energy = entry
            name = str(name)
            energy = float(energy)
        by_class[name] = by_class.get(name, 0.0) + energy
    return by_class


def diagnose_positions(context: robosample.Context, decimals: int = 6) -> list[str]:
    """
    Sanity-check the coordinates the C++ side will hand to OpenMM.

    Infinite energy almost always means two particles share a position (a 1/r
    term with r == 0). This reports NaNs, atoms piled at the origin, and exact
    coordinate duplicates -- the usual culprits.
    """
    st = context.system_topology
    n = int(st.num_atoms)
    x = np.asarray(list(st.atoms_x), dtype=float)
    y = np.asarray(list(st.atoms_y), dtype=float)
    z = np.asarray(list(st.atoms_z), dtype=float)

    issues: list[str] = []
    if not (len(x) == len(y) == len(z) == n):
        issues.append(
            f"length mismatch: num_atoms={n} but len(x,y,z)=({len(x)},{len(y)},{len(z)})"
        )
        print("positions: " + issues[-1])
        return issues

    pts = np.stack([x, y, z], axis=1)

    if np.isnan(pts).any():
        issues.append(
            f"{int(np.isnan(pts).any(axis=1).sum())} atom(s) have NaN coordinates"
        )

    zero_rows = int(np.all(pts == 0.0, axis=1).sum())
    if zero_rows > 0:
        issues.append(f"{zero_rows} atom(s) sit exactly at the origin (0,0,0)")

    # Exact coordinate duplicates -> r == 0 -> infinite nonbonded energy.
    _, counts = np.unique(np.round(pts, decimals), axis=0, return_counts=True)
    dup_sites = int((counts > 1).sum())
    if dup_sites > 0:
        dup_atoms = int(counts[counts > 1].sum())
        issues.append(
            f"{dup_atoms} atom(s) share {dup_sites} coordinate(s) "
            "(overlap -> infinite 1/r). Most likely multiple instances of a "
            "molecule type all placed at the prototype's coordinates."
        )

    print(
        f"positions: n={n}, "
        f"x[{x.min():.3f},{x.max():.3f}] "
        f"y[{y.min():.3f},{y.max():.3f}] "
        f"z[{z.min():.3f},{z.max():.3f}] nm"
    )
    if issues:
        print("PROBLEMS:")
        for msg in issues:
            print("  - " + msg)
    else:
        print("  positions look sane (no NaN, no origin pile-up, no exact duplicates)")
    return issues


def compare(
    context: robosample.Context,  # robosample.Context (already loaded via load_amber)
    prmtop_path: str | os.PathLike[str],
    inpcrd_path: str | os.PathLike[str],
    abs_tol_kj_per_mol: float = 1.0e-3,
    rel_tol: float = 1.0e-4,
    platform_name: str = "Reference",
    verbose: bool = True,
) -> tuple[bool, float, float]:
    """
    Build the C++ OpenMM system, compute both energies, and compare.

    Returns
    -------
    (ok, cpp_pe, ref_pe) : the pass/fail flag and both energies [kJ/mol].
    """
    if not context.initialize_openmm():
        raise RuntimeError("Context.initialize_openmm() failed")

    cpp_pe = context.calc_openmm_potential_energy()
    ref_pe = reference_potential_energy(
        context, prmtop_path, inpcrd_path, platform_name
    )

    if not math.isfinite(cpp_pe):
        print(f"C++ OpenMM energy is not finite ({cpp_pe}); diagnosing positions:")
        diagnose_positions(context)

    abs_diff = abs(cpp_pe - ref_pe)
    rel_diff = abs_diff / max(abs(ref_pe), 1.0e-12)
    ok = math.isfinite(cpp_pe) and (
        abs_diff <= abs_tol_kj_per_mol or rel_diff <= rel_tol
    )

    if verbose:
        print("OpenMM potential energy comparison [kJ/mol]")
        print(f"  C++ (robosample) : {cpp_pe:.6f}")
        print(f"  Python reference : {ref_pe:.6f}")
        print(f"  abs diff         : {abs_diff:.3e}")
        print(f"  rel diff         : {rel_diff:.3e}")
        print(f"  result           : {'PASS' if ok else 'FAIL'}")

    return ok, cpp_pe, ref_pe


def compare_by_force_group(
    context: robosample.Context,  # robosample.Context (already loaded via load_amber)
    prmtop_path: str | os.PathLike[str],
    inpcrd_path: str | os.PathLike[str],
    abs_tol_kj_per_mol: float = 1.0e-3,
    rel_tol: float = 1.0e-4,
    platform_name: str = "Reference",
    verbose: bool = True,
) -> tuple[bool, float, float, dict[str, tuple[float, float]]]:
    """
    Per-force-group energy comparison between the C++ build and the reference.

    Turns on separate force groups on the C++ side, initializes, and asks for the
    total PE plus the per-group breakdown; builds the matching per-group
    breakdown on the Python reference side; then compares both the total and each
    force class.

    The *total* PE is the authoritative pass/fail (a sum is order- and
    grouping-invariant). The per-class rows are diagnostic: a class present on
    only one side, or one whose energies disagree, points straight at the term to
    investigate. A class that appears on both sides and disagrees beyond
    tolerance also fails the result.

    Returns
    -------
    (ok, cpp_total, ref_total, per_class)
        ``per_class`` maps force-class name -> (cpp_pe, ref_pe); a side that
        lacks a class records ``nan`` for it.
    """
    # Force groups are assigned while the System is built, so the flag must be
    # set before initialize_openmm(). Surface a clear error if the binding
    # predates this feature (i.e. the extension hasn't been rebuilt).
    if not hasattr(context, "set_separate_force_groups"):
        raise AttributeError(
            "context.set_separate_force_groups is missing -- rebuild the robosample "
            "extension after adding OpenMMContext::setSeparateForceGroups and its "
            "binding."
        )
    if not hasattr(context, "calc_openmm_potential_energy_by_group"):
        raise AttributeError(
            "context.calc_openmm_potential_energy_by_group is missing -- rebuild the "
            "robosample extension after adding "
            "OpenMMContext::computePotentialEnergyByGroup and its binding."
        )

    context.set_separate_force_groups(True)
    if not context.initialize_openmm():
        raise RuntimeError("Context.initialize_openmm() failed")

    cpp_total, cpp_groups_raw = context.calc_openmm_potential_energy_by_group()
    cpp_total = float(cpp_total)
    cpp_by_class = _normalize_cpp_groups(cpp_groups_raw)

    ref_total, ref_by_class = reference_potential_energy_by_group(
        context, prmtop_path, inpcrd_path, platform_name
    )

    if not math.isfinite(cpp_total):
        print(f"C++ OpenMM energy is not finite ({cpp_total}); diagnosing positions:")
        diagnose_positions(context)

    # ---- Total comparison (authoritative) --------------------------------
    total_abs = abs(cpp_total - ref_total)
    total_rel = total_abs / max(abs(ref_total), 1.0e-12)
    total_ok = math.isfinite(cpp_total) and (
        total_abs <= abs_tol_kj_per_mol or total_rel <= rel_tol
    )

    # ---- Per-class comparison (diagnostic) -------------------------------
    names = sorted(set(cpp_by_class) | set(ref_by_class))
    per_class: dict[str, tuple[float, float]] = {}
    groups_ok = True
    for name in names:
        c = cpp_by_class.get(name, math.nan)
        r = ref_by_class.get(name, math.nan)
        per_class[name] = (c, r)
        if math.isnan(c) or math.isnan(r):
            groups_ok = False  # present on only one side
            continue
        a = abs(c - r)
        rel = a / max(abs(r), 1.0e-12)
        if not (a <= abs_tol_kj_per_mol or rel <= rel_tol):
            groups_ok = False

    ok = total_ok and groups_ok

    if verbose:
        print("OpenMM per-force-group comparison [kJ/mol]")
        print(f"  {'force class':<24} {'C++':>16} {'reference':>16} {'abs diff':>12}")
        print(f"  {'-' * 24} {'-' * 16} {'-' * 16} {'-' * 12}")
        for name in names:
            c, r = per_class[name]
            if math.isnan(c) or math.isnan(r):
                diff_str = "MISSING"
            else:
                diff_str = f"{abs(c - r):.3e}"
            c_str = "  --  " if math.isnan(c) else f"{c:.6f}"
            r_str = "  --  " if math.isnan(r) else f"{r:.6f}"
            print(f"  {name:<24} {c_str:>16} {r_str:>16} {diff_str:>12}")
        print(f"  {'-' * 24} {'-' * 16} {'-' * 16} {'-' * 12}")
        print(
            f"  {'TOTAL':<24} {cpp_total:>16.6f} {ref_total:>16.6f} {total_abs:>12.3e}"
        )
        print(f"  total rel diff   : {total_rel:.3e}")
        print(f"  per-group match  : {'PASS' if groups_ok else 'FAIL'}")
        print(f"  result           : {'PASS' if ok else 'FAIL'}")

    if not ok:
        raise ValueError(
            f"Energy mismatch: C++ OpenMM = {cpp_total:.6f} kJ/mol, Python reference = {ref_total:.6f} kJ/mol"
        )

    return ok, cpp_total, ref_total, per_class
