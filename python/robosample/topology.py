"""
topology.py
"""

from typing import NamedTuple

from attr import dataclass


@dataclass(frozen=True)
class _FieldSpec:
    """
    Declares how one list field is copied from a MoleculePrototype into the
    flat system topology.

    Attributes
    ----------
    proto_attr : str
        Attribute name on MoleculePrototype.
    sys_attr : str
        Attribute name on the system topology object.  Differs from proto_attr
        when naming conventions diverge (e.g. ``exclusions_j`` -> ``exclusion_j``
        or ``urey_bradley_j`` -> ``urey_bradley_k``).
    atom_offset : bool
        When True every element is shifted by the cumulative atom count at the
        start of the current molecule instance.  Set this for every field that
        stores atom indices (bonds, angles, torsions, z-matrix, ...).
    """

    proto_attr: str
    sys_attr: str
    atom_offset: bool = False


class _RangeSpec(NamedTuple):
    """
    Declares one begin / end range pair that records where a molecule
    instance's data starts and ends in the corresponding flat array.

    Attributes
    ----------
    counter : str
        Key in the local ``counters`` dict that accumulates the running total.
    proto_count : str
        Attribute on MoleculePrototype that gives this instance's increment.
        Usually equals counter, but z_matrix reuses ``num_atoms`` while keeping
        a separate counter so z_matrix_begin / end stay independent.
    begin : str
        Name of the begin list on the system topology.
    end : str
        Name of the end list on the system topology.
    """

    counter: str
    proto_count: str
    begin: str
    end: str


# ---------------------------------------------------------------------------
# Range schema -- one entry per interaction category
# ---------------------------------------------------------------------------

_RANGE_SPECS: tuple[_RangeSpec, ...] = (
    _RangeSpec("num_atoms", "num_atoms", "atoms_begin", "atoms_end"),
    _RangeSpec("num_bonds", "num_bonds", "bonds_begin", "bonds_end"),
    _RangeSpec("num_angles", "num_angles", "angles_begin", "angles_end"),
    _RangeSpec(
        "num_periodic_torsions",
        "num_periodic_torsions",
        "periodic_torsions_begin",
        "periodic_torsions_end",
    ),
    _RangeSpec(
        "num_harmonic_torsions",
        "num_harmonic_torsions",
        "harmonic_torsions_begin",
        "harmonic_torsions_end",
    ),
    _RangeSpec(
        "num_urey_bradley", "num_urey_bradley", "urey_bradley_begin", "urey_bradley_end"
    ),
    _RangeSpec("num_scaling14", "num_scaling14", "scaling14_begin", "scaling14_end"),
    _RangeSpec("num_exclusions", "num_exclusions", "exclusion_begin", "exclusion_end"),
    # z_matrix has exactly one row per atom; proto_count reuses num_atoms but
    # the counter is tracked independently so begin/end stay their own arrays.
    _RangeSpec("num_z_matrix_rows", "num_atoms", "z_matrix_begin", "z_matrix_end"),
)


# ---------------------------------------------------------------------------
# Field schema -- one entry per list attribute
# ---------------------------------------------------------------------------

_FIELD_SPECS: tuple[_FieldSpec, ...] = (
    # -- Atom properties (physical quantities; no index offset) ---------------
    _FieldSpec("atoms_charge", "atoms_charge"),
    _FieldSpec("atoms_mass", "atoms_mass"),
    _FieldSpec("atoms_sigma", "atoms_sigma"),
    _FieldSpec("atoms_epsilon", "atoms_epsilon"),
    _FieldSpec("atoms_radius", "atoms_radius"),
    _FieldSpec("atoms_screen", "atoms_screen"),
    _FieldSpec("atoms_x", "atoms_x"),
    _FieldSpec("atoms_y", "atoms_y"),
    _FieldSpec("atoms_z", "atoms_z"),
    # -- Atom topology indices (compound-local / table indices; NO offset) ----
    # atoms_compound_atom_index is the SimTK Compound::AtomIndex (0..n-1 within
    # the compound) and is provided by MoleculePrototype directly.
    # atoms_class_index / atoms_charged_atom_type_index are stamped onto each
    # prototype by Context.load_amber before the flattening loop runs.  All
    # three are table/compound-local indices, NOT system atom indices, so they
    # must NOT receive the atom offset.
    _FieldSpec(
        "atoms_compound_atom_index", "atoms_compound_atom_index", atom_offset=False
    ),
    _FieldSpec("atoms_class_index", "atoms_class_index", atom_offset=False),
    _FieldSpec(
        "atoms_charged_atom_type_index",
        "atoms_charged_atom_type_index",
        atom_offset=False,
    ),
    # Parallel per-atom NAME strings for the class / charged-type tables.
    # These let the consumer build the DuMM class & charged-atom-type registries
    # (unique index -> name); strings, so no offset.
    _FieldSpec("atoms_class_names", "atoms_class_names", atom_offset=False),
    _FieldSpec(
        "atoms_charged_type_names", "atoms_charged_type_names", atom_offset=False
    ),
    # Element labels and 0-based non-bonded (LJ) TYPE index.  element_* are
    # strings; nonbonded_index is a per-type index, not an atom index -- none
    # receive the atom offset.
    _FieldSpec("atoms_element_name", "atoms_element_name", atom_offset=False),
    _FieldSpec("atoms_element_symbol", "atoms_element_symbol", atom_offset=False),
    _FieldSpec("atoms_nonbonded_index", "atoms_nonbonded_index", atom_offset=False),
    # NOTE: atoms_unique_name is intentionally NOT a field spec.  Its embedded
    # residue/atom numbers are GLOBAL (whole-system, prmtop) values that depend
    # on each molecule instance's position, which a per-prototype array cannot
    # encode, so Context.load_amber assembles it directly in the flattening loop.
    # -- Bonds (i / j are atom indices; stiffness / equilibrium are scalars) --
    _FieldSpec("bonds_i", "bonds_i", atom_offset=True),
    _FieldSpec("bonds_j", "bonds_j", atom_offset=True),
    _FieldSpec("bonds_stiffness", "bonds_stiffness"),
    _FieldSpec("bonds_equilibrium", "bonds_equilibrium"),
    # -- Angles ---------------------------------------------------------------
    _FieldSpec("angles_i", "angles_i", atom_offset=True),
    _FieldSpec("angles_j", "angles_j", atom_offset=True),
    _FieldSpec("angles_k", "angles_k", atom_offset=True),
    _FieldSpec("angles_stiffness", "angles_stiffness"),
    _FieldSpec("angles_equilibrium", "angles_equilibrium"),
    # -- Periodic torsions ----------------------------------------------------
    _FieldSpec("periodic_torsions_i", "periodic_torsions_i", atom_offset=True),
    _FieldSpec("periodic_torsions_j", "periodic_torsions_j", atom_offset=True),
    _FieldSpec("periodic_torsions_k", "periodic_torsions_k", atom_offset=True),
    _FieldSpec("periodic_torsions_l", "periodic_torsions_l", atom_offset=True),
    _FieldSpec("periodic_torsions_n", "periodic_torsions_n"),
    _FieldSpec("periodic_torsions_phase", "periodic_torsions_phase"),
    _FieldSpec("periodic_torsions_stiffness", "periodic_torsions_stiffness"),
    _FieldSpec("periodic_torsions_improper", "periodic_torsions_improper"),
    # -- Harmonic (CHARMM-style improper) torsions ----------------------------
    _FieldSpec("harmonic_torsions_i", "harmonic_torsions_i", atom_offset=True),
    _FieldSpec("harmonic_torsions_j", "harmonic_torsions_j", atom_offset=True),
    _FieldSpec("harmonic_torsions_k", "harmonic_torsions_k", atom_offset=True),
    _FieldSpec("harmonic_torsions_l", "harmonic_torsions_l", atom_offset=True),
    _FieldSpec("harmonic_torsions_stiffness", "harmonic_torsions_stiffness"),
    _FieldSpec("harmonic_torsions_phase", "harmonic_torsions_phase"),
    # -- Urey-Bradley ---------------------------------------------------------
    # MoleculePrototype names the second 1,3-atom "urey_bradley_j".
    # The system topology calls it "urey_bradley_k" (angle-position naming).
    _FieldSpec("urey_bradley_i", "urey_bradley_i", atom_offset=True),
    _FieldSpec("urey_bradley_k", "urey_bradley_k", atom_offset=True),
    _FieldSpec("urey_bradley_stiffness", "urey_bradley_stiffness"),
    _FieldSpec("urey_bradley_equilibrium", "urey_bradley_equilibrium"),
    # -- 1-4 non-bonded scaling pairs -----------------------------------------
    _FieldSpec("scaling14_i", "scaling14_i", atom_offset=True),
    _FieldSpec("scaling14_l", "scaling14_l", atom_offset=True),
    _FieldSpec("scaling14_charge_product", "scaling14_charge_product"),
    _FieldSpec("scaling14_epsilon", "scaling14_epsilon"),
    _FieldSpec("scaling14_sigma", "scaling14_sigma"),
    # -- Exclusions -----------------------------------------------------------
    _FieldSpec("exclusion_i", "exclusion_i", atom_offset=True),
    _FieldSpec("exclusion_j", "exclusion_j", atom_offset=True),
    # -- Z-matrix (sentinel-aware atom-index offset) --------------------------
    # z_matrix_i is always a valid atom index; no sentinel handling needed.
    # z_matrix_j / k / l use -1 as a sentinel for root-triplet padding rows;
    # those entries must be left at -1 and NOT shifted.
    _FieldSpec("z_matrix_i", "z_matrix_i", atom_offset=True),
    _FieldSpec("z_matrix_j", "z_matrix_j", atom_offset=True),
    _FieldSpec("z_matrix_k", "z_matrix_k", atom_offset=True),
    _FieldSpec("z_matrix_l", "z_matrix_l", atom_offset=True),
)
