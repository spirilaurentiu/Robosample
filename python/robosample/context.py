from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import amber_loader, prmtop_reader, topology
from .amber_dihedral_classifier import AmberDihedralClassifier
from .amber_dihedral_types import DihedralType
from .molecule_prototype import MoleculePrototype
from .robo_bindings import (
    AcceptRejectMode,
    JointType,
    NonbondedMethod,
    SystemTopology,
)
from .robo_bindings import Context as _Context
from .secondary_structure import DSSPCode


@dataclass
class RigidBodyAssignment:
    atom_body_idx: np.ndarray  # (n_atoms,) int64  -- body index per atom
    n_bodies: int
    atom_is_inboard: (
        np.ndarray
    )  # (n_atoms,) bool   -- atom is inboard across a PIN joint
    atom_is_outboard: (
        np.ndarray
    )  # (n_atoms,) bool   -- atom is outboard         across a PIN joint
    shake_pairs: list[tuple[int, int]] = field(default_factory=list)


class Context(_Context):
    system_topology: SystemTopology
    dihedral_classifier: AmberDihedralClassifier
    df_bonds: pd.DataFrame

    def __init__(
        self, base_name: str, seed: int, dihedral_classifier: AmberDihedralClassifier
    ) -> None:
        super().__init__(base_name=base_name, seed=seed)
        self.dihedral_classifier = dihedral_classifier

    @property
    def standard_dihedral_bonds(self) -> pd.DataFrame:
        """The per-molecule bond table (``df_bonds``).

        Columns of interest for flexibility selection:
          * ``atom1_idx`` / ``atom2_idx`` -- molecule-LOCAL (prototype) atom
            indices of the bond's two atoms.
          * ``molecule_idx``              -- which molecule instance the bond is in.
          * ``dihedral_type``             -- e.g. phi / psi / omega / chi.
          * ``is_ring_closing``           -- ring-closing bonds stay rigid.

        Filter this by ``dihedral_type`` (and optionally ``molecule_idx``) and
        pass the result to :meth:`build_flexibilities`, e.g.::

            bonds = ctx.standard_dihedral_bonds.loc[
                ctx.standard_dihedral_bonds["dihedral_type"].isin(["phi", "psi"])
            ]
            sele = ctx.build_flexibilities(bonds, rb.JointType.Torsion, False)
        """
        return self.df_bonds

    def build_flexibilities(self, bonds, mobility, flag):
        """Choose which bonds are flexible (given ``mobility``) vs rigid.

        ``bonds`` may be:
          * ``None``            -- every eligible (non-ring, non-terminal) bond.
          * a ``pandas.DataFrame`` of ``standard_dihedral_bonds`` rows -- only
            those bonds become flexible. ``atom1_idx``/``atom2_idx`` are
            molecule-local and are mapped to global atom indices here via
            ``system_topology.atoms_begin[molecule_idx]`` (matching how the C++
            topology stores ``bonds_i``/``bonds_j``). Ring-closing rows are
            skipped.
          * an iterable of ``(i, j)`` GLOBAL atom-index pairs -- passed through.

        The actual rigid/flexible decomposition is done in C++ (the base class):
        ring-closing bonds always stay rigid, and only non-terminal bonds are
        eligible.
        """
        if bonds is None:
            bonds = []

        if isinstance(bonds, pd.DataFrame):
            atoms_begin = list(self.system_topology.atoms_begin)
            pairs: list[tuple[int, int]] = []
            for row in bonds.itertuples(index=False):
                if bool(getattr(row, "is_ring_closing", False)):
                    continue
                off = atoms_begin[int(row.molecule_idx)]
                pairs.append((off + int(row.atom1_idx), off + int(row.atom2_idx)))
        else:
            pairs = [(int(i), int(j)) for (i, j) in bonds]

        return super().build_flexibilities(pairs, mobility, flag)

    def load_amber(
        self,
        prmtop_path: str | os.PathLike[str],
        inpcrd_path: str | os.PathLike[str],
        *,
        use_gbsa_obc2: bool | None = None,
        nonbonded_method: NonbondedMethod | None = None,
        nonbonded_cutoff: float | None = None,
        ewald_error_tolerance: float | None = None,
    ) -> None:
        """
        Read an AMBER ``.prmtop`` + coordinate file and populate
        ``self.system_topology``.

        Solvent model (auto-detected, OpenMM-style)
        -------------------------------------------
        The solvent model follows the periodic box the same way OpenMM's
        ``createSystem`` does -- it is read from the topology/coordinates, not
        passed by hand:

          * **Box present** (a solvated system: ``IFBOX`` > 0 / the rst7 carries
            box vectors) -> explicit solvent. PME electrostatics under periodic
            boundary conditions, GBSA off, the box stored as three reduced
            lattice vectors (nm) on ``system_topology.box_vectors``, a 1.0 nm
            cutoff by default, and **every** molecule given a FREE (6-DOF) root
            so solvent can translate and reorient (a welded water would be frozen
            in place).
          * **No box** -> implicit solvent. GBSA-OBC2 with a non-periodic,
            no-cutoff nonbonded treatment (the historical Robosample default).

        ``use_gbsa_obc2`` overrides the implicit-solvent force:
          * ``None`` (default) -- auto: GBSA-OBC2 on when there is no box, off
            when there is one.
          * ``False`` -- force GBSA-OBC2 off (e.g. a gas-phase run with no
            implicit solvent at all).
          * ``True``  -- force GBSA-OBC2 on. This is implicit solvation and is
            mutually exclusive with a periodic box, so if the system carries a
            box (explicit-solvent PME) this raises ``ValueError``.

        The ``nonbonded_method`` / ``nonbonded_cutoff`` / ``ewald_error_tolerance``
        arguments override the auto-selected defaults when given. A periodic
        ``nonbonded_method`` without a box, or a non-periodic one with a box, is a
        contradiction and raises.
        """
        self.system_topology = SystemTopology()
        self.system_topology.gbsa_solute_dielectric = 1.0
        self.system_topology.gbsa_solvent_dielectric = 78.5

        # Parse files and do typing. Fully ParmEd-free as of fast-loader Step
        # 4b (docs/specs/fast-amber-loader.md): dedup (a former Structure.split()
        # call) and the per-prototype MoleculePrototype input (a former sliced
        # ParmEd Structure) were already replaced in Step 4a by
        # amber_loader.partition_molecules + amber_loader.PrototypeTopology
        # below; box-vector construction and virtual-site (ExtraPoint) frame
        # detection -- the last two ParmEd usages in this method -- are
        # replaced below by amber_loader.box_vectors_from_lengths_angles and
        # amber_loader.extract_virtual_sites, both driven off raw_data/
        # atom_coords_nm alone. No topology-file load through ParmEd remains.
        parm_file = prmtop_reader.parse_prmtop(prmtop_path)
        raw_data = parm_file["raw_data"]

        # Direct, ParmEd-free coordinate read (component (A) of the
        # fast-loader spec): a single vectorized (N,3) array in prmtop atom
        # order, instead of per-atom `parm.atoms[i].xx/.xy/.xz` Python
        # attribute access. `rst7_box` -- [a, b, c, alpha, beta, gamma] (nm,
        # rad) or None -- is the raw box GEOMETRY parsed from the coordinate
        # file's box line, if any; whether it is actually USED is decided
        # below from IFBOX (POINTERS), not from its mere presence (see the
        # periodicity block) -- some non-periodic rst7 files (e.g. a stripped
        # CHAMBER system) carry a degenerate placeholder box line that must
        # be ignored when IFBOX == 0.
        atom_coords_nm, rst7_box = amber_loader.read_amber_coordinates(inpcrd_path)
        natom_prmtop = len(raw_data["ATOM_NAME"])
        if atom_coords_nm.shape[0] != natom_prmtop:
            raise ValueError(
                f"Coordinate file '{inpcrd_path}' has {atom_coords_nm.shape[0]} atoms, "
                f"but topology '{prmtop_path}' has {natom_prmtop}."
            )

        # NTYPES (POINTERS[1], 0-based) -- number of distinct LJ atom types.
        self.system_topology.num_nb_types = int(raw_data["POINTERS"][1])
        self.system_topology.a_coef, self.system_topology.b_coef = (
            prmtop_reader.load_lj_coefs(
                raw_data,
                self.system_topology.num_nb_types,
            )
        )

        # NBFIX detection (docs/specs/fast-amber-loader.md §4a + the GAP
        # flagged in docs/specs/loader-feature-matrix.md "Nonbonded" table):
        # has_nbfix_fast reproduces OpenMM's own combining-rule-deviation
        # check on the raw (unconverted) LJ A/B coefficient table. Setting
        # this flag is what makes OpenMMContext.cpp::initialize's existing
        # `if (systemTopology.hasNBfix) { createCustomNonbondedForce(...) }`
        # guard fire (that function unconditionally raises "not implemented
        # yet") -- an NBFIX prmtop must raise there instead of silently
        # getting combining-rule LJ.
        self.system_topology.has_nb_fix = bool(
            prmtop_reader.has_nbfix_fast(
                raw_data["NONBONDED_PARM_INDEX"],
                self.system_topology.num_nb_types,
                raw_data["LENNARD_JONES_ACOEF"],
                raw_data["LENNARD_JONES_BCOEF"],
            )
        )

        # Molecule-instance / prototype dedup (component (B) of the
        # fast-loader spec, docs/specs/fast-amber-loader.md §3(B)/§5):
        # replaces parm.split() with scipy connected-components over the bond
        # graph plus a conservative per-instance fingerprint. Energy is
        # invariant to prototype grouping/order (spec §5), so correctness
        # only requires (a) never merging non-identical molecules and (b)
        # numbering instances the same way ParmEd does -- both are gated by
        # the loader differential test.
        partition = amber_loader.partition_molecules(raw_data, natom_prmtop)

        # Whole-system (prmtop order, 0-based) residue index and residue/atom
        # name lookup. Used below both to build each prototype's
        # PrototypeTopology (residue grouping) AND, after the instance loop,
        # to build atoms_unique_name in one vectorized pass instead of one
        # f-string call per atom. RESIDUE_POINTER holds each residue's 1-based
        # first-atom index; searchsorted on the 0-based starts gives, for
        # every prmtop atom, the 0-based index of the residue it belongs to.
        residue_starts0 = np.asarray(raw_data["RESIDUE_POINTER"], dtype=np.int64) - 1
        residue_of_atom = (
            np.searchsorted(residue_starts0, np.arange(natom_prmtop), side="right") - 1
        )
        resname_by_residue = np.asarray(raw_data["RESIDUE_LABEL"])
        resname_of_atom = resname_by_residue[residue_of_atom]
        atom_name_of_atom = np.asarray(raw_data["ATOM_NAME"])

        # Build one MoleculePrototype per unique prototype (parsing is
        # expensive) from an amber_loader.PrototypeTopology built directly
        # from the whole-system raw prmtop arrays at the representative
        # (first-occurring) instance's atom indices -- component (C) of the
        # fast-loader spec (docs/specs/fast-amber-loader.md §3(C)). Replaces
        # the ParmEd `parm[mask]` slice (Structure.__getitem__ boolean mask),
        # which rebuilds/prunes a full ParmEd object graph over the WHOLE
        # system on every call and was the dominant remaining load-time cost.
        molecule_prototypes = []
        for atom_indices in partition.prototype_representative_atoms:
            proto_topology = amber_loader.PrototypeTopology(
                raw_data,
                atom_indices,
                atom_coords_nm,
                residue_of_atom,
                resname_by_residue,
            )
            molecule_prototypes.append(
                MoleculePrototype(proto_topology, self.dihedral_classifier)
            )

        # Per-prototype BFS/compound->local (prmtop-order) permutation as a
        # numpy array, cached once per UNIQUE prototype rather than rebuilt
        # from the Python list on every instance (a water prototype alone can
        # recur tens of thousands of times).
        compound_to_local_by_proto: list[np.ndarray] = [
            np.asarray(proto.compound_to_local, dtype=np.int64)
            for proto in molecule_prototypes
        ]

        # Flatten into (instance_index, prototype_index), ordered by instance
        # (already ordered by construction; .sort() is a defensive no-op
        # matching the pre-rewrite code's explicit ordering guarantee).
        self.molecules: list[tuple[int, int]] = list(
            enumerate(partition.prototype_of_instance)
        )
        self.molecules.sort(key=lambda x: x[0])

        COLUMNS: dict[str, str | type] = {
            "atom1_idx": np.int32,
            "atom2_idx": np.int32,
            "molecule_idx": np.int32,
            "dihedral_type": pd.CategoricalDtype(categories=list(DihedralType)),
            "is_ring_closing": bool,
            "dss": pd.CategoricalDtype(categories=list(DSSPCode)),
            "residue_name": str,
            "residue_idx": np.int32,
        }
        # Bond-table rows are accumulated column-wise here and assembled into a
        # single DataFrame after the instance loop.  Building the frame
        # incrementally with pd.concat() per bond is O(N^2) -- each concat copies
        # the whole accumulated frame -- and a large solvated system has >100k
        # bonds, which is a dominant part of the multi-second load time.
        bond_cols: dict[str, list] = {col: [] for col in COLUMNS}

        # ---- Local accumulators ----------------------------------------------
        acc: dict[str, list] = {}
        for spec in topology._RANGE_SPECS:
            acc[spec.begin] = []
            acc[spec.end] = []
        for fspec in topology._FIELD_SPECS:
            acc.setdefault(fspec.sys_attr, [])

        counters: dict[str, int] = {spec.counter: 0 for spec in topology._RANGE_SPECS}

        # atoms_unique_name is assembled here rather than via _FIELD_SPECS: it
        # embeds GLOBAL (whole-system, prmtop) residue and atom numbers, which
        # depend on each instance's position and so cannot live on a shared
        # prototype. Built in ONE vectorized pass after the instance loop (see
        # below) from atoms_prmtop_index + the whole-system resname_of_atom/
        # atom_name_of_atom arrays precomputed above -- not per-atom f-strings.
        acc.setdefault("atoms_unique_name", [])
        acc.setdefault("atoms_x", [])
        acc.setdefault("atoms_y", [])
        acc.setdefault("atoms_z", [])

        # BFS/compound position -> 0-based prmtop atom index. Appended in c
        # (compound) order below, so the list position IS the global BFS index;
        # transmitted to C++ to write trajectory output in prmtop atom order.
        acc.setdefault("atoms_prmtop_index", [])

        # Per-molecule arrays (one entry per instance, NOT per atom):
        #   atoms_root_index : the root atom's index in the GLOBAL, BFS-ordered
        #                      atoms vector.  The root is compound index 0 within
        #                      its molecule, so its global index is just
        #                      atom_off + atoms_root_compound_index (== atom_off).
        #                      This is a compound/global index, never a prmtop one.
        #   root_mobilities  : how each molecule's root body attaches to ground.
        # POLICY: every molecule defaults to a freely-moving (6-DOF) root.  This
        # is a simulation choice, not a property of the molecule -- adjust here
        # (per molecule or per type) if some molecules should be welded/pinned.
        acc.setdefault("atoms_root_index", [])
        acc.setdefault("root_mobilities", [])
        default_root_mobility = JointType.Rigid

        for instance_idx, prototype_idx in self.molecules:
            proto: MoleculePrototype = molecule_prototypes[prototype_idx]
            atom_off: int = counters["num_atoms"]

            # 1. Begin markers
            for spec in topology._RANGE_SPECS:
                acc[spec.begin].append(counters[spec.counter])

            # 1b. Per-molecule root: global index of this instance's root atom
            #     (compound 0 of the molecule) and its root mobility.
            acc["atoms_root_index"].append(atom_off + proto.atoms_root_compound_index)
            acc["root_mobilities"].append(default_root_mobility)

            # 2. Field data (includes atoms_compound_atom_index,
            #    atoms_class_index, atoms_charged_atom_type_index once those are
            #    registered in topology._FIELD_SPECS with atom_offset=False).
            for fspec in topology._FIELD_SPECS:
                src: list = getattr(proto, fspec.proto_attr)
                dst: list = acc[fspec.sys_attr]
                if not fspec.atom_offset:
                    dst.extend(src)
                else:
                    dst.extend(x + atom_off for x in src)

            # 2c. Per-instance atom bookkeeping, vectorized over the whole
            #     instance instead of one Python-level iteration per atom:
            #     prmtop index = atom_off + local (prmtop-order) index within
            #     the molecule; global (compound/BFS) index = atom_off + c,
            #     i.e. arange(atom_off, atom_off + num_atoms) since atoms are
            #     appended in compound order.
            local_idx_arr = compound_to_local_by_proto[prototype_idx]
            prmtop_idx_arr = atom_off + local_idx_arr

            # BFS->prmtop permutation for trajectory output (position in this
            # list == global_idx, value == prmtop_idx); also doubles as the
            # per-instance index set used below to gather coordinates and,
            # after the loop, to build atoms_unique_name and
            # prmtop_to_global_index in one vectorized pass.
            acc["atoms_prmtop_index"].extend(prmtop_idx_arr.tolist())

            # Per-INSTANCE coordinates: gathered from the whole-system
            # coordinate array at this instance's prmtop indices, NOT from the
            # shared prototype (which would superimpose every copy of a
            # molecule type and blow the energy up to inf).
            inst_coords_nm = atom_coords_nm[prmtop_idx_arr]
            acc["atoms_x"].extend(inst_coords_nm[:, 0].tolist())
            acc["atoms_y"].extend(inst_coords_nm[:, 1].tolist())
            acc["atoms_z"].extend(inst_coords_nm[:, 2].tolist())

            # 3. Advance counters
            for spec in topology._RANGE_SPECS:
                counters[spec.counter] += getattr(proto, spec.proto_count)

            # 4. End markers
            for spec in topology._RANGE_SPECS:
                acc[spec.end].append(counters[spec.counter])

            # 5. Accumulate bond-table rows column-wise (assembled into df_bonds
            #    once after the loop; see bond_cols above).
            for i in range(proto.num_bonds):
                bond_cols["atom1_idx"].append(proto.bonds_i[i])
                bond_cols["atom2_idx"].append(proto.bonds_j[i])
                bond_cols["molecule_idx"].append(instance_idx)
                bond_cols["dihedral_type"].append(proto.bonds_dihedral_type[i])
                bond_cols["is_ring_closing"].append(proto.bonds_is_ring_closing[i])
                bond_cols["dss"].append(proto.bonds_secondary_structure[i])
                bond_cols["residue_name"].append("UNK")
                bond_cols["residue_idx"].append(-1)

        # Assemble the bond table in one shot.  Building column-wise and
        # constructing the DataFrame once is O(N) in total bonds, versus the
        # O(N^2) per-bond pd.concat it replaces; column contents are identical.
        self.df_bonds = pd.DataFrame(bond_cols)

        # ---- atoms_unique_name + prmtop_to_global_index (vectorized, once) ---
        # acc["atoms_prmtop_index"] is now the full global (compound/BFS-order)
        # -> prmtop-index permutation built by the instance loop above. Format
        # unchanged: "{resname}{res}_{atomname}_{atom}", e.g. "ALA1_N_4", both
        # numbers GLOBAL and 1-based (res = whole-system 0-based residue index
        # of the atom + 1; atom = the atom's 1-based prmtop index).
        atoms_prmtop_index_arr = np.asarray(acc["atoms_prmtop_index"], dtype=np.int64)
        resname_arr = resname_of_atom[atoms_prmtop_index_arr]
        resnum_arr = (residue_of_atom[atoms_prmtop_index_arr] + 1).astype(str)
        atomname_arr = atom_name_of_atom[atoms_prmtop_index_arr]
        atomnum_arr = (atoms_prmtop_index_arr + 1).astype(str)
        unique_names = np.char.add(
            np.char.add(
                np.char.add(np.char.add(resname_arr, resnum_arr), "_"), atomname_arr
            ),
            np.char.add("_", atomnum_arr),
        )
        acc["atoms_unique_name"] = unique_names.tolist()

        # prmtop_to_global_index[prmtop_idx] = global (compound/BFS) index;
        # used below by the CMAP section and the virtual-site block. Every
        # prmtop atom belongs to exactly one instance/compound position, so
        # atoms_prmtop_index_arr is a permutation of range(natom_prmtop) and a
        # plain array inverse is well-defined (see docs/specs/
        # fast-amber-loader.md §3(D): "a dict or an int array is fine").
        self.prmtop_to_global_index = np.empty(natom_prmtop, dtype=np.int64)
        self.prmtop_to_global_index[atoms_prmtop_index_arr] = np.arange(
            atoms_prmtop_index_arr.shape[0], dtype=np.int64
        )

        # ---- Flush accumulators to topology ----------------------------------
        for attr, data in acc.items():
            setattr(self.system_topology, attr, data)
        for counter, value in counters.items():
            setattr(self.system_topology, counter, value)
        self.system_topology.num_molecules = len(self.molecules)

        # ---- CMAP correction maps and torsions -------------------------------
        # Each CMAP map defines a 2D periodic energy surface E(phi, psi) for a
        # pair of backbone dihedrals. The raw prmtop stores energies in kcal/mol
        # with psi varying fastest, starting at -180°. OpenMM expects kJ/mol
        # with phi varying fastest, starting at 0°. We perform both conversions.
        #
        # SoA layout consumed by OpenMMContext::createCMAPTorsionForce:
        #   cmap_grid_size          int            uniform grid resolution (e.g. 24)
        #   cmap_grid_energy        list[float]    all maps concatenated,
        #                                          res*res values per map, phi-fastest
        #   num_cmap_grids          int
        #   num_cmap_torsions       int
        #   cmap_torsions_map_index list[int]      0-based map index per torsion
        #   cmap_torsions_a[i/j/k/l]  list[int]   torsion A atom global indices
        #   cmap_torsions_b[i/j/k/l]  list[int]   torsion B atom global indices
        # (raw_data was already bound near the top of this method, from the
        # same parm_file dict.)
        # CMAP flag names carry a "CHARMM_" prefix in ParmEd/modern AmberTools
        # `chamber` output (see ParmEd's ChamberParm._cmap_prefix == "CHARMM_")
        # -- added to disambiguate a CHAMBER file's CMAP section from ff19SB's
        # own native (non-CHARMM) CMAP corrections, which use the same
        # unprefixed flag names. Older CHAMBER-format prmtop files (e.g.
        # examples/GfcDstrippedMin.prmtop) predate that rename and still use
        # the unprefixed names. Accept either.
        def _cmap_flag(name: str) -> str:
            prefixed = "CHARMM_" + name
            return prefixed if prefixed in raw_data else name

        cmap_resolution: list[int] = raw_data.get(_cmap_flag("CMAP_RESOLUTION"), [])
        num_cmap_grids = len(cmap_resolution)

        cmap_grid_size: int = 0
        cmap_grid_energy: list[float] = []

        for grid_idx in range(num_cmap_grids):
            res = int(cmap_resolution[grid_idx])
            if cmap_grid_size == 0:
                cmap_grid_size = res
            elif cmap_grid_size != res:
                raise ValueError(
                    f"Mixed CMAP grid resolutions ({cmap_grid_size} vs {res}); "
                    "this system requires per-map resolution support not yet implemented."
                )

            key = _cmap_flag(f"CMAP_PARAMETER_{grid_idx + 1:02d}")
            raw_cmap: list[float] = raw_data[key]

            # Amber/ParmEd layout: psi varies fastest, phi starts at -180°.
            # OpenMM layout:       phi varies fastest, phi starts at   0°.
            # The shift of half the grid (res//2 bins) rotates 0° → 180°.
            half = res // 2
            flat_omm = [0.0] * (res * res)
            for phi_omm in range(res):
                for psi_omm in range(res):
                    phi_amber = (phi_omm + half) % res
                    psi_amber = (psi_omm + half) % res
                    old_index = phi_amber * res + psi_amber  # psi-fastest in Amber
                    new_index = phi_omm + res * psi_omm  # psi-fastest in OpenMM
                    flat_omm[new_index] = raw_cmap[old_index] * 4.184  # kcal→kJ

            cmap_grid_energy.extend(flat_omm)

        # Build the per-torsion arrays.  The prmtop CMAP_INDEX section stores
        # 6 integers per torsion:
        #   [atom1, atom2, atom3, atom4, atom5, map_index]
        # where atoms 1-4 define torsion A (phi) and atoms 2-5 define torsion B
        # (psi), with atoms 2-4 shared.  All atom indices are 1-based in prmtop.
        # We convert to global (BFS) indices via prmtop_to_global_index.
        cmap_index_raw: list[int] = raw_data.get(_cmap_flag("CMAP_INDEX"), [])
        num_cmap_torsions = len(cmap_index_raw) // 6

        cmap_torsions_map_index: list[int] = []
        cmap_torsions_ai: list[int] = []
        cmap_torsions_aj: list[int] = []
        cmap_torsions_ak: list[int] = []
        cmap_torsions_al: list[int] = []
        cmap_torsions_bi: list[int] = []
        cmap_torsions_bj: list[int] = []
        cmap_torsions_bk: list[int] = []
        cmap_torsions_bl: list[int] = []

        def _g(prmtop_1based: int) -> int:
            """Convert 1-based prmtop atom index to global (BFS) index."""
            return int(self.prmtop_to_global_index[prmtop_1based - 1])

        for i in range(0, len(cmap_index_raw), 6):
            a1, a2, a3, a4, a5, map_idx = cmap_index_raw[i : i + 6]
            if map_idx < 1 or map_idx > num_cmap_grids:
                raise ValueError(
                    f"CMAP torsion {i // 6}: map index {map_idx} out of range "
                    f"[1, {num_cmap_grids}]"
                )
            cmap_torsions_map_index.append(map_idx - 1)  # 0-based map index
            cmap_torsions_ai.append(_g(a1))
            cmap_torsions_aj.append(_g(a2))
            cmap_torsions_ak.append(_g(a3))
            cmap_torsions_al.append(_g(a4))
            cmap_torsions_bi.append(_g(a2))
            cmap_torsions_bj.append(_g(a3))
            cmap_torsions_bk.append(_g(a4))
            cmap_torsions_bl.append(_g(a5))

        # Attribute names mirror TopologyElements.hpp camelCase -> snake_case:
        #   cmapGridSize        -> cmap_grid_size
        #   cmapGridEnergy      -> cmap_grid_energy
        #   cmapTorsionMapIndex -> cmap_torsion_map_index   (singular "Torsion")
        #   cmapTorsionA1..A4   -> cmap_torsion_a1..a4     (1-based numeric suffix)
        #   cmapTorsionB1..B4   -> cmap_torsion_b1..b4
        self.system_topology.cmap_grid_size = cmap_grid_size
        self.system_topology.cmap_grid_energy = cmap_grid_energy
        self.system_topology.cmap_torsion_map_index = cmap_torsions_map_index
        self.system_topology.cmap_torsion_a1 = cmap_torsions_ai
        self.system_topology.cmap_torsion_a2 = cmap_torsions_aj
        self.system_topology.cmap_torsion_a3 = cmap_torsions_ak
        self.system_topology.cmap_torsion_a4 = cmap_torsions_al
        self.system_topology.cmap_torsion_b1 = cmap_torsions_bi
        self.system_topology.cmap_torsion_b2 = cmap_torsions_bj
        self.system_topology.cmap_torsion_b3 = cmap_torsions_bk
        self.system_topology.cmap_torsion_b4 = cmap_torsions_bl

        # NOTE: Urey-Bradley (ureyBradleyI / ureyBradleyK / stiffness /
        # equilibrium) is per-molecule data sourced from MoleculePrototype and
        # is already populated by the _FIELD_SPECS loop above, with correct
        # atom offsets applied. No inline block needed here.

        # --------------------------------------------------------------------
        #  Virtual sites (extra points). 4-point waters (OPC/TIP4P) carry a
        #  massless EP whose position is a 3-particle affine average of the real
        #  atoms. amber_loader.extract_virtual_sites reads each EP's frame +
        #  weights directly from raw bond/angle arrays (fast-loader Step 4b,
        #  replacing ParmEd's ExtraPoint/ThreeParticleExtraPointFrame) in
        #  prmtop-index space; map through the same prmtop->global reindexing
        #  the rest of the SoA uses, so the indices match the OpenMM particle
        #  order. Done unconditionally: an EP must be declared whenever it
        #  exists, regardless of the solvent flag. Any frame type other than
        #  the 3-particle in-plane average (e.g. TIP5P's out-of-plane frame)
        #  raises amber_loader.UnsupportedTopologyFeature there, naming the
        #  pattern -- never silently approximated (spec §4a).
        # --------------------------------------------------------------------
        vs_records = amber_loader.extract_virtual_sites(raw_data, atom_coords_nm)
        g = self.prmtop_to_global_index
        vs_site = [int(g[r.site]) for r in vs_records]
        vs_a1 = [int(g[r.atom1]) for r in vs_records]
        vs_a2 = [int(g[r.atom2]) for r in vs_records]
        vs_a3 = [int(g[r.atom3]) for r in vs_records]
        vs_w1 = [r.weight1 for r in vs_records]
        vs_w2 = [r.weight2 for r in vs_records]
        vs_w3 = [r.weight3 for r in vs_records]

        self.system_topology.num_virtual_sites = len(vs_site)
        self.system_topology.vs_site = vs_site
        self.system_topology.vs_atom1 = vs_a1
        self.system_topology.vs_atom2 = vs_a2
        self.system_topology.vs_atom3 = vs_a3
        self.system_topology.vs_weight1 = vs_w1
        self.system_topology.vs_weight2 = vs_w2
        self.system_topology.vs_weight3 = vs_w3

        # --------------------------------------------------------------------
        #  Solvent / periodicity configuration.
        #
        #  OpenMM-style auto-detection: periodicity is decided from IFBOX
        #  (POINTERS[27], 0-based) -- exactly how OpenMM's own AmberPrmtopFile
        #  reader does it -- NOT from the mere presence of a box line in the
        #  coordinate file. This matters: some non-periodic (IFBOX == 0)
        #  systems carry a degenerate placeholder box line in their rst7
        #  (e.g. "0 0 0 90 90 90" on a stripped CHAMBER system), which must
        #  be ignored rather than treated as "has a box" -- the pre-Step-4b
        #  loader had exactly this bug (ParmEd's `parm.box = f.box` setter
        #  re-derives IFBOX from the coordinate file's box line alone,
        #  overwriting the prmtop's own IFBOX == 0). A periodic box in the
        #  topology selects explicit solvent (PME under PBC, GBSA off); its
        #  absence selects implicit solvent (GBSA-OBC2, non-periodic). The
        #  caller never passes an "explicit_solvent" flag -- it follows IFBOX.
        # --------------------------------------------------------------------
        periodic_methods = (
            NonbondedMethod.CutoffPeriodic,
            NonbondedMethod.Ewald,
            NonbondedMethod.PME,
        )
        ifbox = int(raw_data["POINTERS"][27])
        has_box = ifbox > 0
        if has_box and rst7_box is None:
            # No box line in the coordinate file (e.g. a restart written
            # without one): fall back to the prmtop's own legacy
            # BOX_DIMENSIONS section -- [beta_deg, a, b, c] (Angstrom),
            # ALL THREE angles taken equal to beta -- exactly the fallback
            # both ParmEd's amber.readparm.LoadParm ("if all else fails, set
            # the box from the prmtop file") and OpenMM's own
            # PrmtopLoader.getBoxBetaAndDimensions/computePeriodicBoxVectors
            # use. Only reached when IFBOX > 0 but the rst7/inpcrd carries no
            # box line (the b1-1n example exercises this).
            if "BOX_DIMENSIONS" not in raw_data:
                raise ValueError(
                    f"Topology '{prmtop_path}' declares a periodic box (IFBOX={ifbox}), but "
                    f"coordinate file '{inpcrd_path}' carries no box line and the topology has "
                    "no BOX_DIMENSIONS fallback section."
                )
            beta_deg, a_ang, b_ang, c_ang = (
                float(x) for x in raw_data["BOX_DIMENSIONS"][:4]
            )
            beta_rad = beta_deg * amber_loader.DEG_TO_RAD
            rst7_box = np.array(
                [
                    a_ang * amber_loader.ANG_TO_NM,
                    b_ang * amber_loader.ANG_TO_NM,
                    c_ang * amber_loader.ANG_TO_NM,
                    beta_rad,
                    beta_rad,
                    beta_rad,
                ],
                dtype=np.float64,
            )

        # Guard the two ways an explicit override contradicts the detected box.
        if not has_box and nonbonded_method in periodic_methods:
            raise ValueError(
                "A periodic nonbonded_method (CutoffPeriodic/Ewald/PME) was requested, but the "
                "coordinate/topology has no periodic box. Provide a solvated inpcrd/rst7 (with box "
                "information), or drop the nonbonded_method override to use implicit solvent."
            )
        if has_box and nonbonded_method is not None and nonbonded_method not in periodic_methods:
            raise ValueError(
                "The system carries a periodic box (explicit solvent), but a non-periodic "
                "nonbonded_method was requested. Explicit solvent requires a periodic method "
                "(PME/Ewald/CutoffPeriodic); drop the override to use the PME default."
            )

        if has_box:
            # Explicit solvent: PME by default, GBSA off, real box, FREE roots.
            # Implicit GBSA-OBC2 is mutually exclusive with a periodic box: real
            # waters carry the solvation, so forcing it on here is an error.
            if use_gbsa_obc2:
                raise ValueError(
                    "use_gbsa_obc2=True (implicit GBSA-OBC2 solvent) is incompatible with a "
                    "periodic box: this system is solvated, so it runs as explicit-solvent PME "
                    "under PBC and implicit solvation must be off. Remove use_gbsa_obc2=True."
                )
            self.system_topology.use_gbsa_obc2 = False
            method = (
                nonbonded_method
                if nonbonded_method is not None
                else NonbondedMethod.PME
            )
            self.system_topology.nonbonded_method = method
            self.system_topology.nonbonded_cutoff = (
                nonbonded_cutoff if nonbonded_cutoff is not None else 1.0
            )
            if ewald_error_tolerance is not None:
                self.system_topology.ewald_error_tolerance = ewald_error_tolerance

            # Box: reduce the rst7's (a, b, c, alpha, beta, gamma) geometry to
            # three lattice vectors (nm), replicating ParmEd's box_vectors
            # (ParmEd's geometry.box_lengths_and_angles_to_vectors) exactly --
            # see amber_loader.box_vectors_from_lengths_angles. Flattened
            # row-major [a.xyz b.xyz c.xyz].
            box_nm = amber_loader.box_vectors_from_lengths_angles(
                rst7_box[:3], rst7_box[3:]
            )
            self.system_topology.box_vectors = [
                float(component) for vec in box_nm for component in vec
            ]

            # Every molecule (solute and each solvent molecule) gets a FREE 6-DOF
            # root so it can translate and reorient. A welded solvent molecule
            # would be frozen in place, which is unphysical for explicit solvent.
            self.system_topology.root_mobilities = [
                JointType.Free
            ] * self.system_topology.num_molecules
        else:
            # Implicit / non-periodic path. GBSA-OBC2 is on by default and can be
            # switched off with use_gbsa_obc2=False (gas phase). Honour explicit
            # method/cutoff overrides; otherwise leave the historical defaults
            # (GBSA-OBC2, NoCutoff) intact.
            self.system_topology.use_gbsa_obc2 = (
                True if use_gbsa_obc2 is None else bool(use_gbsa_obc2)
            )
            if nonbonded_method is not None:
                self.system_topology.nonbonded_method = nonbonded_method
            if nonbonded_cutoff is not None:
                self.system_topology.nonbonded_cutoff = nonbonded_cutoff
            if ewald_error_tolerance is not None:
                self.system_topology.ewald_error_tolerance = ewald_error_tolerance

        return

    def add_docking_world(
        self, ligand_molecule_indices, mass_scale=None, reversibility_check_every=0
    ):
        """Add a rigid-body docking world.

        Parameters
        ----------
        ligand_molecule_indices : int | Iterable[int]
            Molecule index/indices treated as ligands. Their root bodies become
            Free (6 external DOF); every other molecule is welded to Ground and
            rigid. The binding-site centre is the centroid of all non-ligand
            (receptor) atoms. Pass the binding-sphere radius via
            ``.add_sampler(sphere_radius=...)``.
        reversibility_check_every : int
            Cadence of the non-destructive HMC reversibility probe (THEORY 5.7);
            0 = OFF (default), N > 0 runs it every N rounds. See
            ``add_robotic_world`` for the full description.

        Returns
        -------
        World
            The new world, for chaining ``.add_sampler(...)``.

        Notes
        -----
        Root mobility here is a *world* property: it is built from
        ``ligand_molecule_indices`` and does NOT modify
        ``self.system_topology.root_mobilities``, so the same Context can host a
        docking world alongside ordinary Cartesian/torsional worlds.
        """
        if isinstance(ligand_molecule_indices, int):
            ligand_molecule_indices = [ligand_molecule_indices]
        world = super().add_docking_world([int(i) for i in ligand_molecule_indices])
        self._apply_mass_scale(world, mass_scale)
        if reversibility_check_every:
            world.set_reversibility_check(int(reversibility_check_every))
        return world

    def _apply_mass_scale(self, world, mass_scale) -> None:
        """Opt-in kinetic-metric mass scaling. Off (physical) when mass_scale is None.

        mass_scale may be:
          * None              -> off (default; bodies keep their physical inertia)
          * float             -> applied to all Free-root bodies (the canonical
                                 explicit-solvent / water-libration case)
          * {JointType: float}-> per-joint-type scale factors
        """
        if mass_scale is None:
            return
        items = (
            mass_scale.items()
            if isinstance(mass_scale, dict)
            else [(JointType.Free, float(mass_scale))]
        )
        for joint_type, scale in items:
            if scale is None or float(scale) == 1.0:
                continue  # 1.0 == physical == off
            world.set_mass_scale_by_joint(joint_type, float(scale))

    def add_robotic_world(
        self, selection, mass_scale=None, reversibility_check_every=0
    ):
        """Add a torsional (internal-coordinate) world.

        Parameters
        ----------
        selection : str
            Atom selection that becomes mobile (its torsions are sampled).
        mass_scale : None | float | dict[JointType, float]
            Opt-in kinetic-metric mass scaling (off / physical when None).
        reversibility_check_every : int
            Cadence of the non-destructive HMC reversibility probe (THEORY 5.7).
            0 = OFF (default; zero overhead). N > 0 runs
            ``RobotEngine::checkReversibility`` every N rounds (round 0 included,
            so it also serves as a startup check), logging the relative
            round-trip residual and warning if dt is too large for the current
            geometry. It certifies dt only for the configuration it runs from;
            the always-on guard remains the per-step corrector throw.
        """
        world = super().add_robotic_world(selection)
        self._apply_mass_scale(world, mass_scale)
        if reversibility_check_every:
            world.set_reversibility_check(int(reversibility_check_every))
        return world

    def add_torsional_world(
        self, selection, mass_scale=None, reversibility_check_every=0
    ):
        # Alias kept in sync with add_robotic_world (mass_scale + reversibility_check_every).
        return self.add_robotic_world(
            selection,
            mass_scale=mass_scale,
            reversibility_check_every=reversibility_check_every,
        )

    def add_cartesian_world(self, mass_scale=None, reversibility_check_every=0):
        world = super().add_cartesian_world()
        self._apply_mass_scale(world, mass_scale)
        if reversibility_check_every:
            world.set_reversibility_check(int(reversibility_check_every))
        return world

    @staticmethod
    def _is_water_composition(beg, end, mass, atomic_number) -> bool:
        """True iff atoms [beg, end) are one O + two H (+ massless extra points).

        ``mass`` / ``atomic_number`` are the WHOLE-system per-atom arrays (fetch
        them ONCE by the caller -- each ``system_topology.*`` access copies the
        full vector, so indexing the property inside a loop is O(N^2)).
        """
        n_o = n_h = n_other_massive = 0
        for a in range(beg, end):
            if float(mass[a]) <= 0.0:
                continue  # massless virtual site / extra point (TIP4P/TIP5P EP)
            z = int(atomic_number[a])
            if z == 8:
                n_o += 1
            elif z == 1:
                n_h += 1
            else:
                n_other_massive += 1
        return n_o == 1 and n_h == 2 and n_other_massive == 0

    def is_solvent(self, molecule_index: int) -> bool:
        """Return True iff molecule ``molecule_index`` is a water molecule.

        Assumes AMBER TIP*/SPC water models. A molecule is classified as water
        when its *massive* atoms are exactly one oxygen and two hydrogens; any
        remaining atoms must be massless virtual sites / extra points (the EP of
        TIP4P, TIP4P-Ew, OPC, TIP5P). Detection is by composition, not residue
        name, so it is robust to tleap naming (WAT / HOH / ...) and deliberately
        matches ONLY water -- ions and other small molecules (even O/H species
        like hydroxide or hydronium) are reported as non-solvent (solute).
        """
        st = self.system_topology
        return self._is_water_composition(
            int(st.atoms_begin[molecule_index]),
            int(st.atoms_end[molecule_index]),
            st.atoms_mass,
            st.atoms_atomic_number,
        )

    def add_ncmc_world(
        self,
        selection,
        timestep,
        ncmc_steps,
        hold_fraction=0.0,
        accept_reject_mode=None,
        use_fixman=None,
        mass_scale=None,
        ncmc_teleport=False,
        relax_solvent=False,
    ):
        """Add an NCMC torsional world with ALL solvent welded to Ground.

        Solute vs solvent
        -----------------
        Every non-water molecule is a *solute*. ``molecule_index`` is computed
        here as the list of all solute (non-solvent) molecule indices (water is
        detected by :meth:`is_solvent`). During each move the intermolecular
        nonbonded between the solute atoms and every other molecule is
        alchemically softened over a lambda:1->0->1 switch, so the solute strides
        through the cage of contacting molecules near lambda=0 and is recoupled
        with the accumulated work paid in the acceptance. Intramolecular physics
        is untouched. The solute atoms must form one contiguous block (true for a
        tleap-built system: solute first, then solvent); this is checked.

        Solvent is WELDed (root mobility -- a per-world property)
        --------------------------------------------------------
        Root mobility is set on THIS world only (via ``world.set_root_mobilities``;
        the shared topology is untouched). In this block **ALL solvent is welded to
        Ground** (Rigid root, zero DOF), exactly as the solutes' roots are: only
        the solute's internal (torsional) DOF move, plus the alchemical lambda
        stride-through. The solvent does NOT relax inside this move and it does not
        need to -- solvent relaxation and overall ergodicity are owned by the
        SEPARATE full-atom MD Gibbs block (composition of pi-preserving kernels);
        alchemy alone provides the stride-through of the welded environment (incl.
        docking). Welding all solvent is both simpler and removes a latent bias: a
        freed rigid-quaternion shell water would need its own orientation-Jacobian
        (logSineSqr) term, which historically was applied to molecule 0 only.

        Explicit solvent is supported: the alchemical decoupling has a PME-exact
        path (OpenMMContext::createAlchemyDecouplingForces). The one unsupported
        case is an NBFIX-corrected force field, which raises at initialize() by
        design. ``ncmc_steps`` is the protocol length (more steps => smoother
        switch, higher acceptance, slower); ``hold_fraction`` holds lambda=0 for
        that fraction of the protocol (the uncaged stride).

        Acceptance / the energy-pump caveat
        -----------------------------------
        The internal-coordinate integrator is only approximately symplectic, so
        long/large-dt trajectories can PUMP energy: every step then climbs in H and
        Metropolis rejects it, and the trajectory barely moves. The per-move
        ``[ncmc] ... dH=`` value is the acceptance driver (acc ~ exp(-dH/RT)); the
        ``gap=`` (work - dH) is diagnostic only. Levers, in rough order of effect:
        a smaller ``timestep``, fewer ``ncmc_steps``, and ``mass_scale`` (a
        kinetic-metric scaling that raises the stable dt ~sqrt(scale) with NO
        configurational bias -- Fixman cancels it). ``mass_scale=None`` (default)
        is physical / off.
        """
        st = self.system_topology
        num_mol = int(st.num_molecules)

        # Fetch the per-atom/per-molecule arrays ONCE (each property access
        # copies the whole vector -- never index the property inside a loop).
        atoms_begin = list(st.atoms_begin)
        atoms_end = list(st.atoms_end)
        mass = list(st.atoms_mass)
        atomic_number = list(st.atoms_atomic_number)

        # Partition molecules. molecule_index = all solute (non-water) indices.
        solvent_flags = [
            self._is_water_composition(
                atoms_begin[m], atoms_end[m], mass, atomic_number
            )
            for m in range(num_mol)
        ]
        molecule_index = [m for m in range(num_mol) if not solvent_flags[m]]
        if not molecule_index:
            raise ValueError(
                "add_ncmc_world: no solute (non-water) molecules found; nothing to sample."
            )

        # NCMC decouples the solute atoms from everything else over the lambda
        # switch. configure_ncmc takes a single contiguous [begin, end) range, so
        # the solutes must occupy a contiguous atom block; fail loudly otherwise.
        for prev, cur in zip(molecule_index, molecule_index[1:]):
            if atoms_end[prev] != atoms_begin[cur]:
                raise ValueError(
                    "add_ncmc_world: solute molecules are not contiguous in atom order "
                    f"(molecule {prev} ends at atom {atoms_end[prev]} but molecule {cur} "
                    f"begins at atom {atoms_begin[cur]}). The NCMC alchemical range "
                    "requires a single contiguous solute block."
                )
        beg = atoms_begin[molecule_index[0]]
        end = atoms_end[molecule_index[-1]]

        # Weld EVERYTHING's root to Ground: all solvent AND all solutes get a Rigid
        # (zero-DOF) root. Only the solute's internal (torsional) DOF -- supplied by
        # `selection` -- move, plus the alchemical lambda stride-through. The solvent
        # is deliberately frozen here; its relaxation + the chain's ergodicity are
        # owned by the SEPARATE full-atom MD Gibbs block (composition of
        # pi-preserving kernels). Welding all solvent also removes a latent
        # orientation-Jacobian (logSineSqr) bias that a freed rigid-quaternion shell
        # water would otherwise require.
        root_mob = [JointType.Rigid] * num_mol

        world = super().add_robotic_world(selection)

        # Root mobility is a per-world property: apply it to THIS world BEFORE
        # add_sampler (the rebuild resets per-body sampler state).
        world.set_root_mobilities(root_mob)

        # Mass scaling MUST come AFTER set_root_mobilities: that rebuild resets
        # every body's mass scale to 1.0. A kinetic-metric mass scale raises the
        # stable dt ~sqrt(scale) with ZERO configurational bias (Fixman cancels it
        # exactly). None (default) = physical/off.
        self._apply_mass_scale(world, mass_scale)

        world.add_sampler(
            timeStep=timestep,
            mdSteps=1,  # NCMC drives its own protocol loop, not mdSteps
            acceptRejectMode=(
                accept_reject_mode
                if accept_reject_mode is not None
                else AcceptRejectMode.MetropolisHastings
            ),
            use_nuts=False,
            use_fixman=use_fixman,
        )
        # configure_ncmc sets moveType=NcmcSwitch LAST (after add_sampler).
        world.configure_ncmc(beg, end, int(ncmc_steps), float(hold_fraction))
        # The lambda=0 trough teleport (default off, so omitting it reproduces prior
        # behaviour).
        if ncmc_teleport:
            world.set_ncmc_teleport(True)

        # Solvent-relaxing NCMC (docs/specs/ncmc_solvent_relax.md). With
        # relax_solvent=True the solvent stays WELDED as rigid bodies (no Fixman/
        # Jacobian contribution, exactly as above) but its atoms are additionally
        # advanced in FLAT Cartesian space by velocity-Verlet driven by OpenMM
        # forces INSIDE the proposal -- so the cage relaxes during the lambda
        # stride instead of being an infinite-mass wall, which is what kills
        # acceptance in explicit solvent / dense contact. Off (default) ==
        # welded-frozen solvent, bit-for-bit. Massless EPs are filtered out C++-side.
        if relax_solvent:
            solvent_atoms = []
            for m in range(num_mol):
                if solvent_flags[m]:
                    solvent_atoms.extend(range(atoms_begin[m], atoms_end[m]))
            world.set_cartesian_solvent(solvent_atoms)
        return world
