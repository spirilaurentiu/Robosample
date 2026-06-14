from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import parmed as pmd

from . import prmtop_reader, topology
from .amber_dihedral_classifier import AmberDihedralClassifier
from .amber_dihedral_types import DihedralType
from .molecule_prototype import MoleculePrototype
from .robo_bindings import Context as _Context
from .robo_bindings import RootMobility, SystemTopology
from .secondary_structure import DSSPCode
from .units import ANG_TO_NM


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
            sele = ctx.build_flexibilities(bonds, rb.BondMobility.Torsion, False)
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
        self, prmtop_path: str | os.PathLike[str], inpcrd_path: str | os.PathLike[str]
    ) -> None:
        """
        Read an AMBER ``.prmtop`` + coordinate file and populate
        ``self.system_topology``.
        """
        self.system_topology = SystemTopology()
        self.system_topology.use_gbsa_obc2 = True
        self.system_topology.gbsa_solute_dielectric = 1.0
        self.system_topology.gbsa_solvent_dielectric = 78.5

        # Parse files and do typing
        parm_file = prmtop_reader.parse_prmtop(prmtop_path)
        parm: pmd.amber.AmberParm = pmd.load_file(prmtop_path, xyz=inpcrd_path)

        # Look at NB fix
        self.system_topology.num_nb_types = parm.pointers["NTYPES"]
        self.system_topology.a_coef, self.system_topology.b_coef = (
            prmtop_reader.load_lj_coefs(
                parm_file["raw_data"],
                self.system_topology.num_nb_types,
            )
        )

        # parm.split() groups molecules by topology identity; one prototype per
        # unique molecule type, plus the instance indices of every occurrence.
        parm_prototypes: list[tuple[pmd.amber.AmberParm, list[int]]] = parm.split()

        # Build one MoleculePrototype per unique type (parsing is expensive).
        molecule_prototypes = [
            MoleculePrototype(mol_struct, self.dihedral_classifier)
            for mol_struct, _ in parm_prototypes
        ]

        # Flatten into (instance_index, prototype_index), ordered by instance.
        self.molecules: list[tuple[int, int]] = [
            (instance_index, prototype_index)
            for prototype_index, (_, instance_indices) in enumerate(parm_prototypes)
            for instance_index in instance_indices
        ]
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
        self.df_bonds = pd.DataFrame(columns=COLUMNS)

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
        # prototype.  residue_off tracks the cumulative residue count in prmtop
        # (instance) order, mirroring how counters["num_atoms"] tracks atoms.
        acc.setdefault("atoms_unique_name", [])
        acc.setdefault("atoms_x", [])
        acc.setdefault("atoms_y", [])
        acc.setdefault("atoms_z", [])

        # BFS/compound position -> 0-based prmtop atom index. Appended in c
        # (compound) order below, so the list position IS the global BFS index;
        # transmitted to C++ to write trajectory output in prmtop atom order.
        acc.setdefault("atoms_prmtop_index", [])

        residue_off: int = 0

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
        default_root_mobility = RootMobility.WELD

        # prmtop_to_global_index maps each 0-based prmtop atom index to its
        # position in the BFS/compound-ordered global flat arrays.
        # Built incrementally inside the per-atom loop below.
        # Used after the loop by the CMAP section to convert 1-based prmtop
        # atom indices in CMAP_INDEX to global indices.
        self.prmtop_to_global_index: dict[int, int] = {}

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

            # 2c. Per-atom unique name, stored in compound order like every
            #     other atom array.  Format: "{resname}{res}_{atomname}_{atom}",
            #     e.g. "ALA1_N_4".  Both numbers are GLOBAL and 1-based:
            #       - res  = residue_off + (prototype-local residue idx) + 1
            #       - atom = atom_off    + (prototype-local prmtop idx)  + 1
            #     The trailing atom number is therefore the actual 1-based index
            #     in the prmtop file (the array *position* stays compound order).
            for c in range(proto.num_atoms):
                local_idx = proto.compound_to_local[c]
                a = proto.molecule.atoms[local_idx]

                # global index = position in the compound-ordered flat array
                global_idx: int = atom_off + c

                # prmtop index = atom_off + local (prmtop-order) index within molecule
                prmtop_idx: int = atom_off + local_idx
                self.prmtop_to_global_index[prmtop_idx] = global_idx

                # BFS->prmtop permutation for trajectory output (position in
                # this list == global_idx, value == prmtop_idx).
                acc["atoms_prmtop_index"].append(prmtop_idx)

                # Per-INSTANCE coordinates: read from the full parm at this
                # instance's global atom index (prmtop_idx), NOT from the shared
                # prototype (which would superimpose every copy of a molecule
                # type and blow the energy up to inf).
                inst_atom = parm.atoms[prmtop_idx]
                acc["atoms_x"].append(inst_atom.xx * ANG_TO_NM)
                acc["atoms_y"].append(inst_atom.xy * ANG_TO_NM)
                acc["atoms_z"].append(inst_atom.xz * ANG_TO_NM)
                acc["atoms_unique_name"].append(
                    f"{a.residue.name}{residue_off + a.residue.idx + 1}"
                    f"_{a.name}_{prmtop_idx + 1}"
                )

            # 3. Advance counters
            for spec in topology._RANGE_SPECS:
                counters[spec.counter] += getattr(proto, spec.proto_count)

            # 4. End markers
            for spec in topology._RANGE_SPECS:
                acc[spec.end].append(counters[spec.counter])

            # 5. Update bond DataFrame
            for i in range(proto.num_bonds):
                entry = {
                    "atom1_idx": [proto.bonds_i[i]],
                    "atom2_idx": [proto.bonds_j[i]],
                    "molecule_idx": [instance_idx],
                    "dihedral_type": [proto.bonds_dihedral_type[i]],
                    "is_ring_closing": [proto.bonds_is_ring_closing[i]],
                    "dss": [proto.bonds_secondary_structure[i]],
                    "residue_name": ["UNK"],
                    "residue_idx": [-1],
                }
                self.df_bonds = pd.concat(
                    [self.df_bonds, pd.DataFrame(entry)],
                    ignore_index=True,
                )

            # Advance the global residue offset for the next instance.
            residue_off += proto.num_residues

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
        raw_data = parm_file["raw_data"]

        cmap_resolution: list[int] = raw_data.get("CMAP_RESOLUTION", [])
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

            key = f"CMAP_PARAMETER_{grid_idx + 1:02d}"
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
        cmap_index_raw: list[int] = raw_data.get("CMAP_INDEX", [])
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
            return self.prmtop_to_global_index[prmtop_1based - 1]

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

        return

    def add_docking_world(self, ligand_molecule_indices):
        """Add a rigid-body docking world.

        Parameters
        ----------
        ligand_molecule_indices : int | Iterable[int]
            Molecule index/indices treated as ligands. Their root bodies become
            Free (6 external DOF); every other molecule is welded to Ground and
            rigid. The binding-site centre is the centroid of all non-ligand
            (receptor) atoms. Pass the binding-sphere radius via
            ``.add_sampler(sphere_radius=...)``.

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
        return super().add_docking_world([int(i) for i in ligand_molecule_indices])
