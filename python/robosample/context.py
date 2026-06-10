from __future__ import annotations

import os
from collections import defaultdict
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

    def load_amber(
        self, prmtop_path: str | os.PathLike[str], inpcrd_path: str | os.PathLike[str]
    ) -> None:
        """
        Read an AMBER ``.prmtop`` + coordinate file and populate
        ``self.system_topology``.
        """
        self.system_topology = SystemTopology()

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

        # -----------------------------------------------------------------
        # Synthetic atom classes (transitional).
        # Computed in their own pass over the prototypes -- intentionally NOT
        # merged into the flattening loop below.
        # -----------------------------------------------------------------
        proto_class_names = self._generate_synthetic_atom_classes(molecule_prototypes)

        # Global atom-class index map (1-based, matching the historical i + 1).
        unique_atom_classes = sorted({n for names in proto_class_names for n in names})
        self.atom_class_indices = {
            name: i + 1 for i, name in enumerate(unique_atom_classes)
        }

        # Charged atom type = class name + partial charge (0-based index).
        proto_charged_keys: list[list[str]] = [
            [f"{names[c]}:{proto.atoms_charge[c]}" for c in range(proto.num_atoms)]
            for proto, names in zip(molecule_prototypes, proto_class_names)
        ]
        unique_charged = sorted({k for keys in proto_charged_keys for k in keys})
        self.charged_atom_type_indices = {
            key: i for i, key in enumerate(unique_charged)
        }

        # Stamp per-prototype, compound-ordered atom-class arrays for the loop
        # to emit: both the integer table indices and the parallel name strings.
        # (atoms_compound_atom_index is intrinsic to the prototype and provided
        # by MoleculePrototype directly.)
        for proto, names, keys in zip(
            molecule_prototypes, proto_class_names, proto_charged_keys
        ):
            proto.atoms_class_index = [self.atom_class_indices[n] for n in names]
            proto.atoms_charged_atom_type_index = [
                self.charged_atom_type_indices[k] for k in keys
            ]
            proto.atoms_class_names = list(names)
            proto.atoms_charged_type_names = list(keys)

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
        default_root_mobility = RootMobility.FREE

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
                acc["atoms_unique_name"].append(
                    f"{a.residue.name}{residue_off + a.residue.idx + 1}"
                    f"_{a.name}_{atom_off + local_idx + 1}"
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

        cmap_data = prmtop_reader.load_cmap(parm_file["raw_data"], parm)
        for attr, val in cmap_data.items():
            setattr(self.system_topology, attr, val)

        return

    def _generate_synthetic_atom_classes(
        self, molecule_prototypes: list[MoleculePrototype]
    ) -> list[list[str]]:
        """
        Assign a synthetic atom-class name to every atom of every prototype.

        Transitional: this will be removed once atom classes are derived
        elsewhere.  It walks each prototype's ParmEd molecule
        (``proto.molecule``), so it must be given the per-type prototypes from
        ``parm.split()``; its loop is deliberately kept separate from the
        flattening loop in :meth:`load_amber`.

        A class is an atom's unique "parameter environment": its atom type,
        mass, and the identities (``id()``) of every bond/angle/dihedral type
        incident on it.  Because each prototype is a distinct ParmEd structure
        with its own parameter-type objects, identical environments in
        *different* prototypes get distinct classes.  That over-splitting is
        safe -- more classes only means more (redundant) definitions, never an
        incorrect merge -- and is acceptable for a soon-to-be-removed helper.

        Returns
        -------
        list[list[str]]
            ``names[p][c]`` is the class name of compound atom ``c`` of
            prototype ``p`` (compound order, aligned with the prototype's
            per-atom arrays).
        """
        # ---- Capture each atom's parameter environment ---------------------
        signatures = defaultdict(
            lambda: {
                "type": None,
                "mass": None,
                "bonds": [],
                "angles": [],
                "dihedrals": [],
            }
        )

        for proto in molecule_prototypes:
            mol = proto.molecule

            for atom in mol.atoms:
                signatures[atom]["type"] = atom.type
                signatures[atom]["mass"] = atom.mass

            for bond in mol.bonds:
                bt_id = id(bond.type)
                signatures[bond.atom1]["bonds"].append(bt_id)
                signatures[bond.atom2]["bonds"].append(bt_id)

            for angle in mol.angles:
                at_id = id(angle.type)
                signatures[angle.atom1]["angles"].append(at_id)
                signatures[angle.atom2]["angles"].append(at_id)
                signatures[angle.atom3]["angles"].append(at_id)

            for dihed in mol.dihedrals:
                dt_id = id(dihed.type)
                signatures[dihed.atom1]["dihedrals"].append(dt_id)
                signatures[dihed.atom2]["dihedrals"].append(dt_id)
                signatures[dihed.atom3]["dihedrals"].append(dt_id)
                signatures[dihed.atom4]["dihedrals"].append(dt_id)

        # ---- Collapse to unique class names (compound order per prototype) --
        unique_sig_to_name: dict[tuple, str] = {}
        type_counters: dict[str, int] = defaultdict(int)
        proto_class_names: list[list[str]] = []

        for proto in molecule_prototypes:
            names: list[str] = []
            for local_index in proto.compound_to_local:  # compound order
                atom = proto.molecule.atoms[local_index]
                sig = (
                    signatures[atom]["type"],
                    signatures[atom]["mass"],
                    tuple(sorted(signatures[atom]["bonds"])),
                    tuple(sorted(signatures[atom]["angles"])),
                    tuple(sorted(signatures[atom]["dihedrals"])),
                )
                if sig not in unique_sig_to_name:
                    base_type = signatures[atom]["type"]
                    type_counters[base_type] += 1
                    unique_sig_to_name[sig] = f"{base_type}_{type_counters[base_type]}"
                names.append(unique_sig_to_name[sig])
            proto_class_names.append(names)

        return proto_class_names
