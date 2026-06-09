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
from .robo_bindings import SystemTopology
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
    )  # (n_atoms,) bool   -- atom is outboard across a PIN joint
    shake_pairs: list[tuple[int, int]] = field(default_factory=list)


class Context(_Context):
    system_topology: SystemTopology
    dihedral_classifier: AmberDihedralClassifier
    df_bonds: pd.DataFrame

    def __init__(self, dihedral_classifier: AmberDihedralClassifier) -> None:
        super().__init__()
        self.dihedral_classifier = dihedral_classifier

    def load_amber(
        self, prmtop_path: str | os.PathLike[str], inpcrd_path: str | os.PathLike[str]
    ) -> None:
        """
        Reads an AMBER parameter/topology file (`.prmtop`) and coordinate file (`.inpcrd` or `.rst7`) and creates a `SystemTopology` object that can be used to initialize a `Context`.

        Parameters
        ----------
        prmtop_path : str or os.PathLike
            Path to the AMBER parameter/topology file (`.prmtop`).
        inpcrd_path : str or os.PathLike
            Path to the AMBER coordinate file (`.inpcrd` or `.rst7`).

        Returns
        -------
        None
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

        # parm.split() groups molecules by topology identity.
        # Each entry is a 2-tuple:
        #
        #   (prototype_structure, instance_indices)
        #
        #   prototype_structure : parmed.Structure
        #       One representative copy of this molecule type, used to parse
        #       atoms, bonds, angles, and torsions.
        #
        #   instance_indices : list[int]
        #       Positions of every occurrence of this molecule type in the full
        #       system.  For example, [0, 3, 5] means copies 0, 3, and 5 in the
        #       topology belong to this type.
        #
        # Example for a box of 2 protein chains + 150 water molecules:
        #
        #   parm_prototypes[0] -> (protein_struct, [0, 1])
        #   parm_prototypes[1] -> (water_struct,   [2, 3, ..., 151])
        parm_prototypes: list[tuple[pmd.amber.AmberParm, list[int]]] = parm.split()

        # Build one MoleculePrototype per unique type.  Parsing bonds/angles/
        # torsions is expensive, so we do it once here and reuse across all
        # instances of the same type.
        molecule_prototypes = [
            MoleculePrototype(mol_struct, self.dihedral_classifier)
            for mol_struct, _ in parm_prototypes
        ]

        # Flatten into a (instance_index, prototype_index) list, then sort by
        # instance_index so the order matches the original topology.
        molecules: list[tuple[int, int]] = [
            (instance_index, prototype_index)
            for prototype_index, (_, instance_indices) in enumerate(parm_prototypes)
            for instance_index in instance_indices
        ]
        molecules.sort(key=lambda x: x[0])

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

        # ---- Local accumulators -------------------------------------------------
        # One plain Python list per output attribute.  Everything is built here,
        # then pushed to the topology at the end via setattr.
        acc: dict[str, list] = {}
        for spec in topology._RANGE_SPECS:
            acc[spec.begin] = []
            acc[spec.end] = []
        for fspec in topology._FIELD_SPECS:
            acc.setdefault(fspec.sys_attr, [])

        counters: dict[str, int] = {spec.counter: 0 for spec in topology._RANGE_SPECS}

        for instance_idx, prototype_idx in molecules:
            proto: MoleculePrototype = molecule_prototypes[prototype_idx]
            atom_off: int = counters["num_atoms"]

            # 1. Begin markers
            for spec in topology._RANGE_SPECS:
                acc[spec.begin].append(counters[spec.counter])

            # 2. Field data
            for fspec in topology._FIELD_SPECS:
                src: list = getattr(proto, fspec.proto_attr)
                dst: list = acc[fspec.sys_attr]

                if not fspec.atom_offset:
                    dst.extend(src)
                else:
                    dst.extend(x + atom_off for x in src)

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

        # ---- Flush accumulators to topology -------------------------------------
        # pybind11 def_readwrite uses copy-in / copy-out semantics: the getter
        # returns a Python copy of the C++ vector, not a reference.  Any
        # modification to that copy (append, extend, +=) is silently discarded.
        # The only reliable write path is a direct assignment that triggers
        # __setattr__, i.e. ``topology.attr = value`` or equivalently setattr().
        # Since SystemTopology() default-constructs with all vectors empty,
        # we can set each attribute directly to the accumulated list.
        for attr, data in acc.items():
            setattr(self.system_topology, attr, data)
        for counter, value in counters.items():
            setattr(self.system_topology, counter, value)
        self.system_topology.num_molecules = len(molecules)

        cmap_data = prmtop_reader.load_cmap(parm_file["raw_data"], parm)
        for attr, val in cmap_data.items():
            setattr(self.system_topology, attr, val)

        return
