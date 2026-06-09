from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum, auto


class DSSPCode(IntEnum):
    """DSSP secondary structure assignment codes.

    Each code represents a local backbone conformation determined by
    hydrogen-bond pattern and backbone torsion angles. The string
    inheritance allows direct comparison with raw MDTraj output.

    Attributes:
        ALPHA_HELIX:      Right-handed alpha helix (i+4 H-bond).
        ISOLATED_BRIDGE:  Isolated beta-bridge (single inter-strand H-bond pair).
        EXTENDED_STRAND:  Beta-ladder strand (>=2 consecutive bridge residues).
        HELIX_310:        3/10 helix (i+3 H-bond).
        PI_HELIX:         Pi helix (i+5 H-bond); rare.
        TURN:             H-bonded turn (i+1 to i+4, not helix).
        BEND:             Local curvature > 70 deg; no H-bond criterion.
        LOOP:             Unassigned loop or irregular element.
        NOT_ASSIGNED:     Residue not assessed (non-protein or missing coords).
    """

    ALPHA_HELIX = auto()
    ISOLATED_BRIDGE = auto()
    EXTENDED_STRAND = auto()
    HELIX_310 = auto()
    PI_HELIX = auto()
    TURN = auto()
    BEND = auto()
    LOOP = auto()
    NOT_ASSIGNED = auto()

    @classmethod
    def from_mdtraj(cls, code: str) -> "DSSPCode":
        """Parse a single raw MDTraj DSSP string into a DSSPCode.

        Args:
            code: Single-character DSSP code or 'NA' as returned by
                  ``mdtraj.compute_dssp``.

        Returns:
            Corresponding DSSPCode member.

        Raises:
            ValueError: If *code* does not match any known assignment.
        """
        if code == "H":
            return cls.ALPHA_HELIX
        if code == "B":
            return cls.ISOLATED_BRIDGE
        if code == "E":
            return cls.EXTENDED_STRAND
        if code == "G":
            return cls.HELIX_310
        if code == "I":
            return cls.PI_HELIX
        if code == "T":
            return cls.TURN
        if code == "S":
            return cls.BEND
        if code == " ":
            return cls.LOOP
        if code == "NA":
            return cls.NOT_ASSIGNED

        raise ValueError(
            f"Unknown DSSP code {code!r}. "
            "Expected one of: 'H', 'B', 'E', 'G', 'I', 'T', 'S', ' ', 'NA'"
        )

    @property
    def is_helix(self) -> bool:
        """True for any helical conformation (H, G, I)."""
        return self in (self.ALPHA_HELIX, self.HELIX_310, self.PI_HELIX)

    @property
    def is_sheet(self) -> bool:
        """True for beta-strand conformations (B, E)."""
        return self in (self.ISOLATED_BRIDGE, self.EXTENDED_STRAND)

    @property
    def is_coil(self) -> bool:
        """True for non-regular structure (T, S, loop, or unassigned)."""
        return not (self.is_helix or self.is_sheet)

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return f"DSSPCode.{self.name}"


@dataclass(frozen=True)
class BondDSSP:
    """Per-bond secondary structure assignment derived from residue-level DSSP codes.

    DSSP assigns secondary structure per residue, but many analyses require
    a per-bond descriptor (e.g. coarse-grained network models, bond-order
    parameters, or visualisation pipelines). For intra-residue bonds the
    assignment is unambiguous. For inter-residue bonds -- most notably the
    peptide bond and disulfide bridges -- the two flanking residues may carry
    different codes, producing a boundary that must be resolved explicitly.

    Attributes:
        atom1_code: DSSP code of the residue containing the first atom.
        atom2_code: DSSP code of the residue containing the second atom.
    """

    atom1_code: DSSPCode
    atom2_code: DSSPCode

    @property
    def is_intra_residue(self) -> bool:
        """True when both atoms belong to residues with the same DSSP code.

        Note: this is a code-level equality check, not a residue index check.
        Two adjacent residues that happen to share the same code (e.g. two
        consecutive alpha-helix residues) will also return True.
        """
        return self.atom1_code == self.atom2_code

    @property
    def is_boundary(self) -> bool:
        """True when the bond spans a secondary structure boundary.

        Boundary bonds are structurally significant: the peptide bond at a
        helix terminus or a disulfide bridge connecting a helix to a strand
        are canonical examples. See ``resolve`` for strategies to collapse
        the two codes to a single assignment.
        """
        return self.atom1_code != self.atom2_code

    def resolve(self, strategy: str = "atom1") -> DSSPCode:
        """Collapse the two per-residue codes to a single bond-level code.

        For intra-residue bonds the result is unambiguous regardless of
        strategy. For boundary bonds the choice of strategy encodes a
        physical assumption about which residue's conformation dominates
        the bond's character.

        Args:
            strategy: Resolution rule applied at boundaries. One of:

                ``atom1``
                    Return the code of the first atom's residue. Suitable
                    when atom ordering follows a chemical convention such as
                    donor-acceptor or N-to-C directionality.

                ``atom2``
                    Return the code of the second atom's residue. Symmetric
                    counterpart to ``atom1``.

                ``priority``
                    Return the more ordered code, preferring helix over
                    sheet over coil. Appropriate when the bond is considered
                    part of the more structured region, e.g. a peptide bond
                    at a helix C-terminus retains helical character.

                ``none``
                    Return ``DSSPCode.NOT_ASSIGNED`` for any boundary bond.
                    Use when downstream code requires unambiguous single-
                    residue assignments and boundary cases must be excluded.

        Returns:
            A single ``DSSPCode`` representing the bond's secondary structure.

        Raises:
            ValueError: If *strategy* is not one of the four recognised values.
        """
        if not self.is_boundary:
            return self.atom1_code

        match strategy:
            case "atom1":
                return self.atom1_code
            case "atom2":
                return self.atom2_code
            case "priority":
                for code in (self.atom1_code, self.atom2_code):
                    if code.is_helix:
                        return code
                for code in (self.atom1_code, self.atom2_code):
                    if code.is_sheet:
                        return code
                return self.atom1_code
            case "none":
                return DSSPCode.NOT_ASSIGNED
            case _:
                raise ValueError(f"Unknown strategy {strategy!r}")
