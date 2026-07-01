"""PyMOL Flex Editor plugin (FlexWizard).

Interactive authoring of Robosample ``.flex`` files from within PyMOL.

Usage
-----
1. In PyMOL: ``run /path/to/FlexWizard.py`` (File -> Run, or the command line).
2. At the PyMOL prompt: ``start_flex_wizard()``.

What it does
------------
* Walks the enabled object's bond graph and starts every bond as ``Rigid``.
* Lets you reassign bonds to a joint type (Torsion-focused for now) either one
  bond at a time (Manual Selection -> pick two bonded atoms) or in bulk
  ("Set Selection" over a PyMOL selection).
* Bonds are keyed by molecule-local (0-based) atom indices, matching how the
  Robosample topology stores ``bonds_i`` / ``bonds_j``.
* "Save Flexibility File" writes the ``.flex``; "Load Flexibility File" reads
  one back for further editing (round-trips with the writer).
* "Detect Dihedrals (AMBER)" generates bond metadata **on the fly** by building
  a throwaway Robosample context in-process from the topology (no file written).
  It pre-colors rigid / ring-closing bonds, locks ring-closing bonds against
  becoming torsions, enables one-click named dihedral selections (phi, psi,
  backbone, sidechain chi), and -- with "Apply to all copies" on -- mirrors each
  assignment onto every other copy of the same molecular prototype (e.g. all
  waters or repeated chains). Pass the topology up front with
  ``start_flex_wizard(prmtop=..., inpcrd=...)`` to run detection automatically.
  This requires Robosample to be importable in PyMOL's Python (run PyMOL from
  the Robosample conda env); otherwise the wizard stays in manual mode.
* An optional Joint Control Panel oscillates the selected torsions live so you
  can eyeball the motion each torsion unlocks.
"""

from __future__ import annotations

import math
import os

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets
from pymol.wizard import Wizard

# User-configurable joint types + colors.
JOINT_TYPES = ["Rigid", "Torsion", "BallF", "BallM", "Cartesian"]

JOINT_COLORS = {
    "Rigid": "0x696969",
    "Torsion": "0x006400",
    "BallM": "0x00BFFF",
    "BallF": "0x00FFBF",
    "Cartesian": "0xFFFF66",
}

JOINT_NAMES = {
    "Rigid": "Rigid Joint",
    "Torsion": "Torsion Joint",
    "BallF": "Ball Joint (Fixed Body)",
    "BallM": "Ball Joint (Mobile Body)",
    "Cartesian": "Cartesian Joint",
}

MANUAL_MODE_NAMES = {0: "Off", 1: "On"}

# Highlight for rotatable, still-unassigned bonds when a metadata sidecar is
# loaded (distinct from every JOINT_COLORS entry).
CANDIDATE_COLOR = "0xFFFFFF"

# DihedralType member names (see robosample.amber_dihedral_types) grouped into
# the named one-click selections the wizard offers. Kept as plain strings so the
# plugin needs no Robosample import.
DIHEDRAL_GROUPS = {
    "phi": {"PROTEIN_PHI"},
    "psi": {"PROTEIN_PSI"},
    "backbone": {"PROTEIN_PHI", "PROTEIN_PSI"},
    "sidechain": {
        "PROTEIN_CHI_1",
        "PROTEIN_CHI_2",
        "PROTEIN_CHI_3",
        "PROTEIN_CHI_4",
        "PROTEIN_CHI_5",
    },
}


class JointControlWindow(QtWidgets.QWidget):
    """Live oscillation of a set of torsions, for visual inspection."""

    def __init__(self, joints, parent=None):
        """Build the control panel.

        Parameters
        ----------
        joints : list[tuple[int, int, int, int]]
            Dihedral atom quadruples. Each atom index is a PyMOL atom ID
            (1-based source id), not a model index.
        """
        super().__init__(parent)
        self.setWindowTitle("Joint Control Panel")
        self.joints = joints
        n = len(joints)

        # Per-joint animation state.
        self.initial_angles = []
        self.sliders = []
        self.angles = [0.0 for _ in range(n)]
        self.centers = [0.0 for _ in range(n)]  # center value for rocking
        self.speeds = [0.0 for _ in range(n)]
        self.amps = [10.0 for _ in range(n)]  # default amplitude
        self.phases = [0.0 for _ in range(n)]  # internal phase offset
        self.speed_boxes = []  # one QDoubleSpinBox per joint
        self.amp_boxes = []

        # Outer layout.
        self.setLayout(QtWidgets.QVBoxLayout())

        # Scroll area holding one row per joint.
        scroll_area = QtWidgets.QScrollArea()
        scroll_area.setWidgetResizable(True)
        self.layout().addWidget(scroll_area)

        self.scroll_widget = QtWidgets.QWidget()
        self.scroll_layout = QtWidgets.QVBoxLayout()
        self.scroll_widget.setLayout(self.scroll_layout)
        scroll_area.setWidget(self.scroll_widget)

        for i, (a1, a2, a3, a4) in enumerate(joints):
            row = QtWidgets.QHBoxLayout()

            label = QtWidgets.QLabel(f"Joint {i + 1}: {a1}-{a2}-{a3}-{a4}")
            label.setFixedWidth(180)

            slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
            slider.setRange(-180, 180)
            initial_angle = cmd.get_dihedral(
                f"id {a1}", f"id {a2}", f"id {a3}", f"id {a4}"
            )
            self.initial_angles.append(initial_angle)
            self.centers[i] = initial_angle
            self.angles[i] = initial_angle
            slider.setValue(int(initial_angle))
            slider.setSingleStep(5)
            slider.setMinimumWidth(200)
            slider.valueChanged.connect(
                lambda val, idx=i: self._set_center(idx, val)
            )

            speed_box = QtWidgets.QDoubleSpinBox()
            speed_box.setRange(0.0, 30.0)
            speed_box.setSingleStep(0.5)
            speed_box.setValue(0.0)
            speed_box.setSuffix(" °/step")
            speed_box.setFixedWidth(90)
            speed_box.valueChanged.connect(
                lambda val, idx=i: self._set_speed(idx, val)
            )

            amp_box = QtWidgets.QDoubleSpinBox()
            amp_box.setRange(0.0, 180.0)
            amp_box.setSingleStep(5.0)
            amp_box.setValue(10.0)
            amp_box.setSuffix(" °")
            amp_box.setFixedWidth(80)
            amp_box.valueChanged.connect(
                lambda val, idx=i: self._set_amplitude(idx, val)
            )

            row.addWidget(label)
            row.addWidget(slider)
            row.addWidget(QtWidgets.QLabel("Speed:"))
            row.addWidget(speed_box)
            row.addWidget(QtWidgets.QLabel("Amp:"))
            row.addWidget(amp_box)
            self.scroll_layout.addLayout(row)

            self.speed_boxes.append(speed_box)
            self.amp_boxes.append(amp_box)
            self.sliders.append(slider)

        # Timer drives the rocking; started on demand via "Start All".
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self._update_angles)

        self.pause_button = QtWidgets.QPushButton("Start All")
        self.pause_button.clicked.connect(self.toggle_pause)
        self.layout().addWidget(self.pause_button)

        commit_button = QtWidgets.QPushButton("Commit Changes")
        commit_button.clicked.connect(self.commit_changes)
        self.layout().addWidget(commit_button)

        close_button = QtWidgets.QPushButton("Close")
        close_button.clicked.connect(self.close)
        self.layout().addWidget(close_button)

        # Batch controls.
        batch_layout = QtWidgets.QHBoxLayout()

        self.all_speed_button = QtWidgets.QPushButton("Set all speeds to")
        self.all_speed_button.clicked.connect(self.set_all_speeds)
        batch_layout.addWidget(self.all_speed_button)
        self.all_speed_spin = QtWidgets.QDoubleSpinBox()
        self.all_speed_spin.setRange(0.0, 30.0)
        self.all_speed_spin.setSingleStep(0.5)
        self.all_speed_spin.setValue(0.0)
        self.all_speed_spin.setSuffix(" °/step")
        batch_layout.addWidget(self.all_speed_spin)

        self.all_amp_button = QtWidgets.QPushButton("Set all amplitudes to")
        self.all_amp_button.clicked.connect(self.set_all_amps)
        batch_layout.addWidget(self.all_amp_button)
        self.all_amp_spin = QtWidgets.QDoubleSpinBox()
        self.all_amp_spin.setRange(0.0, 180.0)
        self.all_amp_spin.setSingleStep(5.0)
        self.all_amp_spin.setValue(10.0)
        self.all_amp_spin.setSuffix(" °")
        batch_layout.addWidget(self.all_amp_spin)

        self.layout().addLayout(batch_layout)

    def commit_changes(self):
        """Commit the current angles as the new reference (copy, not alias)."""
        self.initial_angles = list(self.angles)
        print("[FlexWizard] Current angles committed as new reference.")

    def closeEvent(self, event):
        """Stop the timer and restore the reference angles on close."""
        if hasattr(self, "timer") and self.timer.isActive():
            self.timer.stop()

        for idx, (a1, a2, a3, a4) in enumerate(self.joints):
            cmd.set_dihedral(
                f"id {a1}",
                f"id {a2}",
                f"id {a3}",
                f"id {a4}",
                self.initial_angles[idx],
            )
        print("[FlexWizard] Initial angles restored.")

        if hasattr(self.parent(), "joint_window"):
            self.parent().joint_window = None
        event.accept()

    def toggle_pause(self):
        """Toggle the rocking timer on/off."""
        if hasattr(self, "timer") and self.timer.isActive():
            self.timer.stop()
            self.pause_button.setText("Start All")
        else:
            self.timer.start(30)
            self.pause_button.setText("Stop All")

    def set_all_speeds(self):
        """Apply the batch speed value to every joint."""
        val = self.all_speed_spin.value()
        for i in range(len(self.speeds)):
            self._set_speed(i, val)
            self.speed_boxes[i].blockSignals(True)
            self.speed_boxes[i].setValue(val)
            self.speed_boxes[i].blockSignals(False)

    def set_all_amps(self):
        """Apply the batch amplitude value to every joint."""
        val = self.all_amp_spin.value()
        for i in range(len(self.amps)):
            self._set_amplitude(i, val)
            self.amp_boxes[i].blockSignals(True)
            self.amp_boxes[i].setValue(val)
            self.amp_boxes[i].blockSignals(False)

    # --- GUI event handlers ---

    def _set_center(self, idx, value):
        """User moved a slider manually -> update the rocking center."""
        self.centers[idx] = value
        self.angles[idx] = value
        self._apply_dihedral(idx)

    def _set_speed(self, idx, speed):
        """Speed > 0 enables rocking for this joint."""
        self.speeds[idx] = speed

    def _set_amplitude(self, idx, amp):
        """Set the oscillation amplitude for this joint."""
        self.amps[idx] = amp

    # --- Core update logic ---

    def _update_angles(self):
        """Oscillate each torsion around its center."""
        for idx in range(len(self.joints)):
            speed = self.speeds[idx]
            if speed <= 0.0:
                continue

            self.phases[idx] += speed * 0.05
            new_angle = self.centers[idx] + self.amps[idx] * math.sin(
                self.phases[idx]
            )

            self.angles[idx] = new_angle
            self.sliders[idx].blockSignals(True)
            self.sliders[idx].setValue(int(new_angle))
            self.sliders[idx].blockSignals(False)
            self._apply_dihedral(idx)

    def _apply_dihedral(self, idx):
        """Apply the current angle of joint ``idx`` to PyMOL."""
        a1, a2, a3, a4 = self.joints[idx]
        angle = self.angles[idx]
        cmd.set_dihedral(
            f"id {a1}", f"id {a2}", f"id {a3}", f"id {a4}", float(angle)
        )


class FlexManager:
    """Track flexible-bond assignments and provide toggle/save helpers."""

    def __init__(self):
        self.joints = {}

    def _normalize_bond(self, a1, a2):
        """Order atom indices so a bond is unique regardless of direction."""
        return tuple(sorted((int(a1), int(a2))))

    def set_joint(self, a1, a2, joint_data, bond_map):
        """Set the joint type of the bond between ``a1`` and ``a2``.

        ``joint_data`` of length 1 updates only the type; a longer list
        replaces the full record ``[type, name1, resn1, resi1, name2, resn2,
        resi2]``.

        Returns the dihedral quadruple (list) when a Torsion was accepted,
        ``False`` when a requested Torsion had no valid dihedral (and was
        demoted to Rigid), ``True`` for any other single-type update, and
        ``None`` when the full record was stored.
        """
        bond = self._normalize_bond(a1, a2)
        if len(joint_data) == 1:
            if joint_data[0] == "Torsion":
                dihedral = self.get_dihedral_from_bond(bond[0], bond[1], bond_map)
                if dihedral:
                    self.joints[bond][0] = "Torsion"
                    return dihedral
                self.joints[bond][0] = "Rigid"
                return False
            self.joints[bond][0] = joint_data[0]
            return True
        self.joints[bond] = joint_data
        return None

    def get_dihedral_from_bond(self, a1, a2, bond_map):
        """Return a dihedral ``[n1, a1, a2, n2]`` for the bond, or ``False``.

        Neighbors are chosen deterministically (lowest atom index) so the
        pinned/animated dihedral is stable across runs.
        """
        a1_neighbors = sorted(n for n in bond_map.get(a1, ()) if n != a2)
        a2_neighbors = sorted(n for n in bond_map.get(a2, ()) if n != a1)
        if not a1_neighbors or not a2_neighbors:
            return False
        return [a1_neighbors[0], a1, a2, a2_neighbors[0]]


class FlexWizard(Wizard):
    """Interactive bond-picking wizard.

    Pick two bonded atoms (Manual Selection) or a whole selection to assign the
    current joint type. Selections are cleared after each pair so picks chain.
    """

    def __init__(self, flex_manager, prmtop=None, inpcrd=None):
        print("[FlexWizard] Welcome to the Robosample Flexibility Wizard!")
        # Older PyMOL builds lack undo_disable; guard the call.
        try:
            cmd.undo_disable()
            print(
                "[FlexWizard] Undo has been disabled, as it interferes with "
                "wizard functionality."
            )
        except AttributeError:
            pass
        Wizard.__init__(self)
        # Clean up any leftover pick selections from a prior instance.
        cmd.delete("_pk1")
        cmd.delete("_pk2")

        self.joint_type = self.session.get("default_mode", "Rigid")
        self.manual_mode = 0

        self.flex_manager = flex_manager
        self.atom1 = None
        self.object_name = None
        self.bond_map = None
        self.pin_joints = []
        # (a, b) 0-based sorted -> {"dihedral_type", "is_ring_closing",
        # "is_rigid", "group"}; populated by detect_dihedrals, empty otherwise.
        self.bond_meta = {}
        # group id -> [(a, b), ...] equivalent bonds across prototype copies.
        self.bond_groups = {}
        # When on, an assignment propagates to every copy of the prototype bond.
        self.apply_to_all_copies = 0
        # AMBER files for on-the-fly dihedral detection (optional).
        self.prmtop_path = prmtop
        self.inpcrd_path = inpcrd

        # Build the bond map for the currently enabled object.
        self.build_bond_map()

        # If a topology was supplied, detect dihedrals immediately (no file).
        if self.prmtop_path:
            self.detect_dihedrals(auto=True)

        # Ensure PyMOL is not in edit mode.
        cmd.edit_mode("on")
        cmd.edit_mode("off")

        self.status = 0

        # Build the joint-type dropdown (see measurement.py for the format).
        smm = [[2, "Joint Type", ""]]
        for a in JOINT_TYPES:
            smm.append(
                [1, JOINT_NAMES[a], 'cmd.get_wizard().set_joint_type("' + a + '")']
            )
        self.menu["jointType"] = smm

    def build_bond_map(self):
        """Build the adjacency map for the enabled object and default to Rigid."""
        object_list = cmd.get_names("objects", enabled_only=1)
        if len(object_list) != 1:
            print(
                "[FlexWizard] Please enable an object for which to build the "
                "bond map."
            )
            self.bond_map = None
            return

        self.object_name = object_list[0]
        print(f"[FlexWizard] Built bond map for {self.object_name.upper()}.")
        model = cmd.get_model(f"model {self.object_name}")
        bond_map = {}

        # PyMOL bond indices are model ranks; convert them to atom IDs.
        rank_to_id = {i: atom.id for i, atom in enumerate(model.atom)}
        bond_ids = [
            (rank_to_id[a], rank_to_id[b])
            for (a, b) in (b.index for b in model.bond)
        ]

        for a1, a2 in bond_ids:
            # Store bonds by 0-based atom index (source id is 1-based).
            a1 -= 1
            a2 -= 1

            atom1 = model.atom[a1]
            atom2 = model.atom[a2]

            bond_map.setdefault(a1, set()).add(a2)
            bond_map.setdefault(a2, set()).add(a1)

            # Every bond starts Rigid unless the user changes it.
            self.flex_manager.set_joint(
                a1,
                a2,
                [
                    self.joint_type,
                    atom1.name,
                    atom1.resn,
                    atom1.resi,
                    atom2.name,
                    atom2.resn,
                    atom2.resi,
                ],
                None,
            )

        self.bond_map = bond_map

        # Reset the torsion pin list and any prior coloring.
        self.pin_joints = []
        cmd.unset_bond("stick_color", "all")
        cmd.unset_bond("line_color", "all")

    def _color_bond(self, sel1, sel2, color):
        """Color a bond (given two atom selections) in stick and line reps."""
        cmd.set_bond("stick_color", color, sel1, sel2)
        cmd.set_bond("line_color", color, sel1, sel2)

    def _apply_joint_type(self, atom1, atom2, joint_type, sel1, sel2, propagate=True):
        """Assign ``joint_type`` to a bond (0-based indices) and color it.

        Shared by manual pick, bulk selection, and file load so all three paths
        behave identically: a Torsion that cannot form a dihedral is demoted to
        Rigid and colored as Rigid, and accepted torsions are pinned once. When
        bond metadata is loaded, a Torsion on a ring-closing bond is refused
        (the engine forces such bonds Rigid) so the drawing matches the run.

        When ``propagate`` is true and "Apply to all copies" is on, the same
        assignment is mirrored to every other copy of this prototype bond.
        """
        if joint_type == "Torsion" and self._is_ring_closing(atom1, atom2):
            print(
                f"[FlexWizard] Bond {atom1}-{atom2} is ring-closing; kept Rigid "
                "(the engine cannot make it a torsion)."
            )
            self.flex_manager.set_joint(atom1, atom2, ["Rigid"], self.bond_map)
            self._color_bond(sel1, sel2, JOINT_COLORS["Rigid"])
            self._propagate_to_copies(atom1, atom2, joint_type, propagate)
            return
        result = self.flex_manager.set_joint(
            atom1, atom2, [joint_type], self.bond_map
        )
        if isinstance(result, list):
            pinned = [x + 1 for x in result]
            if pinned not in self.pin_joints:
                self.pin_joints.append(pinned)
            color = JOINT_COLORS[joint_type]
        elif result is False:
            # Requested Torsion had no valid dihedral -> demoted to Rigid.
            color = JOINT_COLORS["Rigid"]
        else:
            color = JOINT_COLORS[joint_type]
        self._color_bond(sel1, sel2, color)
        self._propagate_to_copies(atom1, atom2, joint_type, propagate)

    def _sibling_bonds(self, atom1, atom2):
        """Other bonds sharing this bond's prototype-copy group (excl. self)."""
        key = tuple(sorted((int(atom1), int(atom2))))
        meta = self.bond_meta.get(key)
        if not meta or meta.get("group") is None:
            return []
        return [b for b in self.bond_groups.get(meta["group"], []) if b != key]

    def _propagate_to_copies(self, atom1, atom2, joint_type, propagate):
        """Mirror an assignment to every copy of this prototype bond."""
        if not (propagate and self.apply_to_all_copies):
            return
        for sib1, sib2 in self._sibling_bonds(atom1, atom2):
            if (sib1, sib2) not in self.flex_manager.joints:
                continue
            sel1 = f"model {self.object_name} and id {sib1 + 1}"
            sel2 = f"model {self.object_name} and id {sib2 + 1}"
            self._apply_joint_type(
                sib1, sib2, joint_type, sel1, sel2, propagate=False
            )

    def toggle_apply_all(self):
        """Toggle whether assignments propagate to all prototype copies."""
        self.apply_to_all_copies = 0 if self.apply_to_all_copies else 1
        if self.apply_to_all_copies and not self.bond_groups:
            print(
                "[FlexWizard] Apply-to-all-copies is on, but no metadata is "
                "loaded yet (Load Bond Metadata to enable copy propagation)."
            )
        cmd.refresh_wizard()

    def do_select(self, name):
        if self.manual_mode == 0:
            return
        self.cmd.unpick()
        if self.status == 0:
            cmd.select("_pk1", name)
            cmd.delete(name)
            self.status = 1  # first atom picked
            self.prompt = ["Pick the second atom."]
            cmd.refresh_wizard()
        elif self.status == 1:
            cmd.select("_pk2", name)
            cmd.delete(name)
            self.do_pick()

    def do_pick(self):
        idx1 = cmd.id_atom("_pk1") - 1
        idx2 = cmd.id_atom("_pk2") - 1

        model01 = cmd.index("_pk1")[0][0]
        model02 = cmd.index("_pk2")[0][0]

        if idx1 == idx2:
            print("[FlexWizard] Selected the same atom twice.")
            self.restore()
            return

        if model01 != model02:
            print("[FlexWizard] Atoms from different models.")
            self.restore()
            return

        if idx2 not in self.bond_map.get(idx1, set()):
            print(f"[FlexWizard] Atoms {idx1} and {idx2} are not bonded.")
            self.restore()
            return

        self._apply_joint_type(idx1, idx2, self.joint_type, "_pk1", "_pk2")
        self.restore()

    def write_flex(self):
        """Write the current assignments to a ``.flex`` file."""
        object_list = cmd.get_names("objects", enabled_only=1)
        if object_list:
            self.object_name = object_list[0]

        dialog = QtWidgets.QFileDialog()
        filename, _ = dialog.getSaveFileName(
            None,
            "Save Flexibility File",
            self.object_name,
            "Flex files (*.flex);;All files (*)",
        )
        if not filename:
            return

        # Sort by first atom index so the file is easy to scan.
        data = dict(
            sorted(self.flex_manager.joints.items(), key=lambda item: item[0][0])
        )

        with open(filename, "w") as f:
            for (a1, a2), record in data.items():
                f.write(
                    "{:9s} {:9s} {:5s}\t#{:3s}_{}_{}-{:3s}_{}_{}\n".format(
                        str(a1),
                        str(a2),
                        record[0],
                        record[2],
                        record[1],
                        record[3],
                        record[5],
                        record[4],
                        record[6],
                    )
                )

        print(f"[FlexWizard] Flex file saved to: {filename}.")

    def load_flex(self):
        """Read a ``.flex`` file back in and re-apply its assignments.

        Round-trips with :meth:`write_flex`: the first two whitespace-separated
        fields are 0-based atom indices and the third is the joint type; the
        trailing ``#`` comment is ignored.
        """
        if self.bond_map is None:
            print(
                "[FlexWizard] No bond map built. Enable an object and build the "
                "bond map first."
            )
            return

        dialog = QtWidgets.QFileDialog()
        filename, _ = dialog.getOpenFileName(
            None,
            "Load Flexibility File",
            self.object_name or "",
            "Flex files (*.flex);;All files (*)",
        )
        if not filename:
            return

        applied = 0
        with open(filename) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) < 3:
                    continue
                try:
                    atom1, atom2 = int(fields[0]), int(fields[1])
                except ValueError:
                    continue
                joint_type = fields[2]
                if joint_type not in JOINT_TYPES:
                    print(f"[FlexWizard] Skipping unknown joint type {joint_type!r}.")
                    continue
                bond = tuple(sorted((atom1, atom2)))
                if bond not in self.flex_manager.joints:
                    print(
                        f"[FlexWizard] Skipping {atom1}-{atom2}: not a bond of "
                        f"{self.object_name}."
                    )
                    continue
                sel1 = f"model {self.object_name} and id {atom1 + 1}"
                sel2 = f"model {self.object_name} and id {atom2 + 1}"
                self._apply_joint_type(atom1, atom2, joint_type, sel1, sel2)
                applied += 1

        print(f"[FlexWizard] Loaded {applied} joint assignments from {filename}.")
        cmd.refresh_wizard()

    # --- AMBER dihedral metadata (see robosample.flex_export) ---

    def _is_ring_closing(self, atom1, atom2):
        """True if metadata marks this bond ring-closing (False if unknown)."""
        meta = self.bond_meta.get(tuple(sorted((int(atom1), int(atom2)))))
        return bool(meta and meta["is_ring_closing"])

    def set_topology(self, prmtop=None, inpcrd=None):
        """Record the AMBER files used for on-the-fly dihedral detection."""
        if prmtop is not None:
            self.prmtop_path = prmtop
        if inpcrd is not None:
            self.inpcrd_path = inpcrd

    def _resolve_inpcrd(self, prmtop):
        """Best-effort sibling coordinate file for ``prmtop`` (or None)."""
        if self.inpcrd_path:
            return self.inpcrd_path
        stem = os.path.splitext(str(prmtop))[0]
        for ext in (".rst7", ".inpcrd", ".rst", ".crd", ".ncrst"):
            candidate = stem + ext
            if os.path.exists(candidate):
                return candidate
        return None

    def detect_dihedrals(self, auto=False):
        """Generate AMBER bond metadata on the fly and ingest it (no file).

        Builds a throwaway Robosample context in-process from the topology, so
        nothing is written to disk. Requires Robosample to be importable in
        PyMOL's Python; if it is not, prints how to proceed and leaves the
        wizard in manual mode. When ``auto`` is true (startup path) a missing
        topology is silent rather than a prompt.
        """
        if self.bond_map is None:
            if not auto:
                print(
                    "[FlexWizard] No bond map built. Enable an object and build "
                    "the bond map first."
                )
            return

        prmtop = self.prmtop_path
        if not prmtop and not auto:
            dialog = QtWidgets.QFileDialog()
            prmtop, _ = dialog.getOpenFileName(
                None,
                "Select AMBER topology (.prmtop) for this object",
                self.object_name or "",
                "AMBER prmtop (*.prmtop *.parm7 *.top);;All files (*)",
            )
        if not prmtop:
            if not auto:
                print("[FlexWizard] No topology given; staying in manual mode.")
            return
        self.prmtop_path = prmtop

        inpcrd = self._resolve_inpcrd(prmtop)
        if not inpcrd:
            print(
                "[FlexWizard] Could not find a coordinate file next to "
                f"{prmtop}. Call set_topology(prmtop, inpcrd) with an explicit "
                "coordinate file."
            )
            return

        try:
            from robosample import flex_export
        except ImportError:
            print(
                "[FlexWizard] Robosample is not importable in this PyMOL. Run "
                "PyMOL from the Robosample conda env for on-the-fly detection, "
                "or precompute a sidecar with robosample.flex_export."
            )
            return

        try:
            rows, n_atoms = flex_export.compute_flex_meta_from_files(prmtop, inpcrd)
        except Exception as exc:  # noqa: BLE001 - surface any build failure
            print(f"[FlexWizard] Dihedral detection failed: {exc!r}")
            return

        self._ingest_meta_rows(rows, n_atoms)

    def _ingest_meta_rows(self, rows, declared_atoms):
        """Populate ``bond_meta`` / ``bond_groups`` from metadata rows.

        ``rows`` are ``(atom1, atom2, dihedral_type, is_ring_closing, is_rigid,
        group)`` tuples (0-based prmtop indices). Pre-colors bonds afterward.
        """
        meta = {}
        groups = {}
        for atom1, atom2, dihedral_type, ring, rigid, group in rows:
            key = tuple(sorted((int(atom1), int(atom2))))
            meta[key] = {
                "dihedral_type": dihedral_type,
                "is_ring_closing": bool(ring),
                "is_rigid": bool(rigid),
                "group": None if group is None else int(group),
            }
            if group is not None:
                groups.setdefault(int(group), []).append(key)

        # Guard the index bridge: prmtop indices only line up with PyMOL's atom
        # ids if the enabled object has the same atom order/count.
        object_atoms = len(cmd.get_model(f"model {self.object_name}").atom)
        if declared_atoms is not None and declared_atoms != object_atoms:
            print(
                f"[FlexWizard] WARNING: metadata is for {declared_atoms} atoms "
                f"but {self.object_name} has {object_atoms}. Indices may not "
                "align; the enabled object must be the same topology."
            )

        self.bond_meta = meta
        self.bond_groups = groups
        self._precolor_from_meta()
        n_multi = sum(1 for members in groups.values() if len(members) > 1)
        print(
            f"[FlexWizard] Detected metadata for {len(meta)} bonds "
            f"({n_multi} multi-copy group(s))."
        )
        cmd.refresh_wizard()

    def _precolor_from_meta(self):
        """Tint bonds by their metadata: rigid/ring-closing vs candidate."""
        for (atom1, atom2), info in self.bond_meta.items():
            if tuple(sorted((atom1, atom2))) not in self.flex_manager.joints:
                continue
            sel1 = f"model {self.object_name} and id {atom1 + 1}"
            sel2 = f"model {self.object_name} and id {atom2 + 1}"
            if info["is_ring_closing"] or info["is_rigid"]:
                self._color_bond(sel1, sel2, JOINT_COLORS["Rigid"])
            else:
                self._color_bond(sel1, sel2, CANDIDATE_COLOR)

    def select_dihedral_type(self, group):
        """Set every bond whose metadata type is in ``group`` to Torsion.

        ``group`` is a key of :data:`DIHEDRAL_GROUPS` (phi, psi, backbone,
        sidechain). Ring-closing bonds are refused by :meth:`_apply_joint_type`.
        """
        if not self.bond_meta:
            print("[FlexWizard] Load bond metadata first (Load Bond Metadata).")
            return
        type_names = DIHEDRAL_GROUPS[group]
        count = 0
        for (atom1, atom2), info in self.bond_meta.items():
            if info["dihedral_type"] not in type_names:
                continue
            if tuple(sorted((atom1, atom2))) not in self.flex_manager.joints:
                continue
            sel1 = f"model {self.object_name} and id {atom1 + 1}"
            sel2 = f"model {self.object_name} and id {atom2 + 1}"
            self._apply_joint_type(atom1, atom2, "Torsion", sel1, sel2)
            count += 1
        print(f"[FlexWizard] Set {count} {group} bond(s) to Torsion.")
        cmd.refresh_wizard()

    def restore(self):
        """Return the wizard to its ready-to-pick state."""
        cmd.delete("_pk1")
        cmd.delete("_pk2")
        self.status = 0
        self.toggle_manual_mode()
        self.toggle_manual_mode()
        cmd.refresh_wizard()
        cmd.unpick()

    def set_joint_type(self, joint_type):
        """Set the active joint type used by subsequent picks."""
        if joint_type in JOINT_TYPES:
            self.joint_type = joint_type
        self.restore()

    def toggle_manual_mode(self):
        """Toggle between atomic picking (on) and residue picking (off)."""
        if self.manual_mode == 0:
            self.manual_mode = 1
            cmd.set("mouse_selection_mode", 0)  # 0 = atomic picking
            self.prompt = ["Pick first atom to toggle a bond."]
        else:
            self.manual_mode = 0
            cmd.set("mouse_selection_mode", 1)  # 1 = residue picking
            self.prompt = None
        cmd.refresh_wizard()

    # --- Panel functions ---

    def get_panel(self):
        """Define the PyMOL Wizard control-panel layout."""
        return [
            [1, "FlexWizard - by T. A. Sulea.", "", ""],
            [3, "Joint type: " + JOINT_NAMES[self.joint_type], "jointType"],
            [
                2,
                "Manual Selection " + MANUAL_MODE_NAMES[self.manual_mode],
                "cmd.get_wizard().toggle_manual_mode()",
                "",
            ],
            [2, "Rebuild Bond Map", "cmd.get_wizard().build_bond_map()", ""],
            [2, "Save Flexibility File", "cmd.get_wizard().write_flex()", ""],
            [2, "Load Flexibility File", "cmd.get_wizard().load_flex()", ""],
            [2, "Detect Dihedrals (AMBER)", "cmd.get_wizard().detect_dihedrals()", ""],
            [
                2,
                "Apply to all copies " + MANUAL_MODE_NAMES[self.apply_to_all_copies],
                "cmd.get_wizard().toggle_apply_all()",
                "",
            ],
            [
                2,
                "phi -> Torsion",
                'cmd.get_wizard().select_dihedral_type("phi")',
                "",
            ],
            [
                2,
                "psi -> Torsion",
                'cmd.get_wizard().select_dihedral_type("psi")',
                "",
            ],
            [
                2,
                "Backbone -> Torsion",
                'cmd.get_wizard().select_dihedral_type("backbone")',
                "",
            ],
            [
                2,
                "Sidechain (chi) -> Torsion",
                'cmd.get_wizard().select_dihedral_type("sidechain")',
                "",
            ],
            [2, "Set Selection", "cmd.get_wizard().set_selection()", ""],
            [
                2,
                "Open Joint Control",
                "cmd.get_wizard().open_joint_control_panel()",
                "",
            ],
            [2, "Help Me!", "cmd.get_wizard().show_help()", ""],
            [2, "Done", "stop_flex_wizard()", ""],
        ]

    def open_joint_control_panel(self):
        """Open the live oscillation panel for the pinned torsions."""
        joints = self.pin_joints  # list of [a1, a2, a3, a4] (1-based ids)
        print(joints)
        if len(joints) == 0:
            print("Please select at least one Torsion-type joint.")
            return

        self.joint_window = JointControlWindow(joints)
        self.joint_window.show()

    def set_selection(self):
        """Assign the current joint type to every bond within a selection."""
        if self.bond_map is None:
            print(
                "[FlexWizard] No bond map built. Enable an object and build the "
                "bond map first."
            )
            return

        if len(cmd.get_names("selections", enabled_only=1)) == 0:
            print("[FlexWizard] Please enable a selection.")
            return

        selection_name = cmd.get_names("selections", enabled_only=1)[0]

        # The selection must live entirely in the mapped object.
        for atom in cmd.index(selection_name):
            if atom[0] != self.object_name:
                print(
                    "[FlexWizard] Selection contains atoms from objects other "
                    f"than {self.object_name.upper()}. Please fix."
                )
                return

        atom_ids = [a.id - 1 for a in cmd.get_model(selection_name).atom]
        model = cmd.index(selection_name)[0][0]
        selected = set(atom_ids)

        for atom in atom_ids:
            for neighbor in self.bond_map.get(atom, set()):
                if neighbor not in selected:
                    continue
                sel1 = f"model {model} and id {atom + 1}"
                sel2 = f"model {model} and id {neighbor + 1}"
                self._apply_joint_type(atom, neighbor, self.joint_type, sel1, sel2)

    def show_help(self):
        """Print a short usage reminder."""
        print(80 * "=")
        print("\nFlexWizard Help Message:")
        print(
            "The easiest way to set joint types is to make a selection, then "
            "set the"
        )
        print('desired joint type, then click the "Set Selection" button. After')
        print('the desired joint types have been set, use the "Save Flexibility')
        print('File" button to save the .flex for the current world. If you wish')
        print('to manually modify joint types, enable the "Manual Selection"')
        print("mode and click on atoms forming the bonds you wish to modify.")
        print("By default, all bonds are set to rigid, unless modified.\n")
        print(80 * "=")


# -------------------- public entry point --------------------
flex_manager = FlexManager()


def start_flex_wizard(prmtop=None, inpcrd=None):
    """Activate the bond-picking wizard.

    Pass ``prmtop`` (and optionally ``inpcrd``) to detect AMBER dihedrals on the
    fly at startup -- no sidecar file. If ``inpcrd`` is omitted, a sibling
    coordinate file next to the prmtop is used. Without ``prmtop`` the wizard
    starts in manual mode; you can still click "Detect Dihedrals (AMBER)" later.
    """
    global flex_manager
    cmd.set_wizard(FlexWizard(flex_manager, prmtop, inpcrd))


def stop_flex_wizard():
    """Deactivate the wizard and clear any bond coloring."""
    cmd.unset_bond("stick_color", "all")
    cmd.unset_bond("line_color", "all")
    cmd.set_wizard()
    # Force PyMOL to drop any lingering JointController edit state.
    cmd.edit_mode("on")
    cmd.edit_mode("off")
    print("[FlexWizard] Flexibility picking wizard deactivated.")
