"""
PyMOL Flex Editor plugin (FlexWizard)

How to use
----------
1. In PyMOL: run this file (File -> Run or `run /path/to/pymol_flex_plugin.py`).
2. From the PyMOL Python prompt run:
       start_flex_wizard()


What it does (short)
--------------------
- Lets you define "worlds" (sets of flex assignments).
- While in Picking (3-Button Editing) mode you pick two atoms (pk1/pk2)
  to select a bond. The plugin cycles that bond through joint types for the
  currently selected world.
- Bond identity in .flex files is written using atom source IDs (cmd.id_atom)
  (so the bond is saved by atom indices rather than PyMOL-internal indices).

"""

from pymol import cmd
from pymol.wizard import Wizard
from pymol.Qt import QtWidgets, QtCore

import math

# User-configurable joint types + colors
JOINT_TYPES = ["Rigid", "Pin", "BallF", "BallM", "Cartesian"]
JOINT_COLORS = {
    "Rigid": "0x696969",
    "Pin": "0x006400",
    "BallM": "0x00BFFF",
    "BallF": "0x00FFBF",
    "Cartesian": "0xFFFF66"
}

JOINT_NAMES = {
    "Rigid"     : "Rigid Joint",
    "Pin"       : "Pin Joint",
    "BallF"     : "Ball Joint (Fixed Body)",
    "BallM"     : "Ball Joint (Mobile Body)",
    "Cartesian" : "Cartesian Joint"
}

MANUAL_MODE_NAMES = {
    0: "Off",
    1: "On"
}

POLL_INTERVAL_MS = 300

class JointControlWindow(QtWidgets.QWidget):
    def __init__(self, joints, parent=None):
        """
        joints: list of tuples (a1, a2, a3, a4)
        Each atom index is a PyMOL atom ID (not model index)
        """
        super().__init__(parent)
        self.setWindowTitle("Joint Control Panel")
        self.setLayout(QtWidgets.QVBoxLayout())
        self.joints = joints
        n = len(joints)

        # Store parameters
        self.initial_angles = []
        self.sliders = []
        self.angles = [0.0 for _ in range(n)]
        self.centers = [0.0 for _ in range(n)]     # center value for rocking
        self.speeds = [0.0 for _ in range(n)]
        self.amps = [10.0 for _ in range(n)]       # default amplitude
        self.phases = [0.0 for _ in range(n)]      # internal phase offset
        self.speed_boxes = []  # list of QDoubleSpinBox per joint
        self.amp_boxes = []

        # Outer layout
        self.setLayout(QtWidgets.QVBoxLayout())

        # Scroll area
        scroll_area = QtWidgets.QScrollArea()
        scroll_area.setWidgetResizable(True)
        self.layout().addWidget(scroll_area)

        # Inner widget that contains all joint rows
        self.scroll_widget = QtWidgets.QWidget()
        self.scroll_layout = QtWidgets.QVBoxLayout()
        self.scroll_widget.setLayout(self.scroll_layout)
        scroll_area.setWidget(self.scroll_widget)

        # Build GUI
        for i, (a1, a2, a3, a4) in enumerate(joints):
            row = QtWidgets.QHBoxLayout()

            label = QtWidgets.QLabel(f"Joint {i+1}: {a1}-{a2}-{a3}-{a4}")
            label.setFixedWidth(180)

            slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
            slider.setRange(-180, 180)
            initial_angle = cmd.get_dihedral(f"id {a1}", f"id {a2}", f"id {a3}", f"id {a4}")
            self.initial_angles.append(initial_angle)
            self.centers[i] = initial_angle
            self.angles[i] = initial_angle
            slider.setValue(int(initial_angle))
            slider.setSingleStep(5)
            slider.setMinimumWidth(200)
            slider.valueChanged.connect(lambda val, idx=i: self._set_center(idx, val))

            speed_box = QtWidgets.QDoubleSpinBox()
            speed_box.setRange(0.0, 30.0)
            speed_box.setSingleStep(0.5)
            speed_box.setValue(0.0)
            speed_box.setSuffix(" °/step")
            speed_box.setFixedWidth(90)
            speed_box.valueChanged.connect(lambda val, idx=i: self._set_speed(idx, val))

            amp_box = QtWidgets.QDoubleSpinBox()
            amp_box.setRange(0.0, 180.0)
            amp_box.setSingleStep(5.0)
            amp_box.setValue(10.0)
            amp_box.setSuffix(" °")
            amp_box.setFixedWidth(80)
            amp_box.valueChanged.connect(lambda val, idx=i: self._set_amplitude(idx, val))

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

        # Timer is created but not started yet
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self._update_angles)
        # Do not start it yet — wait for "Start All"

        # Pause/Resume button
        self.pause_button = QtWidgets.QPushButton("Start All")
        self.pause_button.clicked.connect(self.toggle_pause)
        self.layout().addWidget(self.pause_button)

        # Allow User to save the new angles.
        commit_button = QtWidgets.QPushButton("Commit Changes")
        commit_button.clicked.connect(self.commit_changes)
        self.layout().addWidget(commit_button)

        # Close button
        close_button = QtWidgets.QPushButton("Close")
        close_button.clicked.connect(self.close)
        self.layout().addWidget(close_button)

        # Batch controls
        batch_layout = QtWidgets.QHBoxLayout()

        # Set all speeds
        self.all_speed_button = QtWidgets.QPushButton("Set all speeds to")
        self.all_speed_button.clicked.connect(self.set_all_speeds)
        batch_layout.addWidget(self.all_speed_button)
        self.all_speed_spin = QtWidgets.QDoubleSpinBox()
        self.all_speed_spin.setRange(0.0, 30.0)
        self.all_speed_spin.setSingleStep(0.5)
        self.all_speed_spin.setValue(0.0)
        self.all_speed_spin.setSuffix(" °/step")
        batch_layout.addWidget(self.all_speed_spin)

        # Set all amplitudes
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
        """Commit the current angles as the new reference."""
        self.initial_angles = self.angles
        print("[FlexWizard] Current angles committed as new reference.")


    def closeEvent(self, event):
        """Called when the window is closed."""
        if hasattr(self, "timer") and self.timer.isActive():
            self.timer.stop()

        ## Restore initial angles
        for idx, (a1, a2, a3, a4) in enumerate(self.joints):
            cmd.set_dihedral(f"id {a1}", f"id {a2}", f"id {a3}", f"id {a4}", self.initial_angles[idx])
        print("[FlexWizard] Initial angles restored.")

        # Clear reference in parent wizard, if exists
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
        val = self.all_speed_spin.value()
        for i in range(len(self.speeds)):
            self._set_speed(i, val)
            self.speed_boxes[i].blockSignals(True)
            self.speed_boxes[i].setValue(val)
            self.speed_boxes[i].blockSignals(False)

    def set_all_amps(self):
        val = self.all_amp_spin.value()
        for i in range(len(self.amps)):
            self._set_amplitude(i, val)
            self.amp_boxes[i].blockSignals(True)
        self.amp_boxes[i].setValue(val)
        self.amp_boxes[i].blockSignals(False)

    # --- GUI event handlers ---

    def _set_center(self, idx, value):
        """User moved slider manually → update the rocking center."""
        self.centers[idx] = value
        self.angles[idx] = value
        self._apply_dihedral(idx)

    def _set_speed(self, idx, speed):
        """Speed > 0 enables rocking."""
        self.speeds[idx] = speed

    def _set_amplitude(self, idx, amp):
        """Set amplitude for local oscillation."""
        self.amps[idx] = amp

    # --- Core update logic ---

    def _update_angles(self):
        """Oscillate each torsion around its center."""
        for idx, (a1, a2, a3, a4) in enumerate(self.joints):
            speed = self.speeds[idx]
            if speed <= 0.0:
                continue

            # advance phase and compute new angle
            self.phases[idx] += speed * 0.05
            amp = self.amps[idx]
            center = self.centers[idx]
            new_angle = center + amp * math.sin(self.phases[idx])

            # apply + sync slider
            self.angles[idx] = new_angle
            self.sliders[idx].blockSignals(True)
            self.sliders[idx].setValue(int(new_angle))
            self.sliders[idx].blockSignals(False)
            self._apply_dihedral(idx)

    def _apply_dihedral(self, idx):
        """Apply torsion to PyMOL."""
        a1, a2, a3, a4 = self.joints[idx]
        angle = self.angles[idx]
        cmd.set_dihedral(f"id {a1}", f"id {a2}", f"id {a3}", f"id {a4}", float(angle))

class FlexManager:
    """
    Keeps track of flexible bonds for multiple 'worlds',
    and provides methods to toggle or save them.
    """

    def __init__(self):
        self.joints = {}

    def _normalize_bond(self, a1, a2):
        """Always order atom indices so bonds are unique regardless of direction."""
        return tuple(sorted((int(a1), int(a2))))

    def set_joint(self, a1, a2, jointData, bond_map):
        """
        Given a jointType, set the bond between two atoms (a1 and a2) to that jointType.
        """
        bond = self._normalize_bond(a1, a2)
        if (len(jointData) == 1):
            ## Only update jointType
            if (jointData[0] == "Pin"):
                if (self.get_dihedral_from_bond(bond[0], bond[1], bond_map)):
                    self.joints[bond][0] = jointData[0]
                    return self.get_dihedral_from_bond(bond[0], bond[1], bond_map)
                else:
                    self.joints[bond][0] = "Rigid"
                    return False
            else: 
                self.joints[bond][0] = jointData[0]
                return True
        else:
            ## Set the entire jointData
            self.joints[bond] = jointData 

    def get_dihedral_from_bond(self, a1, a2, bond_map):
        """
        Given two atoms, check if they can form a dihedral (i.e. if 
        neither of them is terminal). If yes, return it
        """
        a1_neighbors = bond_map.get(a1)
        a2_neighbors = bond_map.get(a2)

        a1_neighbors = [x for x in a1_neighbors if x != a2]
        a2_neighbors = [x for x in a2_neighbors if x != a1]

        if (len(a1_neighbors) == 0) or (len(a2_neighbors) == 0):
            return False
        else:
            return [a1_neighbors[0], a1, a2, a2_neighbors[0]]

class FlexWizard(Wizard):
    """
    Wizard for interactive bond picking.
    Left-click on two atoms to toggle bond type.
    Clears selections after each pair to allow consecutive picks.
    """
    def __init__(self, flex_manager):
        print ("[FlexWizard] Welcome to the Robosample Flexibility Wizard!")
        ## Older versions of PyMOL don't have the undo_disable function,
        ## so we check.
        try:
            cmd.undo_disable()
            print ("[FlexWizard] Undo has been disabled, as it interferes with wizard functionality.")
        except AttributeError:
            pass
        Wizard.__init__(self)
        ## Just in case a failed instantiation of this class was previously used
        cmd.delete("_pk1")
        cmd.delete("_pk2")

        self.jointType = self.session.get('default_mode', 'Rigid')
        self.manualMode = 0

        self.flex_manager = flex_manager
        self.atom1 = None

        ## We build the bond_map when we initialize
        self.build_bond_map()

        # Ensure PyMOL is not in edit mode
        cmd.edit_mode("on")
        cmd.edit_mode("off")

        ## Set status (internal flag)
        self.status=0

        # Set mode to viewing, not editing
        cmd.edit_mode("off")

        ## Define lists for dropdown menus, as per the measurement.py example
        smm = []
        smm.append([2, "Joint Type", ""])
        for a in JOINT_TYPES:
            smm.append([1, JOINT_NAMES[a], 'cmd.get_wizard().set_jointType("'+a+'")'])
        self.menu['jointType'] = smm

    def build_bond_map(self):

        objectList = cmd.get_names("objects", enabled_only=1)
        if (len(objectList) != 1):
            print ("[FlexWizard] Please enable an object for which to build the bond map.")
            self.bond_map = None
            return
        self.objectName = objectList[0]
        print ("[FlexWizard] Built bond map for {}.".format(self.objectName.upper()))
        model = cmd.get_model("model {}".format(self.objectName))
        bond_map = {}

        ## We need to convert indices to atom IDs
        rank_to_id = {i: atom.id for i, atom in enumerate(model.atom)}
        bond_ids = [(rank_to_id[a], rank_to_id[b]) for (a, b) in (b.index for b in model.bond)]

        for b in bond_ids:
            a1, a2 = b

            a1 -= 1
            a2 -= 1

            ## Get the metadata of each atom in bond
            atom1_name = model.atom[a1].name
            atom1_resname = model.atom[a1].resn
            atom1_resid = model.atom[a1].resi

            atom2_name = model.atom[a2].name
            atom2_resname = model.atom[a2].resn
            atom2_resid = model.atom[a2].resi
            
            bond_map.setdefault(a1, set()).add(a2)
            bond_map.setdefault(a2, set()).add(a1)

            ## By default, we set all bonds to "Rigid"
            flex_manager.set_joint(a1, a2, [self.jointType, 
                                            atom1_name, atom1_resname, atom1_resid,
                                            atom2_name, atom2_resname, atom2_resid],
                                            None)

        self.bond_map = bond_map

        ## Set a list, for the joint control panel, which we
        ## re-initialize everytime we rebuild the bond map.
        self.pin_joints = []
        cmd.unset_bond("stick_color", "all")
        cmd.unset_bond("line_color", "all")

        return 

    def do_select(self, name):
        if (self.manualMode == 0):
            return
        self.cmd.unpick()
        if (self.status == 0):
            cmd.select("_pk1", name)
            cmd.delete(name)
            self.status = 1 ## To indicate that the first atom has been picked
            self.prompt = ["Pick the second atom."]
            cmd.refresh_wizard()

        elif (self.status == 1):
            cmd.select("_pk2", name)
            cmd.delete(name)
            self.do_pick()

    def do_pick(self):

        idx1 = cmd.id_atom("_pk1")
        idx2 = cmd.id_atom("_pk2")

        idx1 -= 1
        idx2 -= 1

        model01 = cmd.index("_pk1")[0][0]
        model02 = cmd.index("_pk2")[0][0]

        ## Check if they're not the same atom selected twice
        if (idx1 == idx2):
            print("[FlexWizard] Selected the same atom twice.")
            self.restore()
            return

        ## If they're not from the same molecule, deselect
        if (model01 != model02):
            print("[FlexWizard] Atoms from different models.")
            self.restore()
            return

        ## If they're not bonded, deselect
        if not (idx2 in self.bond_map.get(idx1)):
            print("[FlexWizard] Atoms {} and {} are not bonded.".format(idx1,idx2))
            self.restore()
            return

        ## If all tests are passed, select joint type
        dihed = flex_manager.set_joint(idx1, idx2, [self.jointType], self.bond_map)
        if (type(dihed) != bool) and (dihed not in self.pin_joints):
            self.pin_joints.append([x+1 for x in dihed])
            
        # Apply bond color
        cmd.set_bond("stick_color", JOINT_COLORS[self.jointType], "_pk1", "_pk2")
        cmd.set_bond("line_color", JOINT_COLORS[self.jointType], "_pk1", "_pk2")

        self.restore()
        return

    def write_flex(self):
        """
        Write the flex file.
        """

        ## In case the bond_map has not been built 
        objectList = cmd.get_names("objects", enabled_only=1)
        self.objectName = objectList[0] 

        dialog = QtWidgets.QFileDialog()
        filename, _ = dialog.getSaveFileName(
            None,
            "Save Flexibility File",
            self.objectName,
            "Flex files (*.flex);;All files (*)"
        )

        data = flex_manager.joints
        ## We sort the joints by the first atom, so it's easier to visually scan.
        data = dict(sorted(data.items(), key=lambda item: item[0][0]))

        with open(filename, "w") as f:
            for joint in data.items():
                f.write("{:9s} {:9s} {:5s}\t#{:3s}_{}_{}-{:3s}_{}_{}\n".format(str(joint[0][0]),str(joint[0][1]), 
                                                                  joint[1][0], 
                                                                  joint[1][2],joint[1][1], joint[1][3], 
                                                                  joint[1][5],joint[1][4], joint[1][6]))

        print ("[FlexWizard] Flex file saved to: {}".format(filename))

    def restore(self):
        """
        This function restores the Wizard to its initial state.
        """
        cmd.delete("_pk1")
        cmd.delete("_pk2")
        self.status = 0
        self.toggle_manualMode()
        self.toggle_manualMode()
        cmd.refresh_wizard()
        cmd.unpick()

    def set_jointType(self, jointType):
        """
        Set the joint type variable.
        """
        if jointType in JOINT_TYPES:
            self.jointType = jointType
        self.restore()

    def toggle_manualMode(self):
        if (self.manualMode == 0):
            self.manualMode = 1
            cmd.set("mouse_selection_mode", 0)  # 0 = atomic picking
            self.prompt = ["Pick first atom to toggle a bond."]
        else:
            self.manualMode = 0
            cmd.set("mouse_selection_mode", 1)  # 1 = residue picking
            self.prompt = None
        cmd.refresh_wizard()

    ### PANEL FUNCTIONS ###

    def get_panel(self):
        """
        Defines the PyMOL Wizard control panel layout.
        Each entry is a list of (type, label, command).
        """
        return [
            # Text-only label (non-clickable)
            [1, 'FlexWizard - by T. A. Sulea.', '', ''],
            
            # Menu for selecting joint type
            [3, 'Joint type: '+JOINT_NAMES[self.jointType], 'jointType'],

            # Button for manual selection toggling
            [2, 'Manual Selection ' + MANUAL_MODE_NAMES[self.manualMode], 'cmd.get_wizard().toggle_manualMode()', ''],

            # Button to rebuild bond graph
            [2, 'Rebuild Bond Map', 'cmd.get_wizard().build_bond_map()', ''],

            # Button for "Save Flex"
            [2, 'Save Flexibility File', 'cmd.get_wizard().write_flex()', ''],

            # Button for "Set Selection"
            [2, 'Set Selection', 'cmd.get_wizard().set_selection()', ''],

            # Button to open the joint control panel
            [2, 'Open Joint Control', 'cmd.get_wizard().open_joint_control_panel()', ''],

            # Button to display help message
            [2, 'Help Me!', 'cmd.get_wizard().helpMe()', ''],

            # Button for "Done"
            [2, 'Done', 'stop_flex_wizard()', '']
            
        ]

    def open_joint_control_panel(self):
        # Example: retrieve your stored list of Pin-type torsions
        joints = self.pin_joints  # list of (a1, a2, a3, a4)
        print (joints)
        if (len(joints) == 0):
           print ("Please select at least one Pin-type joint.")
           return 

        self.joint_window = JointControlWindow(joints)
        self.joint_window.show()

    def set_selection(self):
        """
        Given a selection, set all bonds between all atoms in the selection to the current joint type
        """
        
        ## First we check if the bond map has been built
        if (self.bond_map == None):
            print ("[FlexWizard] No bond map built. Please enable an object and build the bond map first.")
            return 

        ## Check if any selections are enabled
        if (len(cmd.get_names("selections", enabled_only=1)) == 0):
            print ("[FlexWizard] Please enable a selection.")
            return


        ## Get all atoms in the selection
        selectionName = cmd.get_names("selections", enabled_only=1)[0]

        ## First, we check if the selection is only contained in the object
        ## for which we built the bond map
        for atom in cmd.index(selectionName):
            if (atom[0] != self.objectName):
                print ("[FlexWizard] Selection contains objects from other objects than {}. Please fix.".format(self.objectName.upper()))
                return

        listOfIxs = [a.id-1 for a in cmd.get_model(selectionName).atom]
        model = cmd.index(selectionName)[0][0]

        for ix in range(len(listOfIxs)):
            for jx in self.bond_map.get(listOfIxs[ix], []):
                if (jx in listOfIxs):
                    dihed = flex_manager.set_joint(listOfIxs[ix], jx, [self.jointType], self.bond_map)
                    if dihed:
                        color = JOINT_COLORS[self.jointType]
                        ## Add to self.pin_joints, if self.jointType is set to "Pin"
                        if (type(dihed) != bool) and (dihed not in self.pin_joints):
                            self.pin_joints.append([x+1 for x in dihed])
                    else: 
                        color = JOINT_COLORS["Rigid"]
                    ## Color the joint accordignly
                    cmd.select("_tempsele1", "model {} and id {}".format(model, listOfIxs[ix]+1))
                    cmd.select("_tempsele2", "model {} and id {}".format(model, jx+1))

                    cmd.set_bond("stick_color", color, "_tempsele1", "_tempsele2")
                    cmd.set_bond("line_color", color, "_tempsele1", "_tempsele2")

                    cmd.delete("_tempsele1")
                    cmd.delete("_tempsele2")

                    


    def helpMe(self):
        print (80*"=")
        print ("\nFlexWizard Help Message:")
        print ("The easiest way to set joint types is to make a selection, then set the ")
        print ("desired joint type, then click the \"Set Selection\" button. After the ")
        print ("desired joint types have been set, use the \"Save Flexibility File\" ")
        print ("button to save the .flex for the current world. If you wish to ")
        print ("manually modify joint types, enable the \"Manual Selection\" , ")
        print ("mode and click on atoms forming the bonds you wish to modify.")
        print ("By default, all bonds are set to rigid, unless modified.\n")
        print (80*"=")


# -------------------- public entry point --------------------
flex_manager = FlexManager()

def start_flex_wizard():
    """Activate the bond-picking wizard."""
    global flex_manager
    cmd.set_wizard(FlexWizard(flex_manager))

def stop_flex_wizard():
    cmd.unset_bond("stick_color", "all")
    cmd.unset_bond("line_color", "all")
    cmd.set_wizard()
    ## In order to remove any elements from the JointController.
    cmd.edit_mode("on")
    cmd.edit_mode("off")
    print("[FlexWizard] Flexibility picking wizard deactivated.")

def set_flexibility(atom_id1, atom_id2, joint_type):
    wizard = cmd.get_wizard()
    cmd.select("_pk1", f"id {atom_id1}")
    cmd.select("_pk2", f"id {atom_id2}")
    old_type = wizard.jointType
    wizard.jointType = joint_type
    wizard.do_pick()
    wizard.jointType = old_type   # restore previous default
