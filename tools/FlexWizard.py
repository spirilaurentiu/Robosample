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
from pymol.Qt import QtWidgets

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
    "BallF"      : "Ball Joint (Fixed Body)",
    "BallM"      : "Ball Joint (Mobile Body)",
    "Cartesian" : "Cartesian Joint"
}

MANUAL_MODE_NAMES = {
    0: "Off",
    1: "On"
}

POLL_INTERVAL_MS = 300

### Helper Functions

def are_atoms_bonded(model, idx1, idx2):
    """Return True if the two atom indices are directly bonded in this model.
    """
    ## Get all atoms that are within 1 bond from idx1
    cmd.select("_sele", "model {} and neighbor index {}".format(model, idx1))
    for idx in cmd.index("_sele"):
        if idx[1] == idx2:
            cmd.delete("_sele")
            return True
    cmd.delete("_sele")
    return False


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

    def set_joint(self, a1, a2, jointData):
        """
        Given a jointType, set the bond between two atoms (a1 and a2) to that jointType.
        """
        bond = self._normalize_bond(a1, a2)
        if (len(jointData) == 1):
            ## Only update jointType
            self.joints[bond][0] = jointData[0] 
        else:
            ## Set the entire jointData
            self.joints[bond] = jointData 

class FlexWizard(Wizard):
    """
    Wizard for interactive bond picking.
    Left-click on two atoms to toggle bond type.
    Clears selections after each pair to allow consecutive picks.
    """
    def __init__(self, flex_manager):
        print ("Welcome to the Robosample Flexibility Wizard!")
        print ("Undo has been disabled, as it interferes with wizard functionality.")
        cmd.undo_disable()
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
        cmd.edit_mode("off")

        ## Set status (internal flag)
        self.status=0

        ## Define lists for dropdown menus, as per the measurement.py example
        smm = []
        smm.append([2, "Joint Type", ""])
        for a in JOINT_TYPES:
            smm.append([1, JOINT_NAMES[a], 'cmd.get_wizard().set_jointType("'+a+'")'])
        self.menu['jointType'] = smm

    def build_bond_map(self):

        objectList = cmd.get_names("objects", enabled_only=1)
        if (len(objectList) != 1):
            print ("Please enable an object for which to build the bond map.")
            self.bond_map = None
            return
        self.objectName = objectList[0]
        print ("Building bond map for {}.".format(self.objectName.upper()))
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

            ## By default, we set all bonds to "Rigid"
            flex_manager.set_joint(a1, a2, [self.jointType, 
                                            atom1_name, atom1_resname, atom1_resid,
                                            atom2_name, atom2_resname, atom2_resid])

        self.bond_map = bond_map

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
            print(f"[FlexWizard] Atoms {idx1} and {idx2} are not bonded.")
            self.restore()
            return

        ## If all tests are passed, select joint type
        flex_manager.set_joint(idx1, idx2, [self.jointType])

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

        print ("Flex file saved to: {}".format(filename))

    def restore(self):
        """
        This function restores the Wizard to its initial state.
        """
        cmd.delete("_pk1")
        cmd.delete("_pk2")
        self.status = 0
        cmd.refresh_wizard()
        cmd.unpick()

    def set_jointType(self, jointType):
        """
        Set the joint type that we want to set upon a bond.
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

            # Button to display help message
            [2, 'Help Me!', 'cmd.get_wizard().helpMe()', ''],

            # Button for "Done"
            [2, 'Done', 'stop_flex_wizard()', '']
            
        ]

    def set_selection(self):
        """
        Given a selection, set all bonds between all atoms in the selection to the current joint type
        """
        
        ## First we check if the bond map has been built
        if (self.bond_map == None):
            print ("No bond map built. Please enable an object and build the bond map first.")
            return 

        ## Check if any selections are enabled
        if (len(cmd.get_names("selections", enabled_only=1)) == 0):
            print ("Please enable a selection.")
            return


        ## Get all atoms in the selection
        selectionName = cmd.get_names("selections", enabled_only=1)[0]

        ## First, we check if the selection is only contained in the object
        ## for which we built the bond map
        for atom in cmd.index(selectionName):
            if (atom[0] != self.objectName):
                print ("Selection contains objects from other objects than {}. Please fix.".format(self.objectName.upper()))
                return

        listOfIxs = [a.id-1 for a in cmd.get_model(selectionName).atom]
        model = cmd.index(selectionName)[0][0]

        for ix in range(len(listOfIxs)):
            for jx in self.bond_map.get(listOfIxs[ix], []):
                if (jx in listOfIxs):
                    flex_manager.set_joint(listOfIxs[ix], jx, [self.jointType])
                    ## Color the joint accordignly
                    cmd.select("_tempsele1", "model {} and id {}".format(model, listOfIxs[ix]+1))
                    cmd.select("_tempsele2", "model {} and id {}".format(model, jx+1))

                    cmd.set_bond("stick_color", JOINT_COLORS[self.jointType], "_tempsele1", "_tempsele2")
                    cmd.set_bond("line_color", JOINT_COLORS[self.jointType], "_tempsele1", "_tempsele2")

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
    print("Flexibility picking wizard deactivated.")
