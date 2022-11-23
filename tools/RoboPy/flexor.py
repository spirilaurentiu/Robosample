import numpy as np
import mdtraj as md
import json, sys


class Flexor:
    def __init__(self, mdtraj_obj):
        '''
        Parameters
        ----------
        mdtraj_obj: mdtraj.core.trajectory.Trajectory

        Notes
        ----------
        This class ties together everything needed to generate appropriate
        flexibility files, compatible with the Robosample software.
        It uses the "graph" representation of the system to get molecular
        topology, from where connectivity data is obtained.
        It uses the compute_contacts function of MDTraj to get inter-residue
        distances, for computing structural neighbors, useful when generating
        flexibility files.

        '''
        self._MDTrajObject = mdtraj_obj
        self._edges = list(self._MDTrajObject.topology.to_bondgraph().edges)
        self._DB = {}

    def loadFlexDB(self, FlexFN):
        '''
        Parameters
        ----------
        FlexFN: string
            Path-like string pointing to file containing flexibility info of
            residues.

        Returns
        ----------
        None

        Notes
        ----------
        The flexibility database files that are loaded are structured in a simple manner.
        All names are based on the AMBER forcefields, since that is what Robosample currently
        supports. Users can create their own databases for Glycans, DNA/RNA or custom residues
        if they wish.
        '''

        with open(FlexFN) as f:
            data = f.read()
            js = json.loads(data)
            print("Flex DB containing \"{}\" loaded! ".format(js["TITLE"]))
            js.pop("TITLE")
            self._DB = {**self._DB, **js}

    def addWorld(self, range, subset, type, jointType, FNOut, distanceCutoff=0, sasa_percentile=0):
        '''
        Parameters
        ----------
        range:  string
            DSL selection string

        distanceCutoff: float
            Include other atoms in the system that are within this distance
            from any atom in the selection defined by range

        subset: list of strings
            Defines what parts of the system to include in the flexibility file
            Options:
                side: sidechains
                phi: phi angle (N-CA)
                psi: psi angle (CA-C)
                omega: peptide bond (N-C)
                rama: macro for phi+psi
                all: macro for all of the above.

        type:  string
            Generate one flexibility file with all lines ("stretch"), or each line in
            a separate file ("roll")

        jointType: string
            Type of joint to use. Currently implemented:
                pin, ball and cartesian

        FNOut: string
            Path-like string to a file where the flexibilities will be written.

        sasa_percentile: float
            If supplied, will filter atoms from range so that only atoms that are in the
            sasa_percentile'th percentile of accesibility will be kept
            Example: if sasa_percentile is 0.9, then only the atoms that are in the top 10%
            accesibility-wise AND in the selected range will be kept.

        Returns
        ----------
        None
        '''

        jointCount = 0
        fileOut = open(FNOut, "w")
        jointType = jointType.lower()
        selIx = list(self._MDTrajObject.topology.select(range))
        ## If sasa_percentile has been supplied, we do the filtering now.
        if (sasa_percentile > 0):
            ## Compute
            sasa = md.shrake_rupley(self._MDTrajObject)[0]

            ## Normalize
            sasa = (sasa - np.min(sasa))/(np.max(sasa) - np.min(sasa))

            ## Filter
            for i,v in enumerate(selIx):
                if sasa[v] < sasa_percentile:
                    selIx.pop(i)
        validJointAtoms = []
        if ("all" in subset):
            subset = ["side", "rama"]
        if ("rama" in subset):
            subset.remove("rama")
            subset.append("phi")
            subset.append("psi")
            subset.append("omega")
        for edge in self._edges:
            if ((edge[0].index in selIx) and (edge[1].index in selIx)):
                if (jointType == "cartesian"):
                    fileOut.write("{:5d} {:5d} Cartesian {} {}\n".format(edge[0].index, edge[1].index, edge[0], edge[1]))
                    validJointAtoms.append(edge[0].index)
                    validJointAtoms.append(edge[1].index)
                    jointCount += 1
                else:
                    for subsetType in subset:
                        if (self.check_joint(edge[0], jointType, subsetType) and
                            self.check_joint(edge[1], jointType, subsetType)):
                            fileOut.write("{:5d} {:5d} {} {} {}\n".format(edge[0].index, edge[1].index, jointType.capitalize(),
                                                          edge[0], edge[1]))
                            validJointAtoms.append(edge[0].index)
                            validJointAtoms.append(edge[1].index)
                            jointCount += 1

        if (distanceCutoff>0):
            ## First we get all atoms that are within 'distanceCutoff' of atoms in joints
            nearbyAtoms = md.compute_neighbors(self._MDTrajObject, distanceCutoff, validJointAtoms)
            ## Then we remove the atoms that are already in the validJointAtoms selection
            validJointAtoms = list(set(validJointAtoms))
            for _ in validJointAtoms:
                ix = np.argwhere(nearbyAtoms==_)
                nearbyAtoms = np.delete(nearbyAtoms,ix)
            ## Then we do the check we did above, but for the new selection
            fileOut.write("### Neighboring Joints ###\n")
            for edge in self._edges:
                if ((edge[0].index in nearbyAtoms) and (edge[1].index in nearbyAtoms)):
                    if (jointType == "cartesian"):
                        fileOut.write(
                            "{:5d} {:5d} Cartesian {} {}\n".format(edge[0].index, edge[1].index, edge[0], edge[1]))
                        jointCount += 1

                    else:
                        for subsetType in subset:
                            if (self.check_joint(edge[0], jointType, subsetType) and
                                    self.check_joint(edge[1], jointType, subsetType)):
                                fileOut.write("{:5d} {:5d} {} {} {}\n".format(edge[0].index, edge[1].index,
                                                                              jointType.capitalize(),
                                                                              edge[0], edge[1]))
                                jointCount += 1

            print ("Flex file {} generated! ({} joints)".format(FNOut, jointCount))
            fileOut.close()
            
    ## Helper functions
    def check_joint(self, atom, jointType, subset):
        '''
        Parameters
        ----------
        atom: mdtraj.core.topology.Atom type
            Atom whose connectivity to check against the databases loaded

        jointType: string
            Same as jointType argument from addWorld function.

        subset: string
            Same as subset argument from addWorld function.

        Returns
        ----------
        bool
            True if supports the proposed joint type.
        '''

        if (subset == "side"):
            if (atom.name in self._DB[str(atom.residue)[:3]]["SC_{}".format(jointType)]):
                return True
        if (subset == "phi"):
            if (atom.name in self._DB[str(atom.residue)[:3]]["phi"]):
                return True
        if (subset == "psi"):
            if (atom.name in self._DB[str(atom.residue)[:3]]["psi"]):
                return True
        if (subset == "omega"):
            if (atom.name in self._DB[str(atom.residue)[:3]]["omega"]):
                return True
