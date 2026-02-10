import numpy as np
import mdtraj as md
import networkx as nx
import os


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
        topology, from which connectivity data is obtained.
        It uses the compute_contacts function of MDTraj to get inter-residue
        distances, for computing structural neighbors, useful when generating
        flexibility files.
        '''
        self._MDTrajObject = mdtraj_obj
        self._edges = list(self._MDTrajObject.topology.to_bondgraph().edges)

    def addWorld(self, range, subset,jointType, FNOut, distanceCutoff=0, sasa_value=-1, rolling=False):
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
                phi: phi angle (N-CA)
                psi: psi angle (CA-C)
                omega: peptide bond (N-C)
                rama: macro for phi+psi+omega
                side: sidechain atoms 
                all: macro for all of the above.
                coils: only selects coil regions of the protein, intersected with
                       whatever range the user made
                
                ligand: DEPRECATED - Automatically uses n_bonds to det. joint types.
                        and check_cycle to determine wether to put a flexibility somewhere.

                        special subset, used for ligands, or other small molecules.
                        made to be used with molecules with non-standard nomenclature,
                        that use n_bonds to determine joint types. Also includes Glycans
                        or DNA molecules.
                
        jointType: string
            Type of joint to use. Currently implemented:
                pin, ball and cartesian
        FNOut: string
            Path-like string to a file where the flexibilities will be written.
        sasa_value: float
            If supplied, will filter atoms from range so that only atoms that have a higher
            SASA than this are kept.
        rolling:  bool
            if True, this world's flexibility will be split across multiple files (one joint
            per file)
        Returns
        ----------
        fileNameIX:list 
            list containing all the flexiblity files generated 
        '''

        jointType = jointType.lower()
        fileNameIX = []
        ## When assigning Pin or Ball flexibility to a joint, we have to account for 
        ## cyclical structures, where 1 DoF is lost. As such, we randomly assign a bond
        ## as "weld" in every cycle we have, and keep track of that in this list.

        ## Proline always has the N-CD bond rigid (due to how the script is written),
        ## so we can skip PRO when doing the Cycle checker (assign -1 to the PRO cycles)

        cycles = nx.cycle_basis(self._MDTrajObject.topology.to_bondgraph())
        cycleVisits = [len(x)-1 if x[0].residue.name == "PRO" else len(x) for x in cycles]
        if ("side" in subset):
            cycleVisits = [-1 if x[0].residue.name == "PRO" else len(x) for x in cycles]
        ## Just so the names of the flex files are shorter 
        if (jointType == "cartesian"):
            jointType = "cart"

        jointCount = 0
        fileCounter = 0
        if (rolling == False):
            fileOut = open("{}.flex".format(FNOut), "w")
            fileNameIX.append(os.path.basename(fileOut.name))
            

        selIx = list(self._MDTrajObject.topology.select(range))
        ## If sasa_value has been supplied, we do the filtering now.
        if (sasa_value >= 0):
            print ("SASA filtering enabled.")
            ## Compute
            sasa = md.shrake_rupley(self._MDTrajObject)[0]

            ## Filter
            buried_ixs = []
            for i,v in enumerate(selIx):
                if sasa[v] <= sasa_value:
                    buried_ixs.append(i)
            for i in buried_ixs:
                selIx.remove(i)
        
        validJointAtoms = []
        if ("coils" in subset):
            if (len(subset) == 1):
                subset = ["coils", "all"]
            else:
                subset.append("coils")
            dssp = md.compute_dssp(self._MDTrajObject, simplified=True)[0]

        if ("all" in subset):
            subset.remove("all")
            subset.append("rama")
            subset.append("side")
        if ("rama" in subset):
            subset.remove("rama")
            subset.append("phi")
            subset.append("psi")
            subset.append("omega")


        for edge in self._edges:
            if ((edge[0].index in selIx) and (edge[1].index in selIx)):
                if ("coils" in subset):
                    if (dssp[edge[0].residue.index] != "C" and dssp[edge[1].residue.index] != "C"):
                        continue
                if (jointType == "cart"):
                    if (rolling == True):
                        fileOut = open("{}.{}.flex".format(FNOut, fileCounter), "w")
                        fileCounter += 1
                    fileOut.write("{:5d} {:5d} Cartesian #{} {}\n".format(edge[0].index, edge[1].index, edge[0], edge[1]))
                    validJointAtoms.append(edge[0].index)
                    validJointAtoms.append(edge[1].index)
                    jointCount += 1
                else:        
                    for subsetType in subset:

                        if ((self.check_joint(edge[0], jointType, subsetType) and
                            self.check_joint(edge[1], jointType, subsetType)) and not
                            (self.check_cycle(edge[0],edge[1], cycleVisits)) or 
                            (self.isProline(edge[0], edge[1], "N", "CD") and subsetType == "side")):
                            if (rolling == True):
                                fileOut = open("{}.{}.flex".format(FNOut, fileCounter), "w")
                                fileNameIX.append(os.path.basename(fileOut.name))
                                fileCounter += 1
                            fileOut.write("{:5d} {:5d} {} #{} {}\n".format(edge[0].index, edge[1].index, jointType.capitalize(),
                                                            edge[0], edge[1]))
                            validJointAtoms.append(edge[0].index)
                            validJointAtoms.append(edge[1].index)
                            jointCount += 1

        if (distanceCutoff>0):
            # First we get all atoms that are within 'distanceCutoff' of atoms in joints
            nearbyAtoms = md.compute_neighbors(self._MDTrajObject, distanceCutoff, validJointAtoms)
            # Then we remove the atoms that are already in the validJointAtoms selection
            validJointAtoms = list(set(validJointAtoms))
            for _ in validJointAtoms:
                ix = np.argwhere(nearbyAtoms==_)
                nearbyAtoms = np.delete(nearbyAtoms,ix)
            # Then we do the check we did above, but for the new selection
            fileOut.write("### Neighboring Joints ###\n")
            for edge in self._edges:
                if ((edge[0].index in nearbyAtoms) and (edge[1].index in nearbyAtoms)):
                    if (jointType == "cart"):
                        fileOut.write(
                            "{:5d} {:5d} Cartesian #{} {}\n".format(edge[0].index, edge[1].index, edge[0], edge[1]))
                        jointCount += 1

                    else:
                        for subsetType in subset:
                            if ((self.check_joint(edge[0], jointType, subsetType) and
                            self.check_joint(edge[1], jointType, subsetType)) and not
                            (self.check_cycle(edge[0], edge[1], cycleVisits))):
                                fileOut.write("{:5d} {:5d} {} #{} {}\n".format(edge[0].index, edge[1].index,
                                                                                jointType.capitalize(),
                                                                                edge[0], edge[1]))
                                jointCount += 1

        ## There is a possibility that we do Rolling selection and 
        ## there are no valid atoms to select, in which case no fileOut
        ## will actually be opened. In that case, we don't need to close 
        ## it and will return an empty fileNameIX
        if (rolling == True and jointCount == 0):
            return (fileNameIX)    

        fileOut.close()
        print ("Flex file {} generated! ({} joints)".format(fileOut.name, jointCount))
        return (fileNameIX)


    # Helper functions
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

        phi   = ["N", "CA"]
        psi   = ["CA", "C"]
        omega = ["N", "C"]

 
        if (atom.n_bonds > 1):
            if ((subset == "phi") and (atom.name in phi)):
                return True
            if ((subset == "psi") and (atom.name in psi)):
                return True
            if ((subset == "omega") and (atom.name in omega)):
                return True
            if ((subset == "side") and (atom.name not in ["N", "C"])):
                return True
        return False

    def check_cycle(self, atom1, atom2, cycleList):
        '''
        Parameters
        ----------
        atom1: mdtraj.core.topology.Atom type
            first atom of the pair to check
        atom2: mdtraj.core.topology.Atom type
            second atom of the pair to check
        Returns
        ----------
        bool
            True if atoms are in the same cycle, False if not.
            
        Notes
        ---------

        Added in v2: the addition of a cycle means the loss of one
        Degree of Freedom. As such, we randomly "Weld" one of the bonds,
        but leave the rest as Pin. We choose the welded bond randomly, evaluating
        each bond in a cycle and welding it with a probability of 1/(remaining bonds),
        if a bond to weld was not yet selected.
       
        In the special case of Proline, we want to leave the backbone bond flexible, so
        that one has the value of "remaining bonds" set to 4 by default and the N-CA bond
        is always left flexible.

        This function will get tripped up by cysteine bonds, since
        those close a large cycle, via the backbone. As such, we check
        if the cycles in the nx.cycle_basis contain ANY SG atoms before
        doing any other check. If it does, we skip it.

        '''
        cycles = nx.cycle_basis(self._MDTrajObject.topology.to_bondgraph())
        for cycle in cycles:
            if ([x.name for x in cycle if x.name=='SG']):
                continue
            if ((atom1 in cycle) and (atom2 in cycle)):
                ## Weld randomly if needed:
                cycleIx = cycles.index(cycle)
                if (cycleList[cycleIx] != -1): ## Meaning none of the bonds in this cycle have been welded
                    ## Check if it's the N-CA PRO bond, if so then skip (leave it flexible)
                    if (self.isProline(atom1,atom2,"N","CA")):
                        return False
                    weldProbability = 1/cycleList[cycleIx]
                    randomNormalValue = np.random.default_rng().uniform()
                    if (randomNormalValue < weldProbability):
                        ## Weld this one
                        cycleList[cycleIx] = -1
                        return True
                    else:
                        cycleList[cycleIx] += -1
                return False
        return False

    def isProline(self, atom1, atom2, atom1Name, atom2Name):
        if ((atom1.residue.name == "PRO" and atom2.residue.name == "PRO") and
            (atom1Name in [atom1.name, atom2.name] and atom2Name in [atom1.name, atom2.name])):
            return True
        else:
            return False
