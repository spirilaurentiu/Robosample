import numpy as np
import mdtraj as md
import networkx as nx
import os
import robosample as robosample


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
        We also remove the CYS-CYS bonds, since those trip up the 
        check_cycle function.
        '''
        self._MDTrajObject = mdtraj_obj
        self._edges = list(self._MDTrajObject.topology.to_bondgraph().edges)

    def create(self, range, subset, jointType, distanceCutoff=0, sasa_value=-1):
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
                rama: macro for phi+psi
                all: macro for all of the above.
                coils: only selects coil regions of the protein, intersected with
                       whatever range the user made
                
                ligand: DEPRECATED - Automatically uses n_bonds to det. joint types.
                        and check_cycle to determine wether to put a flexibility somewhere.

                        special subset, used for ligands, or other small molecules.
                        made to be used with molecules with non-standard nomenclature,
                        that use n_bonds to determine joint types. Also includes Glycans
                        or DNA molecules.
                
                side: To be DEPRECATED - Use the "sidechain" keyword to the range to include
                      or exclude the sidechain atoms.
                    
        jointType: string
            Type of joint to use. Currently implemented:
                pin, ball and cartesian
        sasa_value: float
            If supplied, will filter atoms from range so that only atoms that have a higher
            SASA than this are kept.
        Returns
        ----------
        None
        '''

        jointType = jointType.lower()
        flexibilities = []

        ## Just so the names of the flex files are shorter 
        if (jointType == "cartesian"):
            jointType = "cart"

        jointCount = 0

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


        for edge in self._edges:
            if ((edge[0].index in selIx) and (edge[1].index in selIx)):
                if ("coils" in subset):
                    if (dssp[edge[0].residue.index] != "C" and dssp[edge[1].residue.index] != "C"):
                        continue
                if (jointType == "cart"):
                    flexibilities.append(robosample.BondFlexibility(edge[0].index, edge[1].index, robosample.BondMobility.Translation))
                    # flexibilities.append("{{{}, {}, BondMobility::Mobility::Translation}},".format(edge[0].index, edge[1].index))
                    validJointAtoms.append(edge[0].index)
                    validJointAtoms.append(edge[1].index)
                    jointCount += 1
                else:        
                    for subsetType in subset:
                        if ((self.check_joint(edge[0], jointType, subsetType) and
                            self.check_joint(edge[1], jointType, subsetType)) and not
                            (self.check_cycle(edge[0],edge[1]))):

                            if jointType.capitalize() == 'Pin':
                                flexibilities.append(robosample.BondFlexibility(edge[0].index, edge[1].index, robosample.BondMobility.Torsion))
                                # flexibilities.append("{{{}, {}, BondMobility::Mobility::Torsion}},".format(edge[0].index, edge[1].index))
                            else:
                                raise ValueError("Joint type not implemented.")
                            # if jointType.capitalize() == 'Ball':
                            #     flexibilities.append(robosample.BondFlexibility(edge[0].index, edge[1].index, robosample.BondMobility.Ball))

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
            for edge in self._edges:
                if ((edge[0].index in nearbyAtoms) and (edge[1].index in nearbyAtoms)):
                    if (jointType == "cart"):
                        flexibilities.append(robosample.BondFlexibility(edge[0].index, edge[1].index, robosample.BondMobility.Translation))
                        # flexibilities.append("{{{}, {}, BondMobility::Mobility::Translation}},".format(edge[0].index, edge[1].index))
                        jointCount += 1
                    else:
                        for subsetType in subset:
                            if ((self.check_joint(edge[0], jointType, subsetType) and
                            self.check_joint(edge[1], jointType, subsetType)) and not
                            (self.check_cycle(edge[0], edge[1]))):
                                if jointType.capitalize() == 'PIN':
                                    flexibilities.append(robosample.BondFlexibility(edge[0].index, edge[1].index, robosample.BondMobility.Pin))
                                    # flexibilities.append("{{{}, {}, BondMobility::Mobility::Pin}},".format(edge[0].index, edge[1].index))
                                if jointType.capitalize() == 'BALL':
                                    flexibilities.append(robosample.BondFlexibility(edge[0].index, edge[1].index, robosample.BondMobility.Ball))
                                    # flexibilities.append("{{{}, {}, BondMobility::Mobility::Ball}},".format(edge[0].index, edge[1].index))
                                jointCount += 1

        print ("Flex generated! ({} joints)".format(jointCount))
        return flexibilities


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
            if ((subset == "phi") and (atom.name in phi) and (atom.residue.name != "PRO")):
                return True
            if ((subset == "psi") and (atom.name in psi)):
                return True
            if ((subset == "omega") and (atom.name in omega)):
                return True
            if ((subset == "side") and (atom.name not in ["N", "C"])):
                return True
        return False

    def check_cycle(self, atom1, atom2):
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
        This function will get tripped up by cysteine bonds, since
        those close a large cycle, via the backbone. As such, we check
        if the cycles in the nx.cycle_basis contain ANY SG atoms before
        doing any other check. If it does, we skip it.

        '''
        for cycle in nx.cycle_basis(self._MDTrajObject.topology.to_bondgraph()):
            if ([x.name for x in cycle if x.name=='SG']):
                continue
            if ((atom1 in cycle) and (atom2 in cycle)):
                return True
        return False