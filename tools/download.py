import requests
import os
import re
from tqdm import tqdm
import csv
import json
import subprocess
import random
import pandas as pd
from Bio.PDB import PDBParser, Superimposer
import numpy as np

from rcsbsearchapi.search import AttributeQuery

from openmm.app import *
from omm import *
from openmm.unit import *
import parmed as pmd

def extract_backbone_atoms(structure):
    models = list(structure.get_models())
    backbone_models = []

    for model in models:
        atoms = []
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":  # Exclude hetero/water
                    continue
                for atom_name in ("N", "CA", "C", "O"):
                    if atom_name in residue:
                        atoms.append(residue[atom_name])
        backbone_models.append(atoms)

    return backbone_models

def compute_rmsd_matrix(backbone_models, trim=20):
    n = len(backbone_models)
    rmsd_matrix = np.zeros((n, n))
    super_imposer = Superimposer()

    for i in range(n):
        for j in range(i + 1, n):
            atoms_i = backbone_models[i][trim:-trim]
            atoms_j = backbone_models[j][trim:-trim]

            if len(atoms_i) != len(atoms_j):
                raise ValueError(f"Model {i} and {j} have unequal backbone lengths after trimming")

            super_imposer.set_atoms(atoms_i, atoms_j)
            rms = super_imposer.rms
            rmsd_matrix[i, j] = rms
            rmsd_matrix[j, i] = rms

    return rmsd_matrix

def compute_rmsf(backbone_models, trim=20):
    n_models = len(backbone_models)
    n_atoms = len(backbone_models[0])

    if n_atoms <= 2 * trim:
        raise ValueError("Too few atoms after trimming. Reduce 'trim' value.")

    coords = np.array([[atom.coord for atom in model] for model in backbone_models])  # (models, atoms, 3)
    trimmed_coords = coords[:, trim:-trim, :]  # exclude first and last 20 atoms

    mean_coords = np.mean(trimmed_coords, axis=0)  # (atoms_trimmed, 3)
    fluctuations = trimmed_coords - mean_coords
    squared = np.square(fluctuations)
    rmsf = np.sqrt(np.mean(np.sum(squared, axis=2), axis=0))  # (atoms_trimmed,)

    return rmsf

def calculate_statistics(pdb_filemname, trim=20):
    # Load the PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("multi_model", pdb_filemname)

    # Extract backbone atoms and compute RMSD matrix
    backbone_models = extract_backbone_atoms(structure)
    rmsd_matrix = compute_rmsd_matrix(backbone_models, trim)

    # Exclude diagonal and duplicate entries
    triu_indices = np.triu_indices_from(rmsd_matrix, k=1)
    rmsd_values = rmsd_matrix[triu_indices]

    # RMSD statistics
    mean_rmsd = np.mean(rmsd_values)
    stdev_rmsd = np.std(rmsd_values)
    max_rmsd = np.max(rmsd_values)
    rmsd_span = max_rmsd - np.min(rmsd_values)

    # RMSF statistics
    rmsf = compute_rmsf(backbone_models, trim)
    median = np.median(rmsf)
    stdev = np.std(rmsf)
    threshold = median + stdev

    rmsf_count = np.sum(rmsf > threshold)
    rmsf_fraction = rmsf_count / len(rmsf)

    # Create a dictionary to hold the statistics
    statistics = {
        "mean_rmsd": mean_rmsd,
        "stdev_rmsd": stdev_rmsd,
        "max_rmsd": max_rmsd,
        "rmsd_span": rmsd_span,
        "median_rmsf": median,
        "stdev_rmsf": stdev,
        "rmsf_fraction": rmsf_fraction
    }

    return statistics

def extract_constraints(str_file, distance_filename, dihedral_filename):
    distance_pattern = re.compile(r"""
        ^\s*
        (?P<id>\S+) \s+
        (?P<member_id>\S+) \s+
        (?P<member_logic_code>\S+) \s+
        (?P<assembly_atom_id_1>\S+) \s+
        (?P<entity_assembly_id_1>\S+) \s+
        (?P<entity_id_1>\S+) \s+
        (?P<comp_index_id_1>\S+) \s+
        (?P<seq_id_1>\S+) \s+
        (?P<comp_id_1>[a-zA-Z]{3}) \s+
        (?P<atom_id_1>\S+) \s+
        (?P<atom_type_1>\S+) \s+
        (?P<atom_isotope_number_1>\S+) \s+
        (?P<resonance_id_1>\S+) \s+
        (?P<assembly_atom_id_2>\S+) \s+
        (?P<entity_assembly_id_2>\S+) \s+
        (?P<entity_id_2>\S+) \s+
        (?P<comp_index_id_2>\S+) \s+
        (?P<seq_id_2>\S+) \s+
        (?P<comp_id_2>[a-zA-Z]{3}) \s+
        (?P<atom_id_2>\S+) \s+
        (?P<atom_type_2>\S+) \s+
        (?P<atom_isotope_number_2>\S+) \s+
        (?P<resonance_id_2>\S+) \s+
        (?P<intensity_val>\S+) \s+
        (?P<intensity_lower_val_err>\S+) \s+
        (?P<intensity_upper_val_err>\S+) \s+
        (?P<distance_val>\S+) \s+
        (?P<distance_lower_bound_val>\S+) \s+
        (?P<distance_upper_bound_val>\S+) \s+
        (?P<contribution_fractional_val>\S+) \s+
        (?P<spectral_peak_id>\S+) \s+
        (?P<spectral_peak_list_id>\S+) \s+
        (?P<pdb_record_id_1>\S+) \s+
        (?P<pdb_model_num_1>\S+) \s+
        (?P<pdb_strand_id_1>\S+) \s+
        (?P<pdb_ins_code_1>\S+) \s+
        (?P<pdb_residue_no_1>\S+) \s+
        (?P<pdb_residue_name_1>\S+) \s+
        (?P<pdb_atom_name_1>\S+) \s+
        (?P<pdb_record_id_2>\S+) \s+
        (?P<pdb_model_num_2>\S+) \s+
        (?P<pdb_strand_id_2>\S+) \s+
        (?P<pdb_ins_code_2>\S+) \s+
        (?P<pdb_residue_no_2>\S+) \s+
        (?P<pdb_residue_name_2>\S+) \s+
        (?P<pdb_atom_name_2>\S+) \s+
        (?P<auth_entity_assembly_id_1>\S+) \s+
        (?P<auth_asym_id_1>\S+) \s+
        (?P<auth_chain_id_1>\S+) \s+
        (?P<auth_seq_id_1>\S+) \s+
        (?P<auth_comp_id_1>\S+) \s+
        (?P<auth_atom_id_1>\S+) \s+
        (?P<auth_alt_id_1>\S+) \s+
        (?P<auth_atom_name_1>\S+) \s+
        (?P<auth_entity_assembly_id_2>\S+) \s+
        (?P<auth_asym_id_2>\S+) \s+
        (?P<auth_chain_id_2>\S+) \s+
        (?P<auth_seq_id_2>\S+) \s+
        (?P<auth_comp_id_2>\S+) \s+
        (?P<auth_atom_id_2>\S+) \s+
        (?P<auth_alt_id_2>\S+) \s+
        (?P<auth_atom_name_2>\S+) \s+
        (?P<entry_id>\S+) \s+
        (?P<gen_dist_constraint_list_id>\S*)
    """, re.VERBOSE)
        
    dihedral_pattern = re.compile(r"""
        \s*
        (?P<torsion_angle_constraint_id>\S+) \s+
        (?P<torsion_angle_name>\S+) \s+
        (?P<assembly_atom_id_1>\S+) \s+
        (?P<entity_assembly_id_1>\S+) \s+
        (?P<entity_id_1>\S+) \s+
        (?P<comp_index_id_1>\S+) \s+
        (?P<seq_id_1>\S+) \s+
        (?P<comp_id_1>[a-zA-Z]{3}) \s+
        (?P<atom_id_1>\S+) \s+
        (?P<atom_type_1>\S+) \s+
        (?P<resonance_id_1>\S+) \s+
        (?P<assembly_atom_id_2>\S+) \s+
        (?P<entity_assembly_id_2>\S+) \s+
        (?P<entity_id_2>\S+) \s+
        (?P<comp_index_id_2>\S+) \s+
        (?P<seq_id_2>\S+) \s+
        (?P<comp_id_2>[a-zA-Z]{3}) \s+
        (?P<atom_id_2>\S+) \s+
        (?P<atom_type_2>\S+) \s+
        (?P<resonance_id_2>\S+) \s+
        (?P<assembly_atom_id_3>\S+) \s+
        (?P<entity_assembly_id_3>\S+) \s+
        (?P<entity_id_3>\S+) \s+
        (?P<comp_index_id_3>\S+) \s+
        (?P<seq_id_3>\S+) \s+
        (?P<comp_id_3>[a-zA-Z]{3}) \s+
        (?P<atom_id_3>\S+) \s+
        (?P<atom_type_3>\S+) \s+
        (?P<resonance_id_3>\S+) \s+
        (?P<assembly_atom_id_4>\S+) \s+
        (?P<entity_assembly_id_4>\S+) \s+
        (?P<entity_id_4>\S+) \s+
        (?P<comp_index_id_4>\S+) \s+
        (?P<seq_id_4>\S+) \s+
        (?P<comp_id_4>[a-zA-Z]{3}) \s+
        (?P<atom_id_4>\S+) \s+
        (?P<atom_type_4>\S+) \s+
        (?P<resonance_id_4>\S+) \s+
        (?P<angle_lower_bound_val>\S+) \s+
        (?P<angle_upper_bound_val>\S+) \s+
        (?P<source_experiment_id>\S+) \s+
        (?P<pdb_record_id_1>\S+) \s+
        (?P<pdb_model_num_1>\S+) \s+
        (?P<pdb_strand_id_1>\S+) \s+
        (?P<pdb_ins_code_1>\S+) \s+
        (?P<pdb_residue_no_1>\S+) \s+
        (?P<pdb_residue_name_1>\S+) \s+
        (?P<pdb_atom_name_1>\S+) \s+
        (?P<pdb_record_id_2>\S+) \s+
        (?P<pdb_model_num_2>\S+) \s+
        (?P<pdb_strand_id_2>\S+) \s+
        (?P<pdb_ins_code_2>\S+) \s+
        (?P<pdb_residue_no_2>\S+) \s+
        (?P<pdb_residue_name_2>\S+) \s+
        (?P<pdb_atom_name_2>\S+) \s+
        (?P<pdb_record_id_3>\S+) \s+
        (?P<pdb_model_num_3>\S+) \s+
        (?P<pdb_strand_id_3>\S+) \s+
        (?P<pdb_ins_code_3>\S+) \s+
        (?P<pdb_residue_no_3>\S+) \s+
        (?P<pdb_residue_name_3>\S+) \s+
        (?P<pdb_atom_name_3>\S+) \s+
        (?P<pdb_record_id_4>\S+) \s+
        (?P<pdb_model_num_4>\S+) \s+
        (?P<pdb_strand_id_4>\S+) \s+
        (?P<pdb_ins_code_4>\S+) \s+
        (?P<pdb_residue_no_4>\S+) \s+
        (?P<pdb_residue_name_4>\S+) \s+
        (?P<pdb_atom_name_4>\S+) \s+
        (?P<auth_entity_assembly_id_1>\S+) \s+
        (?P<auth_asym_id_1>\S+) \s+
        (?P<auth_chain_id_1>\S+) \s+
        (?P<auth_seq_id_1>\S+) \s+
        (?P<auth_comp_id_1>\S+) \s+
        (?P<auth_atom_id_1>\S+) \s+
        (?P<auth_alt_id_1>\S+) \s+
        (?P<auth_atom_name_1>\S+) \s+
        (?P<auth_entity_assembly_id_2>\S+) \s+
        (?P<auth_asym_id_2>\S+) \s+
        (?P<auth_chain_id_2>\S+) \s+
        (?P<auth_seq_id_2>\S+) \s+
        (?P<auth_comp_id_2>\S+) \s+
        (?P<auth_atom_id_2>\S+) \s+
        (?P<auth_alt_id_2>\S+) \s+
        (?P<auth_atom_name_2>\S+) \s+
        (?P<auth_entity_assembly_id_3>\S+) \s+
        (?P<auth_asym_id_3>\S+) \s+
        (?P<auth_chain_id_3>\S+) \s+
        (?P<auth_seq_id_3>\S+) \s+
        (?P<auth_comp_id_3>\S+) \s+
        (?P<auth_atom_id_3>\S+) \s+
        (?P<auth_alt_id_3>\S+) \s+
        (?P<auth_atom_name_3>\S+) \s+
        (?P<auth_entity_assembly_id_4>\S+) \s+
        (?P<auth_asym_id_4>\S+) \s+
        (?P<auth_chain_id_4>\S+) \s+
        (?P<auth_seq_id_4>\S+) \s+
        (?P<auth_comp_id_4>\S+) \s+
        (?P<auth_atom_id_4>\S+) \s+
        (?P<auth_alt_id_4>\S+) \s+
        (?P<auth_atom_name_4>\S+) \s+
        (?P<entry_id>\S+) \s+
        (?P<torsion_angle_constraint_list_id>\S*)
    """, re.VERBOSE)

    def prepare_atom_name(original):
        # Q means symmetrical hydrogen group (all hydrogen atom bound to a certain carbon atom)
        if 'Q' in original:
            # this is the case with QD that should be replaced with HD*
            name = original.replace('Q', 'H') + '*'

            # in the previous example, make sure it is boud to CD carbon atom
            carbon = original.replace('Q', 'C') + '*'

            # final result: "name HD* and bonded name CD"
            return name #+ " and bonded name " + carbon

        # methyl group that contains all hydrogen atoms bound to a single carbon atom
        if 'M' in original:
            name_h = original.replace('M', 'H') + '*'
            name_c = original.replace('M', 'C') + '*'
            
            # final result: "name CD or (name HD* and bonded name CD)"
            # return name_c + " or (name " + name_h + " and bonded name " + name_c + ")"
            return name_c + " or (name " + name_h + ")"

        # nothing to change
        return original

    num_distance = 0
    num_dihedral = 0

    with open(str_file, 'r') as fin, open(distance_filename, mode="w", newline="") as distance_file, open(dihedral_filename, mode="w", newline="") as dihedral_file:
        distance_writer = csv.writer(distance_file)
        distance_writer.writerow(["group1", "group2", "avg_dist", "min_dist", "max_dist", "id1", "id2"])

        dihedral_writer = csv.writer(dihedral_file)
        dihedral_writer.writerow(["group1", "group2", "group3", "group4", "min_angle", "max_angle"])

        # parse all lines
        for line in fin:

            # match distances
            distance_match = distance_pattern.match(line)
            if distance_match:
                name1 = prepare_atom_name(distance_match.group("atom_id_1"))
                group1 = "(chainid {} and resid {}) and (name {})".format(
                    distance_match.group("pdb_strand_id_1"), distance_match.group("seq_id_1"), name1
                )

                name2 = prepare_atom_name(distance_match.group("atom_id_2"))
                group2 = "(chainid {} and resid {}) and (name {})".format(
                    distance_match.group("pdb_strand_id_2"), distance_match.group("seq_id_2"), name2
                )

                average_distance = distance_match.group("distance_val")
                min_distance = distance_match.group("distance_lower_bound_val")
                max_distance = distance_match.group("distance_upper_bound_val")

                id1 = distance_match.group("pdb_strand_id_1") + distance_match.group("seq_id_1")
                id2 = distance_match.group("pdb_strand_id_2") + distance_match.group("seq_id_2")

                try:
                    average_distance = float(average_distance)
                    min_distance = float(min_distance)
                    max_distance = float(max_distance)
                except ValueError:
                    return 0, 0
                
                distance_writer.writerow([group1, group2, average_distance, min_distance, max_distance, id1, id2])
                num_distance += 1
                continue

            # match dihedrals
            dihedral_match = dihedral_pattern.match(line)
            if dihedral_match:
                name1 = prepare_atom_name(dihedral_match.group("atom_id_1"))
                group1 = "(chainid {} and resid {}) and (name {})".format(
                    dihedral_match.group("pdb_strand_id_1"), dihedral_match.group("seq_id_1"), name1
                )

                name2 = prepare_atom_name(dihedral_match.group("atom_id_2"))
                group2 = "(chainid {} and resid {}) and (name {})".format(
                    dihedral_match.group("pdb_strand_id_2"), dihedral_match.group("seq_id_2"), name2
                )

                name3 = prepare_atom_name(dihedral_match.group("atom_id_3"))
                group3 = "(chainid {} and resid {}) and (name {})".format(
                    dihedral_match.group("pdb_strand_id_3"), dihedral_match.group("seq_id_3"), name3
                )

                name4 = prepare_atom_name(dihedral_match.group("atom_id_4"))
                group4 = "(chainid {} and resid {}) and (name {})".format(
                    dihedral_match.group("pdb_strand_id_4"), dihedral_match.group("seq_id_4"), name4
                )

                min_angle = dihedral_match.group("angle_lower_bound_val")
                max_angle = dihedral_match.group("angle_upper_bound_val")

                try:
                    min_angle = float(min_angle)
                    max_angle = float(max_angle)
                except ValueError:
                    return 0, 0
                
                dihedral_writer.writerow([group1, group2, group3, group4, min_angle, max_angle])
                num_dihedral += 1

    return num_distance, num_dihedral

def download_file(url, filename):
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for error HTTP statuses

        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if not chunk:
                    break
                f.write(chunk)

        return True

    except:
        return False

def map_resids(pdb_filename, distace_filename):
    """
    Resets residue numbering in a PDB file such that residue numbers start from 0 for each chain.
    
    Args:
        pdb_file (str): Path to the input PDB file.
        output_file (str): Path to the output corrected PDB file.
    """
    residue_map = {}  # To store the mapping between original and new residue identifiers
    new_pdb_lines = []

    with open(pdb_filename, 'r') as infile:
        current_residue = None
        residue_counter = 0
        current_chain = None
        prev_original_residue_number = None
        
        for line in infile:
            # Process ATOM and HETATM records only
            if line.startswith(('ATOM', 'HETATM')):
                residue = line[17:20].strip()  # Residue name
                chain = line[21].strip()  # Chain identifier
                original_residue_number = int(line[22:26].strip())  # Original residue number

                # Reset counter when encountering a new chain
                if chain != current_chain:
                    current_chain = chain
                    current_residue = None
                    residue_counter = 0

                if residue != current_residue or original_residue_number != prev_original_residue_number:
                    current_residue = residue
                    residue_counter += 1
                
                # Update mapping
                residue_map[f"{chain}{residue_counter}"] = original_residue_number
                prev_original_residue_number = original_residue_number
                
                # Replace residue number in the PDB line
                line = f"{line[:22]}{residue_counter:>4}{line[26:]}"
                
            # Write the modified or unmodified line to output
            new_pdb_lines.append(line)

    # # print the residue map
    # for k, v in residue_map.items():
    #     print(f"{k} -> {v}")

    # Write the new PDB file
    with open(pdb_filename, 'w') as outfile:
        outfile.writelines(new_pdb_lines)

    # Read the dataframe and update residue numbers
    df = pd.read_csv(distace_filename)

    # Create two new vectors: 'resid1' and 'resid2'. they will contain the new residue numbers as 'chain1' and 'resid1' with residue_map['chain1resid1']
    df['resid1'] = df.apply(lambda row: residue_map[row['id1']], axis=1)
    df['resid2'] = df.apply(lambda row: residue_map[row['id2']], axis=1)

    # Remove the old columns
    df.drop(columns=['id1', 'id2'], inplace=True)

    # Save the updated dataframe
    df.to_csv(distace_filename, index=False)

def minimize_and_save(prmtop_file, rst7_file, output_rst7_file, inpcrd_min_file):
    # Load the prmtop and rst7 files
    prmtop = AmberPrmtopFile(prmtop_file)
    inpcrd = AmberInpcrdFile(rst7_file)

    # Create a simulation system
    # @TODO Add implicit solvent model - same salt conc as in robosample
    system = prmtop.createSystem(
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=1.0*nanometers,
        constraints=HBonds,
        implicitSolvent=OBC2
    )

    # Set up an integrator
    integrator = LangevinIntegrator(
        300*kelvin,      # Temperature (dummy, not needed for minimization)
        1.0/picoseconds, # Friction coefficient (dummy)
        0.002*picoseconds # Time step (dummy)
    )

    # Create a simulation object with CPU platform
    platform = Platform.getPlatformByName('CPU')
    simulation = Simulation(prmtop.topology, system, integrator, platform)

    # Set the positions from the rst7 file
    simulation.context.setPositions(inpcrd.positions)

    # If box vectors are available, set them as well
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    # Perform energy minimization
    print("Minimizing energy...")
    simulation.minimizeEnergy()
    print("Energy minimization complete.")

    # Extract minimized positions and box vectors
    state = simulation.context.getState(getPositions=True)
    minimized_positions = state.getPositions()
    # box_vectors = state.getPeriodicBoxVectors()

    # Use ParmEd to write the minimized structure to an rst7 file
    structure = pmd.load_file(prmtop_file, rst7_file)
    structure.positions = minimized_positions
    # if box_vectors is not None:
    #     structure.box = [
    #         box_vectors[0][0] / angstroms, box_vectors[1][1] / angstroms, box_vectors[2][2] / angstroms,
    #         box_vectors[0][1] / angstroms, box_vectors[0][2] / angstroms, box_vectors[1][2] / angstroms
    #     ]

    structure.save(output_rst7_file, format="rst7", overwrite=True)
    structure.save(inpcrd_min_file, format="rst7", overwrite=True)

    print(f"Minimized structure saved to {output_rst7_file}")

if __name__ == '__main__':
    basedir = 'data-raw-optimized'
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    with open('cath-b-newest-names', 'r') as f:
        content = f.read()
    pattern = r'^(\d+\.\d+)\s+(.*)$'
    regex = re.compile(pattern, re.MULTILINE)
    matches = regex.findall(content)

    # Hold the data11
    pdbids = []
    cath_lineages = []
    cath_architecture_names = []

    # for match in tqdm(matches, desc="Processing CATH matches"):
    #     cath_lineage = match[0] # Matches class and architecture
    #     name = match[1] # Name of the architecture
            
    #     # # Search for this architecture in the RCSB database
    #     # # In 3.40.50.300, 3 stand for class, 40 for architecture, 50 for topology (fold) and 300 for homologous superfamily
    #     q1 = AttributeQuery("rcsb_polymer_instance_annotation.annotation_lineage.id", "exact_match", cath_lineage)
    #     q2 = AttributeQuery("rcsb_polymer_instance_annotation.type", "exact_match", "CATH")

    #     # NMR methods
    #     q3 = AttributeQuery("exptl.method", "exact_match", "SOLUTION NMR")
    #     q4 = AttributeQuery("exptl.method", "exact_match", "SOLID-STATE NMR")
    #     q_nmr = q3.or_(q4)

    #     # Get single chain proteins
    #     # This enforces the number of modeled polymer chains in the asymmetric unit to be exactly one
    #     q5 = AttributeQuery("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 1)

    #     # Protein lenght is at least 50 residues
    #     q6 = AttributeQuery("entity_poly.rcsb_sample_sequence_length", "greater_or_equal", 50)

    #     # Length and quality filters
    #     q7 = AttributeQuery("rcsb_assembly_info.unmodeled_polymer_monomer_count", "equals", 0)
    #     q8 = AttributeQuery("rcsb_entry_info.deposited_nonpolymer_entity_instance_count", "equals", 0)

    #     # Combine all queries
    #     query = q1.and_(q2).and_(q_nmr).and_(q5).and_(q6).and_(q7).and_(q8)
            
    #     # get all pdbids
    #     all_pdbids = list(query())
    #     pdbids += all_pdbids
    #     cath_lineages = [cath_lineage] * len(all_pdbids)
    #     cath_architecture_names = [name] * len(all_pdbids)

    #     break

    pdbids = ['1APQ', '1AF8']
    cath_lineages = ['1.1', '2.1']
    cath_architecture_names = ['ana', 'cimpanzini bananini']

    df = pd.DataFrame(columns=['pdbid', 'cath_lineage', 'cath_architecture_name', 'mean_rmsd', 'stdev_rmsd', 'max_rmsd', 'rmsd_span', 'median_rmsf', 'stdev_rmsf', 'rmsf_fraction'])

    # Use tqdm and enumerate to traverse all pdbids
    # pdbids = ['1APQ']
    # cath_architecture_names = ['CATH Architecture Example'] * len(pdbids)  # Dummy names for the example
    #     cath_lineages = [cath_lineage] * len(all_pdbids)
    for i, (pdbid, cath_lineage, cath_architecture_name) in enumerate(tqdm(zip(pdbids, cath_lineages, cath_architecture_names), desc="Processing PDB IDs")):

        # if pdbid != '1AF8':
        #     continue

        mr_url = f'https://files.rcsb.org/download/{pdbid}_mr.str'
        mr_filename = f'{basedir}/{pdbid}_mr.str'

        distance_filename = f'{basedir}/{pdbid}.distance.csv'
        dihedral_filename = f'{basedir}/{pdbid}.dihedral.csv'

        pdb_url = f'https://files.rcsb.org/download/{pdbid}.pdb'
        pdb_filename = f'{basedir}/{pdbid}.pdb'
        pdb_protein_only_filename = f'{basedir}/{pdbid}.protein_only.pdb'

        pdb4amber_nonprot_filename = f'{basedir}/{pdbid}.protein_only_nonprot.pdb'
        pdb4amber_renum_filename = f'{basedir}/{pdbid}.protein_only_renum.txt'
        pdb4amber_sslink_filename = f'{basedir}/{pdbid}.protein_only_sslink'

        prmtop_file = f'{basedir}/{pdbid}.prmtop'
        rst7_file = f'{basedir}/{pdbid}.rst7'
        rst7_min_file = f'{basedir}/{pdbid}_min.rst7'
        inpcrd_min_file = f'{basedir}/{pdbid}_min.inpcrd'

        def cleanup(all=False):
            # TODO rst7_file inpcrd_min_file
            remove = [mr_filename, pdb_filename, pdb4amber_nonprot_filename, pdb4amber_renum_filename, pdb4amber_sslink_filename, 'tleap.txt']
            if all:
                remove += [distance_filename, dihedral_filename, pdb_protein_only_filename, prmtop_file, rst7_min_file, inpcrd_min_file]
            [os.remove(file) for file in remove if os.path.exists(file)]
            
        # def do():
        try:
            # Download the .mr.str file
            download_file(mr_url, mr_filename)
            if not os.path.exists(mr_filename):
                raise ValueError("No mr file")

            # Check if it has distance and dihedral constraints
            num_distance, num_dihedral = extract_constraints(mr_filename, distance_filename, dihedral_filename)
            if num_distance == 0 or num_dihedral == 0:
                raise ValueError("No distance or dihedral constraints")
                
            # Download the pdb file
            download_file(pdb_url, pdb_filename)

            # Perform statistics here
            stats = calculate_statistics(pdb_filename)

            # Prepare for amber
            pdb4amber = f"pdb4amber -i {pdb_filename} -o {pdb_protein_only_filename}"
            subprocess.run(pdb4amber, shell=True, check=True, text=True, capture_output=False)

            # check if the pdb file contains HETATM records
            with open(pdb_protein_only_filename, 'r') as f:
                content = f.read()
                if 'HETATM' in content:
                    raise ValueError("HETATM records found")
                    
            # Write tleap to file.txt
            with open("tleap.txt", "w") as file:
                tleap = f"""
source leaprc.protein.ff19SB
protein = loadpdb {pdb_protein_only_filename}
saveAmberParm protein {prmtop_file} {rst7_file}
savepdb protein {pdb_protein_only_filename}
quit
        """
                file.write(tleap.strip())

            # Execute tleap with the file
            subprocess.run("tleap -f tleap.txt", shell=True, check=True, text=True, capture_output=False)

            # Minimize the structure
            minimize_and_save(prmtop_file, rst7_file, rst7_min_file, inpcrd_min_file)

            # # Renumber the pdb file from A1, A2, A3, ..., An, Bn+1, Bn+2, ..., Bn+m to A1, A2, A3, ..., An, B1, B2, ..., Bm
            # # This processes only ATOM records and ignoring HETATM and TER records
            # # TODO WARNING THIS DOES NOT TOUCH THE DISTANCE_FILE NOW
            # # WARNING THIS DOES NOT WORK WITH MULTIPLE CHAINS
            # map_resids(pdb_protein_only_filename, distance_filename)

            # Delete all temporary files
            cleanup()

            new_row = pd.DataFrame([{
                'pdbid': pdbid,
                'cath_lineage': cath_lineage,
                'cath_architecture_name': cath_architecture_name,
                **stats
            }])

            df = pd.concat([df, new_row], ignore_index=True)
        except Exception as e:
            cleanup(True)
            print(f"Error processing {pdbid}: {e}")
            continue

        # do()

        # # Do only one for now
        # break

    # print(df)

    # save df to file
    df.to_csv(f'{basedir}/statistics.csv', index=False)
