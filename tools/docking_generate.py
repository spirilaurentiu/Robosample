import argparse
import mdtraj as md
import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
import tempfile
import os

def sample_points_around_convex_hull(points, num_samples, distance):
    # Create the convex hull
    hull = ConvexHull(points)
    
    # Get the triangles (simplices) of the convex hull
    triangles = hull.simplices
    
    # Calculate the areas of the triangles
    areas = np.zeros(len(triangles))
    for i, simplex in enumerate(triangles):
        a, b, c = points[simplex]
        areas[i] = np.linalg.norm(np.cross(b-a, c-a)) / 2
    
    # Normalize the areas to get probabilities
    probs = areas / areas.sum()
    
    # Sample triangles based on their areas
    chosen_triangles = np.random.choice(len(triangles), num_samples, p=probs)
    
    sampled_points = []
    for triangle_idx in chosen_triangles:
        # Get the vertices of the chosen triangle
        a, b, c = points[triangles[triangle_idx]]
        
        # Generate random barycentric coordinates
        r1, r2 = np.random.random(2)
        r1 = np.sqrt(r1)
        barycentric_coords = [1 - r1, r1 * (1 - r2), r1 * r2]
        
        # Calculate the point on the triangle
        point = barycentric_coords[0] * a + barycentric_coords[1] * b + barycentric_coords[2] * c
        
        # Calculate the normal vector of the triangle
        normal = np.cross(b - a, c - a)
        normal = normal / np.linalg.norm(normal)
        
        # Ensure the normal points outward
        centroid = np.mean(points, axis=0)
        if np.dot(point - centroid, normal) < 0:
            normal = -normal
        
        # Move the point outward by the specified distance
        point += distance * normal
        
        sampled_points.append(point)
    
    return np.array(sampled_points)

def translate_mol2(input_mol2_file, point, output_mol2_file):
    # Load the mol2 file
    mol = Chem.MolFromMol2File(input_mol2_file, removeHs=False)
    
    # Get the conformer
    conf = mol.GetConformer()
    
    # Calculate the centroid of all atom positions
    centroid = np.zeros(3)
    num_atoms = mol.GetNumAtoms()
    for atom_idx in range(num_atoms):
        pos = conf.GetAtomPosition(atom_idx)
        centroid += np.array([pos.x, pos.y, pos.z])
    centroid /= num_atoms

    # Center the atoms to the origin
    for atom_idx in range(num_atoms):
        pos = conf.GetAtomPosition(atom_idx)
        centered_pos = np.array([pos.x, pos.y, pos.z]) - centroid
        conf.SetAtomPosition(atom_idx, centered_pos)

    # Add the point to all atom positions
    for atom_idx in range(num_atoms):
        pos = conf.GetAtomPosition(atom_idx)
        new_pos = np.array([pos.x + point[0], pos.y + point[1], pos.z + point[2]])
        conf.SetAtomPosition(atom_idx, new_pos)
    
    # Save the modified molecule to a new mol file
    temp_mol_file = 'temp.mol'
    Chem.MolToMolFile(mol, temp_mol_file)

    # Convert the mol file to mol2
    subprocess.run(['obabel', '-imol', temp_mol_file, '-O', output_mol2_file])
    os.remove(temp_mol_file)


def main():
    # parser = argparse.ArgumentParser(description="Process PDB and MOL2 files with a specified distance.")
    # parser.add_argument('pdb_file', type=str, help='Path to the PDB file (protein)')
    # parser.add_argument('mol2_file', type=str, help='Path to the MOL2 file (ligand)')
    # parser.add_argument('distance', type=float, help='Distance in nanometers')

    # args = parser.parse_args()

    # pdb_file = args.pdb_file
    # distance = args.distance

    pdb_file = 'rage.noh.pdb'
    mol2_file = 'fps.mol2'
    distance = 1

    traj = md.load_pdb(pdb_file)

    pts = traj.xyz[0]
    hull = ConvexHull(pts)

    # Sample random points on the surface of the convex hull and project them outward
    num_points = 10  # Number of points to sample
    sampled_points = sample_points_around_convex_hull(pts, num_points, distance)

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection="3d")

    # # Plot defining corner points
    # ax.plot(pts.T[0], pts.T[1], pts.T[2], "ko")

    # # Plot sampled points
    # ax.plot(sampled_points.T[0], sampled_points.T[1], sampled_points.T[2], "bo")

    # # Plot the convex hull
    # for s in hull.simplices:
    #     s = np.append(s, s[0])  # Here we cycle back to the first coordinate
    #     ax.plot(pts[s, 0], pts[s, 1], pts[s, 2], "r-")

    # # Make axis label
    # for i in ["x", "y", "z"]:
    #     eval("ax.set_{:s}label('{:s}')".format(i, i))

    # plt.show()

    mol = Chem.MolFromMolFile(mol2_file)
    for i, point in enumerate(sampled_points):
        translated_mol2 = f'translated_{i}.mol2'
        amber_mol2 = f'amber_{i}.mol2'
        prepi_file = f'fps_{i}.prepi'
        frcmod_file = f'fps_{i}.frcmod'
        prmtop = f'rage.fps.{i}.prmtop'
        rst7 = f'rage.fps.{i}.rst7'

        # convert point around protein from nm to angstrom and translate the ligand
        translate_mol2(mol2_file, point * 10, translated_mol2)

        # prepare amber input
        subprocess.run(['antechamber', '-i', translated_mol2, '-fi', 'mol2', '-o', amber_mol2, '-fo', 'mol2', '-c', 'gas', '-at', 'amber', '-pf', 'y'])
        subprocess.run(['antechamber', '-i', amber_mol2, '-fi', 'mol2', '-o', prepi_file, '-fo', 'prepi', '-c', 'gas', '-at', 'amber', '-pf', 'y'])
        subprocess.run(['parmchk2', '-i', amber_mol2, '-f', 'mol2', '-o', frcmod_file])

        tleap_input = f"""
        source leaprc.protein.ff19SB
        source leaprc.gaff
        loadamberprep {prepi_file}
        loadamberparams {frcmod_file}

        protein = loadpdb {pdb_file}
        ligand = loadmol2 {amber_mol2}
        complex = combine {{protein ligand}}

        saveamberparm complex {prmtop} {rst7}
        quit
        """
        # Write tleap_input to a temporary file
        with tempfile.NamedTemporaryFile('w', delete=False) as tleap_file:
            tleap_file.write(tleap_input)
            tleap_file_path = tleap_file.name

        # Execute tleap with the input file
        subprocess.run(['tleap', '-f', tleap_file_path])

        # delete temporary files with python
        os.remove(translated_mol2)
        os.remove(amber_mol2)
        os.remove(prepi_file)
        os.remove(frcmod_file)
        os.remove('leap.log')


if __name__ == "__main__":
    main()
