import openmm
import openmm.app as app
import openmm.unit as unit
from sys import stdout
import argparse

def minimize_energy(input_pdb, output_pdb):
    # Load PDB file
    pdb = app.PDBFile(input_pdb)
    
    # Create force field
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    # Create system from topology and force field
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    
    # Set up simulation environment
    integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    simulation = app.Simulation(pdb.topology, system, integrator)
    
    # Set initial positions
    simulation.context.setPositions(pdb.positions)
    
    # Minimize energy
    print("Minimizing energy...")
    simulation.minimizeEnergy()
    
    # Get minimized positions
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    
    # Save minimized structure to output PDB file
    with open(output_pdb, 'w') as f:
        app.PDBFile.writeFile(simulation.topology, positions, f)
    
    print(f"Energy minimization completed. Output saved to {output_pdb}.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Energy minimize a PDB file using OpenMM.')
    parser.add_argument('--input', required=True, type=str, help='Input PDB file path.')
    parser.add_argument('--output', required=True, type=str, help='Output PDB file path.')
    
    args = parser.parse_args()
    
    minimize_energy(args.input, args.output)
