# Molecule information
MOLECULES ala10                  # Molecules to load
PRMTOP ligand.prmtop       # Parameter file
INPCRD ligand.rst7      # Coordinate / Restart file
RBFILE ligand.rb    # Rigid bodies definition file 
FLEXFILE ligand.all.flex  # Flexible bonds definition file
ROOT_MOBILITY Cartesian  # Ground to Compound mobilizer
OUTPUT_DIR temp # Output directory
MEMBRANE 0 0 0 0 # Membrane x y z width and resolution
CONTACTS -1  # Membrane-atom contact
CONSTRAINTS -1  # Constrained atoms

# Simulation
RUN_TYPE Normal # normal HMC or Non-Eq HMC
ROUNDS 1                   # Gibbs sampling rounds
ROUNDS_TILL_REBLOCK 10 
RANDOM_WORLD_ORDER FALSE 
WORLDS R0 # Worlds unique names with 2 letters 

ROOTS 0 # Atoms representing the first body
SAMPLER VV      # Sampler
TIMESTEPS 0.001              # Timesteps 
MDSTEPS 1     # Number of MD trial steps
BOOST_MDSTEPS 1   
SAMPLES_PER_ROUND 1  # Number of acc-rej steps within a Gibbs round
REPRODUCIBLE TRUE
SEED 9999


ROUNDS_TIL_REBLOCK 1

# Thermodynamics
THERMOSTAT Andersen     # Thermostat 
TEMPERATURE_INI  300    # Initial temperature
TEMPERATURE_FIN  300    # Final temperature
BOOST_TEMPERATURE   1     # Guidance Hamiltonian temperature
FFSCALE AMBER        # Force field scale factors
GBSA 0.0         # GBSA scale factor

# Generalized coordinates related
FIXMAN_POTENTIAL FALSE  # Use Fixman potential
FIXMAN_TORQUE FALSE        # Use Fixman torque

# Output
#VISUAL TRUE TRUE              # Use the visualizer
VISUAL FALSE           # Use the visualizer
PRINT_FREQ 1 # Output frequency
WRITEPDBS 0    # Write pdbs frequency
GEOMETRY FALSE             # Calculate geometric features

DISTANCE 2 17
DIHEDRAL 16 14 0 6 14 0 6 8

# Software specs
THREADS 0
OPENMM TRUE
OPENMM_CalcOnlyNonbonded FALSE

NONBONDED_METHOD    0
NONBONDED_CUTOFF    0
