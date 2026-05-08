import argparse

import openmm as mm
from openmm import app, unit

# python3 python/robosample/simulate_openmm.py ffar1 examples/ffar1.prmtop examples/ffar1.rst7 6000

# Create the parser
parser = argparse.ArgumentParser(description="Process PDB code and seed.")

# Add the arguments
parser.add_argument("name", type=str, help="Name of the simulation.")
parser.add_argument("prmtop", type=str, help="Relative path to the .prmtop file.")
parser.add_argument("inpcrd", type=str, help="Relative path to the .inpcrd file.")
parser.add_argument("seed", type=int, help="The seed.")

# Parse the arguments
args = parser.parse_args()

if __name__ == "__main__":
    prmtop = app.AmberPrmtopFile(args.prmtop)
    inpcrd = app.AmberInpcrdFile(args.inpcrd)
    system = prmtop.createSystem(
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1.2 * unit.nanometer,
        implicitSolvent=app.OBC2,
        constraints=None,
    )

    for force in system.getForces():
        if isinstance(force, mm.HarmonicBondForce):
            force.setForceGroup(0)
        elif isinstance(
            force, (mm.HarmonicAngleForce, mm.PeriodicTorsionForce, mm.RBTorsionForce)
        ):
            force.setForceGroup(1)
        else:
            # NonbondedForce, GBSAOBCForce, CustomGBForce, CMMotionRemover
            force.setForceGroup(2)

    thermostat = mm.AndersenThermostat(300 * unit.kelvin, 1.0 / unit.picosecond)
    thermostat.setForceGroup(2)  # couple to the slowest group
    system.addForce(thermostat)

    timestep = 2 * unit.femtoseconds  # outer step
    integrator = mm.MTSIntegrator(
        timestep,
        [
            (2, 1),  # group 2 once per outer step  → 2 fs effective
            (1, 2),  # group 1 twice per group-2 step → 1 fs effective
            (0, 4),
        ],  # group 0 four times per group-1 step → 0.25 fs effective
    )

    platform = mm.Platform.getPlatformByName("CUDA")
    simulation = app.Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)

    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

    # Reporters
    steps_per_frame = 5_000
    output_dcd = f"{args.name}_omm_{args.seed}.dcd"
    simulation.reporters.append(app.DCDReporter(output_dcd, steps_per_frame))

    simulation.runForClockTime(2.5 * unit.hour)

# python3 python/robosample/simulate_openmm.py ffar1 examples/ffar1.prmtop examples/ffar1.rst7 6000
# python3 python/robosample/simulate_openmm.py ffar1 examples/ffar1.prmtop examples/ffar1.rst7 6001
# python3 python/robosample/simulate_openmm.py ffar1 examples/ffar1.prmtop examples/ffar1.rst7 6002
# python3 python/robosample/simulate_openmm.py ffar1 examples/ffar1.prmtop examples/ffar1.rst7 6003
# python3 python/robosample/simulate_openmm.py ffar1 examples/ffar1.prmtop examples/ffar1.rst7 6004
