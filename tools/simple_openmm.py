import openmm as mm
import openmm.app as app
import openmm.unit as unit

# 1. Create a tiny system (2 atoms)
system = mm.System()
system.addParticle(12.0 * unit.amu) # Atom 0 (Carbon-ish)
system.addParticle(12.0 * unit.amu) # Atom 1

# 2. Add a simple harmonic bond
force = mm.HarmonicBondForce()
# bond between atom 0 and 1, length 1.0 nm, k=100 kJ/mol/nm^2
force.addBond(0, 1, 1.0 * unit.nanometer, 100.0 * unit.kilojoules_per_mole/unit.nanometer**2)
system.addForce(force)

# 3. Set positions
positions = [[0, 0, 0], [1.1, 0, 0]] * unit.nanometer # Stretched by 0.1nm

# 4. Initialize Integrator and Platform
integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)

try:
    # Try to use OpenCL first
    platform = mm.Platform.getPlatformByName('OpenCL')
    # If you want to force it to use a specific device index:
    # properties = {'OpenCLPrecision': 'single', 'OpenCLDeviceIndex': '0'}
    context = mm.Context(system, integrator, platform)
except:
    # Fallback to Reference
    print("OpenCL failed or not found. Falling back to Reference platform.")
    platform = mm.Platform.getPlatformByName('Reference')
    context = mm.Context(system, integrator, platform)

context.setPositions(positions)

# 5. Calculate and Print Energy
state = context.getState(getEnergy=True)
print(f"Using Platform: {context.getPlatform().getName()}")
print(f"Potential Energy: {state.getPotentialEnergy()}")

# 6. Clean up
del context
