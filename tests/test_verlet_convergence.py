import robosample
import numpy as np
import mdtraj as md
import pytest

name = '1APQ.test.rigid'
seed = 42
prmtop = 'examples/1APQ.prmtop'
xyz = 'examples/1APQ.rst7'
write_freq = 1
equil_steps = 0
prod_steps = 2

context = robosample.Context(name=name, seed=seed, prmtop=prmtop, inpcrd=xyz, write_freq=write_freq, testing=True)

# Put torsions on middle bonds in standard dihedrals and integrate with 10 fs time step for 1 ps
# Torsions are protein phi, psi and chi1
sele = context.getDefaultBonds('standard')
context.addTorsionalWorld(sele).addSampler(timeStep=0.01, mdSteps=2, boostMDSteps=2)

# Add one replica at 300 K
context.initialize([300.0])

for num_steps in range [1, 2, 4, 8]:
    time_step = 0.001 # 1 fs
    violations = robosample.RigidBodyViolations()
    context.getWorld(0).has_rigid_body_violations(time_step, num_steps, violations)

    print(violations)
