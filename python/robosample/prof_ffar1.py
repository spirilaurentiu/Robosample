"""Dedicated, self-contained profiling workload for FFAR1.

One robotic (internal-coordinate) world sampling the phi/psi backbone torsions
of the target GPCR at a 5 fs timestep, on top of Free (6 external DOF) root
bodies for the target and the lipid/solvent molecules. This is deliberately a
*minimal* proxy for `run_ffar1.py`'s six-world production setup: it isolates the
ABA recursion spine + the OpenMM per-atom force bridge that the `optimizer`
agent cares about, without the noise of five extra worlds.

    python3 python/robosample/prof_ffar1.py ffar1 examples/ffar1.prmtop examples/ffar1.rst7 6000 0 100 20

Invoked by `nox -s profile`. Root-mobility layout mirrors run_ffar1.py so the
body count (and thus the spine load) is representative of the real system:
  molecule 0        -> Free  (target GPCR external DOFs)
  molecules 1, 2    -> Weld  (outer proteins, rigid)
  molecules 3..N-1  -> Free  (lipids + solvent, each a free rigid body)
"""

import argparse
import time

import robosample

parser = argparse.ArgumentParser(description="FFAR1 profiling workload.")
parser.add_argument("name", type=str, help="Name of the simulation.")
parser.add_argument("prmtop", type=str, help="Relative path to the .prmtop file.")
parser.add_argument("inpcrd", type=str, help="Relative path to the .inpcrd file.")
parser.add_argument("seed", type=int, help="The seed.")
parser.add_argument("equil_steps", type=int, help="The number of equilibration steps.")
parser.add_argument("prod_steps", type=int, help="The number of production steps.")
parser.add_argument("write_freq", type=int, help="CSV and DCD write frequency.")
args = parser.parse_args()

# Single replica at 300 K: profiling wants steady-state force evaluation, not a
# temperature ladder.
T0 = 300.0

context = robosample.Context(
    name=args.name,
    seed=args.seed,
    prmtop=args.prmtop,
    inpcrd=args.inpcrd,
    write_freq=args.write_freq,
    testing=False,
)

# ---- External free DOFs (root mobility) -------------------------------------
# Mirror run_ffar1.py so the profiled body count reflects the real system.
context.set_root_mobility(0, robosample.rb.RootMobility.Free)  # target GPCR
context.set_root_mobility(1, robosample.rb.RootMobility.Weld)  # outer protein
context.set_root_mobility(2, robosample.rb.RootMobility.Weld)  # outer protein
for mol_ix in range(3, context.get_num_molecules()):
    context.set_root_mobility(mol_ix, robosample.rb.RootMobility.Free)  # lipids/solvent

# ---- One robotic world: phi/psi backbone torsions of the target at 5 fs ------
bonds = context.standard_dihedral_bonds.loc[
    context.standard_dihedral_bonds["dihedral_type"].isin(["phi", "psi"])
    & (context.standard_dihedral_bonds["molecule_index"] == 0)
]
sele = context.build_flexibilities(bonds)
context.add_robotic_world(sele).add_sampler(
    timeStep=0.005,  # 5 fs
    mdSteps=500,
    boostMDSteps=500,
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
)

context.initialize([T0])

start_time = time.perf_counter()
context.run_rex(args.equil_steps, args.prod_steps, args.write_freq, True)
duration = time.perf_counter() - start_time

# Parsed by nox -s profile into MANIFEST.md; keep this exact wording.
print(f"run_rex() took {duration:.4f} seconds")
