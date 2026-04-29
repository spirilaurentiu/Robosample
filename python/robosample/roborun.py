import argparse
import time

import robosample

# python3 roborun.py 2but ../examples/2but.prmtop ../examples/2but.rst7 6000 0 10 1

# python3 roborun.py 2ala_test ./data-raw/2ala.prmtop ./data-raw/2ala.inpcrd 6000 1 1 1
# python3 roborun.py 2ala.2ala.test ./data-raw/2ala.2ala.prmtop ./data-raw/2ala.2ala.inpcrd 6000 0 10 1

# python3 roborun.py 2ala_test ./data-raw/2ala.prmtop ./data-raw/2ala.inpcrd 6000 0 1000 1
# vmd -f data-raw/2ala.prmtop 2ala_test_6000.repl0.dcd

# python3 roborun.py gfcd ./data-raw/GfcD121A-1min.prmtop ./data-raw/GfcD121A-1min.inpcrd 6000 0 10 1
# python3 roborun.py gfcd ./data-raw/FFAR1_nanodisc.AMBER.prmtop ./data-raw/FFAR1_nanodisc.AMBER.min.rst7 6000 0 10 1

# python3 roborun.py 2KTA_TEST ./data-raw/2KTA.prmtop ./data-raw/2KTA_min.inpcrd 6000 0 10 1
# python3 roborun.py 2KQ8_TEST ./data-raw/2KQ8.prmtop ./data-raw/2KQ8_min.inpcrd 6000 0 10 1
# python3 roborun.py 2HHI_TEST ./data-raw/2HHI.prmtop ./data-raw/2HHI_min.inpcrd 6000 0 10 1
# python3 roborun.py 2KNQ_TEST ./data-raw/2KNQ.prmtop ./data-raw/2KNQ_min.inpcrd 6000 0 10 1
# python3 roborun.py 2EZA_TEST ./data-raw/2EZA.prmtop ./data-raw/2EZA_min.inpcrd 6000 0 10 1
# python3 roborun.py 2LGJ_TEST ./data-raw/2LGJ.prmtop ./data-raw/2LGJ_min.inpcrd 6000 0 10 1
# python3 roborun.py 2KCU_TEST ./data-raw/2KCU.prmtop ./data-raw/2KCU_min.inpcrd 6000 0 10 1
# python3 roborun.py 2LT2_TEST ./data-raw/2LT2.prmtop ./data-raw/2LT2_min.inpcrd 6000 0 10 1
# python3 roborun.py 2KPU_TEST ./data-raw/2KPU.prmtop ./data-raw/2KPU_min.inpcrd 6000 0 10 1
# python3 roborun.py 2JOZ_TEST ./data-raw/2JOZ.prmtop ./data-raw/2JOZ_min.inpcrd 6000 0 10 1
# python3 roborun.py 1L1I_TEST ./data-raw/1L1I.prmtop ./data-raw/1L1I_min.inpcrd 6000 0 10 1
# python3 roborun.py 1BLA_TEST ./data-raw/1BLA.prmtop ./data-raw/1BLA_min.inpcrd 6000 0 10 1
# python3 roborun.py 2L3B_TEST ./data-raw/2L3B.prmtop ./data-raw/2L3B_min.inpcrd 6000 0 10 1
# python3 roborun.py 2JR6_TEST ./data-raw/2JR6.prmtop ./data-raw/2JR6_min.inpcrd 6000 0 10 1
# python3 roborun.py 1APQ_TEST ./data-raw/1APQ.prmtop ./data-raw/1APQ_min.inpcrd 6000 0 10 1
# python3 roborun.py 1A5E_TEST ./data-raw/1A5E.prmtop ./data-raw/1A5E_min.inpcrd 6000 0 10 1

# python3 python/robosample/roborun.py ala-dipeptide examples/ala-dipeptide.prmtop examples/ala-dipeptide.rst7 6000 0 10 1
# python3 python/robosample/roborun.py ala-dipeptide examples/ala-dipeptide.prmtop examples/ala-dipeptide.rst7 6000 0 500 1

# perf record -e cycles:u -j any,u --call-graph fp --no-bpf-event -o perf.data -- python3 python/robosample/roborun.py 1APQ examples/1APQ.prmtop examples/1APQ.rst7 6000 0 100 1
# perf inject -j -i perf.data -o perf.lbr.data
# create_gcov --binary=python/robosample/robo_bindings.cpython-312-x86_64-linux-gnu.so --profile=perf.lbr.data --gcov=autofdo.afdo


# perf stat -e cycles,instructions,branches,branch-misses,cache-misses,task-clock,context-switches,cpu-migrations python3 python/robosample/roborun.py 1apq examples/1APQ.prmtop examples/1APQ.rst7 6000 0 100 1


# PYTHONPERFSUPPORT=1 perf record -o profile-data/perf.data -g --symbol-filter='Context::RunREX(int, int)' -e cycles:u -j any,u -- python3 python/robosample/roborun.py 1apq examples/1APQ.prmtop examples/1APQ.rst7 6000 0 100 1
# hotspot profile-data/perf.data

"""
# actual profiling
PYTHONPERFSUPPORT=1 \
  perf record \
    -F 999 \
    -e cycles:u \
    -j any,u \
    -g \
    --call-graph fp \
    -o profile-data/perf.1APQ.cycles.data \
python3 python/robosample/roborun.py 1apq examples/1APQ.prmtop examples/1APQ.rst7 6000 0 500 20


python3 python/robosample/roborun.py 1apq examples/1APQ.prmtop examples/1APQ.rst7 6000 0 10 1






perf record -o profile-data/perf.data -e cycles:u -j any,u -- python3 python/robosample/roborun.py GfcD examples/GfcDstrippedMin.prmtop examples/GfcDstrippedMin.rst7 6000 0 20 1

# this has no file output, it just prints some text to the terminal

perf stat -e cycles,instructions,branches,branch-misses \
-e L1-dcache-loads,L1-dcache-load-misses \
-e l2_cache_req_stat.all,l2_cache_req_stat.dc_access_in_l2,l2_cache_req_stat.dc_hit_in_l2 \
-e LLC-loads,LLC-load-misses \
-e dTLB-loads,dTLB-load-misses \
-e fp_ret_sse_avx_ops.all,node-load-misses,node-stores \
-e task-clock,context-switches,page-faults \
python3 python/robosample/roborun.py 1apq examples/1APQ.prmtop examples/1APQ.rst7 6000 0 100 1



perf record -o profile-data/perf.1APQ.cache.data \
  -e cycles:u -j any,u \
  -e mem_load_retired.l1_miss:upp \
  -e mem_load_retired.l2_miss:upp \
  -e mem_load_retired.l3_miss:upp \
  -j any,u -g \
python3 python/robosample/roborun.py 1apq examples/1APQ.prmtop examples/1APQ.rst7 6000 0 100 1



# this generates some sort of file
nsys profile --trace=cuda,osrt python3 python/robosample/roborun.py 1apq examples/1APQ.prmtop examples/1APQ.rst7 6000 0 100 1
nsys stats report1.nsys-rep # or CLI

nsys profile -o runrex_trace \
  --trace=cuda,osrt,nvtx \
  --sample=cpu \
python3 python/robosample/roborun.py 1apq examples/1APQ.prmtop examples/1APQ.rst7 6000 0 100 1
"""

# python3 python/robosample/roborun.py 1apq examples/1APQ.prmtop examples/1APQ.rst7 6000 0 100 1
# python3 python/robosample/roborun.py ffar1 examples/ffar1.prmtop examples/ffar1.rst7 6000 0 20 1
# python3 python/robosample/roborun.py GfcD examples/GfcDstrippedMin.prmtop examples/GfcDstrippedMin.rst7 6000 0 20 1


# /usr/bin/time -v python3 python/robosample/roborun.py GfcD examples/GfcDstrippedMin.prmtop examples/GfcDstrippedMin.rst7 6000 0 1 1
# nvidia-smi --query-gpu=memory.used --format=csv -lms 100 > vram_usage.csv

# Create the parser
parser = argparse.ArgumentParser(description="Process PDB code and seed.")

# Add the arguments
parser.add_argument("name", type=str, help="Name of the simulation.")
parser.add_argument("prmtop", type=str, help="Relative path to the .prmtop file.")
parser.add_argument("inpcrd", type=str, help="Relative path to the .inpcrd file.")
parser.add_argument("seed", type=int, help="The seed.")
parser.add_argument("equil_steps", type=int, help="The number of equilibration steps.")
parser.add_argument("prod_steps", type=int, help="The number of production steps.")
parser.add_argument("write_freq", type=int, help="CSV and DCD write frequency.")

# Parse the arguments
args = parser.parse_args()

# Temperature replica exchange parameters
T0 = 300.0
T_MAX = 1000.0
NOF_REPLICAS = 1
R = 1 if NOF_REPLICAS == 1 else (T_MAX / T0) ** (1.0 / (NOF_REPLICAS - 1))

# Mean first passage time to cross an energy barrier
# 3 kcal/mol - sub-picosecond to picosecond transitions (modest barrier, ~5KbT)
# 6 kcal/mol - tens to hundreds of picoseconds (moderate barrier, ~10KbT)
# 10 kcal/mol - nanoseconds or longer (high barrier , ~16KbT)
# 2 ps of MD is enough to explore shallow wells, but not to cross deep barriers without enhanced sampling (e.g., HMC, replica exchange)
TIMESTEP_TD = 0.001  # Torsional dymaics time step is 10 fs
MDSTEPS_TD = 64  # Torsional dynamics block trajectory length 1 ps

TIMESTEP_CARTESIAN = 0.001
MDSTEPS_CARTESIAN = 64

# create robosample context
context = robosample.Context(
    name=args.name,
    seed=args.seed,
    prmtop=args.prmtop,
    inpcrd=args.inpcrd,
    write_freq=args.write_freq,
    testing=False,
)

# [[rb.BondFlexibility(), rb.BondFlexibility(), ...], [...], ...]

# Add cartesian world (will integrate with OpenMM)
context.addCartesianWorld().add_sampler(
    timeStep=TIMESTEP_CARTESIAN,
    mdSteps=MDSTEPS_CARTESIAN,
    boostMDSteps=MDSTEPS_CARTESIAN,
    acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
)

# # Add torsional world with standardized dihedrals
# sele = context.build_flexibilities(context.standard_dihedral_bonds)
# context.addTorsionalWorld(sele).add_sampler(
#     timeStep=TIMESTEP_TD,
#     mdSteps=MDSTEPS_TD,
#     boostMDSteps=MDSTEPS_TD,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
# )

num_residues = context.standard_dihedral_bonds["resid"].max() + 1
for resid in range(num_residues):
    bonds = context.standard_dihedral_bonds[
        (context.standard_dihedral_bonds["resid"] == resid)
        & (
            (context.standard_dihedral_bonds["dihedral_type"] == "phi")
            | (context.standard_dihedral_bonds["dihedral_type"] == "psi")
        )
    ]
    if bonds.empty:
        continue

    sele = context.build_flexibilities(bonds)
    context.addTorsionalWorld(sele).add_sampler(
        timeStep=TIMESTEP_TD,
        mdSteps=MDSTEPS_TD,
        boostMDSteps=MDSTEPS_TD,
        acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
    )

# phi = context.standard_dihedral_bonds[
#     context.standard_dihedral_bonds["dihedral_type"] == "phi"
# ]
# sele = context.build_flexibilities(phi)
# context.addTorsionalWorld(sele).add_sampler(
#     timeStep=TIMESTEP_TD,
#     mdSteps=MDSTEPS_TD,
#     boostMDSteps=MDSTEPS_TD,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
# )

# psi = context.standard_dihedral_bonds[
#     context.standard_dihedral_bonds["dihedral_type"] == "psi"
# ]
# sele = context.build_flexibilities(psi)
# context.addTorsionalWorld(sele).add_sampler(
#     timeStep=TIMESTEP_TD,
#     mdSteps=MDSTEPS_TD,
#     boostMDSteps=MDSTEPS_TD,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
# )

sidechain_bonds = context.standard_dihedral_bonds[
    (context.standard_dihedral_bonds["dihedral_type"] != "phi")
    & (context.standard_dihedral_bonds["dihedral_type"] != "psi")
    & (context.standard_dihedral_bonds["dihedral_type"] != "omega")
]
sele = context.build_flexibilities(sidechain_bonds)
context.addTorsionalWorld(sele).add_sampler(
    timeStep=TIMESTEP_TD,
    mdSteps=MDSTEPS_TD,
    boostMDSteps=MDSTEPS_TD,
    acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
)

# for flex in context.getDefaultBonds('macrocycle'):
# 	context.addTorsionalWorld([flex]).add_sampler(timeStep=TIMESTEP_TD, mdSteps=MDSTEPS_TD, boostMDSteps=MDSTEPS_TD, acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept)

# Add replicas (geometric temperature ladder)
temperatures = []
for i in range(NOF_REPLICAS):
    temperatures.append(T0 * (R**i))

context.initialize(temperatures)

# print(context.calculate_openmm_energy(0))
# print(context.calculate_openmm_energy(1))

# # Test OpenMM energies
# robosample_openmm = context.calculate_openmm_energy(1)
# native_openmm = get_native_openmm_energy()
# if robosample_openmm != native_openmm:
# 	message = f"OpenMM energy calculated by Robosample ({robosample_openmm}) does not match native OpenMM energy ({native_openmm})."
# 	raise RuntimeError(message)

# Run the simulation
start_time = time.perf_counter()

context.run_rex(args.equil_steps, args.prod_steps, args.write_freq, True, True)

end_time = time.perf_counter()
duration = end_time - start_time

print(f"run_rex() took {duration:.4f} seconds")

"""
source leaprc.protein.ff19SB
mol1 = sequence { ACE ALA ALA NME }
mol2 = sequence { ACE ALA ALA NME }
translate mol2 { 15.0 0.0 0.0 }
system = combine { mol1 mol2 }
saveamberparm system 2ala.2ala.prmtop 2ala.2ala.inpcrd
saveamberparm system 2ala.2ala.prmtop 2ala.2ala.rst7
savepdb system 2ala.2ala.pdb
quit
"""

"""
source leaprc.protein.ff19SB
system = sequence { ACE GLY ALA VAL LEU ILE SER THR CYS MET PHE TYR TRP PRO ASP GLU ASN GLN HIS LYS ARG NME }
saveamberparm system all.prmtop all.inpcrd
saveamberparm system all.prmtop all.rst7
savepdb system all.pdb
quit
"""
