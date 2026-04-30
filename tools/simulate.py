import argparse

import batstat
import flexor
import mdtraj as md
import robosample

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
parser.add_argument(
    "temperature_init", type=int, help="Temperature of the first replica."
)
parser.add_argument("pdbid", type=str, help="pdbid")
parser.add_argument("type", type=str, help="type of simulation")
parser.add_argument("spr_method", type=str, help="samples per round method")

# Parse the arguments
args = parser.parse_args()

# prepare flexor generator
mdtrajObj = md.load(args.inpcrd, top=args.prmtop)
flexorObj = flexor.Flexor(mdtrajObj)

# create robosample context
context = robosample.Context(args.name, args.seed, 0, 1, robosample.RunType.REMC, 1, 0)
context.setPdbRestartFreq(0)  # WRITE_PDBS
context.setPrintFreq(args.write_freq)  # PRINT_FREQ
context.setNonbonded(0, 1.2)
context.setGBSA(1)
context.setVerbose(False)

# load system
context.loadAmberSystem(args.prmtop, args.inpcrd)

if args.type == "tdnr":
    # Do torsional dynamics - non-redundant dihedrals
    stats = batstat.BATCorrelations(args.prmtop, args.inpcrd)
    atom_indices = stats.get_dihedral_atom_indices()

    flex = flexorObj.create_from_list(atom_indices, robosample.BondMobility.Torsion)
    context.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)
    context.getWorld(0).add_sampler(
        robosample.SamplerName.HMC,
        robosample.IntegratorType.OMMVV,
        robosample.ThermostatName.ANDERSEN,
        False,
    )

    nof_replicas = 1
    temperature = args.temperature_init
    temperatures = []
    boost_temperatures = []
    for i in range(nof_replicas):
        temperatures.append(temperature + (i * 10))
        boost_temperatures.append(temperature + (i * 10))  # used for openmm velocities

    accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings]
    timesteps = [0.0075]
    worldIndexes = [0]
    world_indexes = [0]
    mdsteps = [10]
    boost_md_steps = mdsteps
    integrators = [robosample.IntegratorType.VERLET]

    distort_options = [0]
    distort_args = ["0"]
    flow = [0]
    work = [0]

    for i in range(nof_replicas):
        context.addReplica(i)
        context.addThermodynamicState(
            i,
            temperatures[i],
            accept_reject_modes,
            distort_options,
            distort_args,
            flow,
            work,
            integrators,
            worldIndexes,
            timesteps,
            mdsteps,
        )

    context.Initialize()
    context.RunREX(args.equil_steps, args.prod_steps)

elif args.type == "tdc":
    # Do torsional dynamics - correlated dihedrals

    # # Cluster the previous simulations
    stats = batstat.BATCorrelations(args.prmtop, args.inpcrd)
    atom_indices = stats.get_dihedral_atom_indices()
    # correlation_file = f"{args.pdbid}_6000_correlation.npy"

    # if os.path.exists(correlation_file):
    # 	corr = np.load(correlation_file)
    # else:
    # 	dcd_files = [f"{args.pdbid}_{i}.dcd" for i in range(5)]
    # 	stats.compute_dihedrals_from_dcd(dcd_files)
    # 	corr = stats.compute_correlations()
    # 	np.save(correlation_file, corr)

    # # Partition the dihedrals into blocks
    # blocks, samples_per_round = stats.dynamic_partitioning(np.mean(np.abs(corr), axis=0))
    # if args.spr_method == 'fixed':
    # 	samples_per_round = [1] * len(blocks)
    # # print("samples_per_round", samples_per_round)

    # python3 simulate.py 1APQ_new data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6000 10 100 1 300 1APQ tdc auto
    blocks = [
        # high correlation blocks
        [0, 2, 3, 4, 8, 10, 11, 12, 13, 53, 55, 57, 58, 64, 65, 66, 89],
        [
            73,
            74,
            75,
            77,
            78,
            79,
            80,
            81,
            82,
            83,
            85,
            86,
            87,
            91,
            92,
            93,
            94,
            96,
            97,
            98,
            99,
            110,
        ],
        [14, 15, 16, 18, 20, 21, 25, 28, 30, 31, 33, 34, 38, 40, 41, 42, 103, 107],
        [49, 50, 51, 52, 67, 68, 71, 108, 109],
        # low correlation block
        [
            1,
            5,
            6,
            9,
            54,
            56,
            59,
            60,
            61,
            62,
            63,
            76,
            84,
            88,
            90,
            95,
            100,
            101,
            102,
            104,
            111,
            7,
            17,
            19,
            22,
            23,
            24,
            26,
            27,
            29,
            32,
            35,
            36,
            37,
            39,
            43,
            44,
            45,
            46,
            47,
            48,
            69,
            70,
            72,
            106,
        ],
    ]

    # # flatten the blocks

    # # add another block that contins all the non-correlated dihedrals
    # uncorrelated_block = []
    # for i in range(stats.num_dihe):
    # 	if i not in [item for sublist in blocks for item in sublist]:
    # 		uncorrelated_block.append(i)

    blocks_as_bond_list = []
    for block in blocks:
        l = []
        for b in block:
            aix1 = stats.atom_indices[b][1]
            aix2 = stats.atom_indices[b][2]
            l.append([aix1, aix2])
        blocks_as_bond_list.append(l)
    blocks = blocks_as_bond_list

    samples_per_round = [
        0.3291,
        0.1206,
        0.3612,
        0.1587,
        0.1206 * 0.5,  # adjust, this is for the non-correlated pairs
    ]

    # nohup python3 simulate.py 1APQ_new data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6000 1 100000 1 300 1APQ tdc auto > /dev/null 2>&1 &
    # nohup python3 simulate.py 1APQ_new data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6001 1 100000 1 300 1APQ tdc auto > /dev/null 2>&1 &
    # nohup python3 simulate.py 1APQ_new data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6002 1 100000 1 300 1APQ tdc auto > /dev/null 2>&1 &
    # nohup python3 simulate.py 1APQ_new data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6003 1 100000 1 300 1APQ tdc auto > /dev/null 2>&1 &
    # nohup python3 simulate.py 1APQ_new data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6004 1 100000 1 300 1APQ tdc auto > /dev/null 2>&1 &
    # nohup python3 simulate.py 1APQ_new data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6005 1 100000 1 300 1APQ tdc auto > /dev/null 2>&1 &

    # normalize such that min(samples_per_round) = 1
    min_samples = min(samples_per_round)
    samples_per_round = [int(round(x / min_samples)) for x in samples_per_round]

    # Add Cartesian flexors (OpenMM)
    flex = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
    context.addWorld(False, 1, robosample.RootMobility.WELD, flex, True, False, 0)

    # Create the flexors from the blocks
    for block, num_samples in zip(blocks, samples_per_round):
        flex = flexorObj.create_from_list(block, robosample.BondMobility.Torsion)
        context.addWorld(
            True, num_samples, robosample.RootMobility.WELD, flex, True, False, 0
        )

    non_correlated_world_index = len(blocks) + 1
    context.getWorld(1).setRollFlexibilities(True)

    # Add samplers
    context.getWorld(0).add_sampler(
        robosample.SamplerName.HMC,
        robosample.IntegratorType.OMMVV,
        robosample.ThermostatName.ANDERSEN,
        False,
    )
    for i in range(1, len(blocks) + 1):
        context.getWorld(i).add_sampler(
            robosample.SamplerName.HMC,
            robosample.IntegratorType.VERLET,
            robosample.ThermostatName.ANDERSEN,
            True,
        )

    nof_replicas = 1
    temperature = args.temperature_init
    temperatures = []
    boost_temperatures = []
    for i in range(nof_replicas):
        temperatures.append(temperature + (i * 10))
        boost_temperatures.append(temperature + (i * 10))  # used for openmm velocities

    accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings] * (
        len(blocks) + 1
    )
    timesteps = [0.0007] + [0.0075] * len(blocks)
    worldIndexes = range(len(blocks) + 1)
    world_indexes = range(len(blocks) + 1)
    mdsteps = [10000] + [150] * len(blocks)
    boost_md_steps = mdsteps
    integrators = [robosample.IntegratorType.OMMVV] + [
        robosample.IntegratorType.VERLET
    ] * len(blocks)

    distort_options = [0] * (len(blocks) + 1)
    distort_args = ["0"] * (len(blocks) + 1)
    flow = [0] * (len(blocks) + 1)
    work = [0] * (len(blocks) + 1)

    for i in range(nof_replicas):
        context.addReplica(i)
        context.addThermodynamicState(
            i,
            temperatures[i],
            accept_reject_modes,
            distort_options,
            distort_args,
            flow,
            work,
            integrators,
            worldIndexes,
            timesteps,
            mdsteps,
        )

    context.Initialize()
    context.RunREX(args.equil_steps, args.prod_steps)

else:
    # # sidechains pins
    # flex = flexorObj.create(range="all", distanceCutoff=0, subset=["all"], jointType="Pin", sasa_value=-1.0)
    # c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

    # # ramachandran pins
    # flex = flexorObj.create(range="all", distanceCutoff=0, subset=["rama"], jointType="Pin", sasa_value=-1.0)
    # c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

    # c.getWorld(1).add_sampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)
    # c.getWorld(2).add_sampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)

    pass


# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6000 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6001 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6002 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6003 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6004 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6005 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6006 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6007 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6008 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6009 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6010 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6011 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6012 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6013 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6014 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
# nohup python3 simulate.py 1APQ_tdnr data-raw/1APQ.prmtop data-raw/1APQ_min.inpcrd 6015 100 1000000 100 300 1APQ tdnr auto > /dev/null 2>&1 &
