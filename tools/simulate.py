# import robosample as robosample
import flexor
import mdtraj as md

# 1i42 2btg 2eej 2k9q 2kbt 2kyy 2l39 2m2f 2oa4 2obu
prmtop = "../build/2eej/2eej.H.capped.prmtop"
inpcrd = "../build/2eej/2eej.H.capped.rst7"

# prepare flexor generator
mdtrajObj = md.load(inpcrd, top=prmtop)
flexorObj = flexor.Flexor(mdtrajObj)

flex = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
# flex = flexorObj.create(range="all", distanceCutoff=0, subset=["all"], jointType="Pin", sasa_value=-1.0)
# flex = flexorObj.create(range="all", distanceCutoff=0, subset=["rama"], jointType="Pin", sasa_value=-1.0)
for f in flex:
    print(f)

# # create robosample context
# c = robosample.Context(300, 300, 42)
# c.setNumThreads(0)
# # c.setPdbPrefix("2but42")
# # c.setOutput("temp") # the log file is created like log.[seed]
# c.setNofRoundsTillReblock(10) # per world?
# c.setRequiredNofRounds(2)
# c.setPdbRestartFreq(0) # WRITE_PDBS
# c.setPrintFreq(1) # PRINT_FREQ
# c.setRunType(robosample.RunType.DEFAULT)

# # load system
# c.loadAmberSystem(prmtop, inpcrd)

# c.appendDCDReporter('pro10.dcd')

# # openmm cartesian
# flex = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
# c.addWorld(False, 1, robosample.RootMobility.WELD, flex, True, False, 0)

# # sidechains pins
# flex = flexorObj.create(range="all", distanceCutoff=0, subset=["side"], jointType="Pin", sasa_value=-1.0)
# c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

# # roll
# flex = flexorObj.create(range="all", distanceCutoff=0, subset=["coils", "rama"], jointType="Pin", sasa_value=-1.0)
# for f in flex:
#     c.addWorld(True, 1, robosample.RootMobility.WELD, [f], True, False, 0)

# # world temperature goes to replica

# # samplers
# sampler = robosample.SamplerName.HMC # rename to type
# generator = robosample.AcceptRejectMode.MC # move to thermodynamic state # accept/reject type
# thermostat = robosample.ThermostatName.ANDERSEN

# # add defaults and move to thermodynamic state
#     # x SimTK::Real timestep,
# 	# x int mdStepsPerSample,
# 	# int mdStepsPerSampleStd, - doar pentru noi
# 	# SimTK::Real boostTemperature,
# 	# int boostMDSteps,
# 	# x int distort,
# 	# x int work,
# 	# x int flow,
# 	# bool useFixmanPotential
#     # x MC - EMPTY
#     # x integrator
# c.getWorld(0).addSampler(sampler, generator, robosample.IntegratorName.OMMVV, thermostat, 0.001, 200, 0, 300, 1, 0, 0, 0, False)
# c.getWorld(1).addSampler(sampler, generator, robosample.IntegratorName.VERLET, thermostat, 0.001, 200, 0, 300, 1, 0, 0, 0, True)

# for i in range(len(flex)):
#     c.getWorld(2 + i).addSampler(sampler, generator, robosample.IntegratorName.VERLET, thermostat, 0.001, 200, 0, 300, 1, 0, 0, 0, True)

# # start the simulation
# c.Run()

# numReplicas = 10
# temperatures = []
# for i in range(numReplicas):
#     temperatures.append(300.0 + i * 10.0)

# samplers = [0, 1, 2, 3]
# distortOptions = [0]
# distortArgs = [0]
# flow = [0]
# work = [0]
# integrators = [0]

# worlds = [0, 1, 2, 3]
# timesteps = [0.001, 0.001, 0.001, 0.001]
# mdsteps = [200, 200, 200, 200]

# for i in range(numReplicas):
#     c.addReplica(i)
#     c.addThermodynamicState(i, temperatures[i], temperatures[i], samplers, distortOptions, distortArgs, flow, work, integrators, worlds, timesteps, mdsteps)



# c.initializeZMatrix()


# # addReplica empty
# # add thermodynamic state with all the parameters - indices for sampler, distort options, distort args, flow, work, integrator, worlds, timesteps, mdsteps

# # setReplicasWorldsParameters