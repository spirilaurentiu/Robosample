import robosample as robosample
import flexor
import mdtraj as md

# 1i42 2btg 2eej 2k9q 2kbt 2kyy 2l39 2m2f 2oa4 2obu
prmtop = "pro10/pro10.2.prmtop"
inpcrd = "pro10/pro10.2.rst7"

# prepare flexor generator
mdtrajObj = md.load(inpcrd, top=prmtop)
flexorObj = flexor.Flexor(mdtrajObj)

# create robosample context
c = robosample.Context(300, 300, 42)
c.setNumThreads(0)
# c.setPdbPrefix("2but42")
# c.setOutput("temp") # the log file is created like log.[seed]
c.setNofRoundsTillReblock(10) # per world?
c.setRequiredNofRounds(2)
c.setPdbRestartFreq(0) # WRITE_PDBS
c.setPrintFreq(1) # PRINT_FREQ
c.setRunType(robosample.RunType.DEFAULT)

# load system
c.loadAmberSystem(prmtop, inpcrd)

c.appendDCDReporter('pro10.dcd')

# openmm cartesian
flex = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
c.addWorld(False, 1, robosample.RootMobility.WELD, flex, True, False, 0)

# sidechains pins
flex = flexorObj.create(range="all", distanceCutoff=0, subset=["side"], jointType="Pin", sasa_value=-1.0)
c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

# roll
flex = flexorObj.create(range="all", distanceCutoff=0, subset=["coils", "rama"], jointType="Pin", sasa_value=-1.0)
for f in flex:
    c.addWorld(True, 1, robosample.RootMobility.WELD, [f], True, False, 0)

# samplers
c.getWorld(0).addSampler(robosample.SamplerName.HMC, robosample.SampleGenerator.MC, robosample.IntegratorName.OMMVV, robosample.ThermostatName.ANDERSEN, 0.001, 200, 0, 300, 1, 0, 0, 0, False)
c.getWorld(1).addSampler(robosample.SamplerName.HMC, robosample.SampleGenerator.MC, robosample.IntegratorName.VERLET, robosample.ThermostatName.ANDERSEN, 0.001, 200, 0, 300, 1, 0, 0, 0, True)

for i in range(len(flex)):
    c.getWorld(2 + i).addSampler(robosample.SamplerName.HMC, robosample.SampleGenerator.MC, robosample.IntegratorName.VERLET, robosample.ThermostatName.ANDERSEN, 0.001, 200, 0, 300, 1, 0, 0, 0, True)

# start the simulation
c.Run()
