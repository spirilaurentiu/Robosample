import robopy_d as robosample
import flexor
import mdtraj as md

prmtop = "2but/ligand.prmtop"
inpcrd = "2but/ligand.rst7"

# prepare flexor generator
mdtrajObj = md.load(inpcrd, top=prmtop)
flexorObj = flexor.Flexor(mdtrajObj)


c = robosample.Context(300, 300, 42)
c.setNumThreads(0)
# c.setPdbPrefix("2but42")
# c.setOutput("temp") # the log file is created like log.[seed]
c.setNofRoundsTillReblock(10) # per world?
c.setRequiredNofRounds(1)
c.setPdbRestartFreq(0) # WRITE_PDBS
c.setPrintFreq(1) # PRINT_FREQ
c.setRunType(robosample.RUN_TYPE.DEFAULT)

c.loadAmberSystem("2but/ligand.prmtop", "2but/ligand.rst7")

# flex_w0 = [
#     robosample.BondFlexibility(3, 1, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(1, 2, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(2, 4, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(1, 0, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(3, 8, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(3, 9, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(3, 10, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(0, 14, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(1, 5, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(2, 6, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(2, 7, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(4, 11, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(4, 12, robosample.BondMobility.Translation),
# 	robosample.BondFlexibility(4, 13, robosample.BondMobility.Translation)
# ]
flex_w0 = flexorObj.create(range="all", distanceCutoff=0, subset=["all"], jointType="Cartesian", sasa_value=-1.0)
c.addWorld(False, 1, robosample.RootMobility.WELD, flex_w0, True, False, 0)

# flex_w1 = [
# 	robosample.BondFlexibility(3, 1, robosample.BondMobility.Torsion),
# 	robosample.BondFlexibility(1, 2, robosample.BondMobility.Torsion),
# 	robosample.BondFlexibility(2, 4, robosample.BondMobility.Torsion),
# 	robosample.BondFlexibility(1, 0, robosample.BondMobility.Torsion)
# ]
flex_w1 = flexorObj.create(range="all", distanceCutoff=0, subset=["all"], jointType="Cartesian", sasa_value=-1.0)
c.addWorld(True, 1, robosample.RootMobility.WELD, flex_w1, True, False, 0)

c.getWorld(0).addSampler(robosample.SamplerName.HMC, robosample.SampleGenerator.MC, robosample.IntegratorName.OMMVV, robosample.ThermostatName.ANDERSEN, 0.001, 200, 0, 300, 1, 0, 0, 0, False)
c.getWorld(1).addSampler(robosample.SamplerName.HMC, robosample.SampleGenerator.MC, robosample.IntegratorName.VERLET, robosample.ThermostatName.ANDERSEN, 0.001, 50, 0, 300, 1, 0, 0, 0, True)

c.Run()
