#include "Context.hpp"

int main() {
    // no need for tini tfin
	Context c("cb8", 42, 0, 1, RUN_TYPE::REMC, 1, 0);

	// c.setNonbonded(0, 1.2); // default - set for each world
	c.setPdbPrefix("cb8"); // will disappear
	c.setOutput("temp"); // the log file is created like log.[seed] - needs refactoring
	c.setRequiredNofRounds(10); // per world? what does it do?
	c.setPdbRestartFreq(1); // WRITE_PDBS
	c.setPrintFreq(1); // PRINT_FREQ

	c.loadAmberSystem("cb8/ligand.prmtop", "cb8/ligand.rst7");

	// World 0 OPENMM
	std::vector<BOND_FLEXIBILITY> flexibilities_w0 = {
		{0, 12, BondMobility::Mobility::Translation},
		{0, 13, BondMobility::Mobility::Translation},
		{0, 3, BondMobility::Mobility::Translation},
		{0, 1, BondMobility::Mobility::Translation},
		{1, 15, BondMobility::Mobility::Translation},
		{1, 14, BondMobility::Mobility::Translation},
		{1, 2, BondMobility::Mobility::Translation},
		{2, 5, BondMobility::Mobility::Translation},
		{2, 10, BondMobility::Mobility::Translation},
		{2, 7, BondMobility::Mobility::Translation},
		{3, 4, BondMobility::Mobility::Translation},
		{3, 6, BondMobility::Mobility::Translation},
		{3, 8, BondMobility::Mobility::Translation},
		{4, 16, BondMobility::Mobility::Translation},
		{4, 17, BondMobility::Mobility::Translation},
		{4, 5, BondMobility::Mobility::Translation},
		{5, 19, BondMobility::Mobility::Translation},
		{5, 18, BondMobility::Mobility::Translation},
		{6, 21, BondMobility::Mobility::Translation},
		{6, 20, BondMobility::Mobility::Translation},
		{6, 7, BondMobility::Mobility::Translation},
		{7, 22, BondMobility::Mobility::Translation},
		{7, 23, BondMobility::Mobility::Translation},
		{8, 25, BondMobility::Mobility::Translation},
		{8, 24, BondMobility::Mobility::Translation},
		{8, 9, BondMobility::Mobility::Translation},
		{9, 26, BondMobility::Mobility::Translation},
		{10, 27, BondMobility::Mobility::Translation},
		{10, 28, BondMobility::Mobility::Translation},
		{10, 11, BondMobility::Mobility::Translation},
		{11, 29, BondMobility::Mobility::Translation},
		{30, 72, BondMobility::Mobility::Translation},
		{31, 76, BondMobility::Mobility::Translation},
		{32, 80, BondMobility::Mobility::Translation},
		{33, 84, BondMobility::Mobility::Translation},
		{34, 88, BondMobility::Mobility::Translation},
		{35, 92, BondMobility::Mobility::Translation},
		{36, 96, BondMobility::Mobility::Translation},
		{37, 128, BondMobility::Mobility::Translation},
		{38, 132, BondMobility::Mobility::Translation},
		{39, 136, BondMobility::Mobility::Translation},
		{40, 140, BondMobility::Mobility::Translation},
		{41, 144, BondMobility::Mobility::Translation},
		{42, 148, BondMobility::Mobility::Translation},
		{43, 152, BondMobility::Mobility::Translation},
		{44, 72, BondMobility::Mobility::Translation},
		{44, 97, BondMobility::Mobility::Translation},
		{44, 100, BondMobility::Mobility::Translation},
		{45, 72, BondMobility::Mobility::Translation},
		{45, 73, BondMobility::Mobility::Translation},
		{45, 102, BondMobility::Mobility::Translation},
		{46, 73, BondMobility::Mobility::Translation},
		{46, 76, BondMobility::Mobility::Translation},
		{46, 104, BondMobility::Mobility::Translation},
		{47, 76, BondMobility::Mobility::Translation},
		{47, 77, BondMobility::Mobility::Translation},
		{47, 106, BondMobility::Mobility::Translation},
		{48, 77, BondMobility::Mobility::Translation},
		{48, 80, BondMobility::Mobility::Translation},
		{48, 108, BondMobility::Mobility::Translation},
		{49, 80, BondMobility::Mobility::Translation},
		{49, 81, BondMobility::Mobility::Translation},
		{49, 110, BondMobility::Mobility::Translation},
		{50, 81, BondMobility::Mobility::Translation},
		{50, 84, BondMobility::Mobility::Translation},
		{50, 112, BondMobility::Mobility::Translation},
		{51, 84, BondMobility::Mobility::Translation},
		{51, 85, BondMobility::Mobility::Translation},
		{51, 114, BondMobility::Mobility::Translation},
		{52, 85, BondMobility::Mobility::Translation},
		{52, 88, BondMobility::Mobility::Translation},
		{52, 116, BondMobility::Mobility::Translation},
		{53, 88, BondMobility::Mobility::Translation},
		{53, 89, BondMobility::Mobility::Translation},
		{53, 118, BondMobility::Mobility::Translation},
		{54, 89, BondMobility::Mobility::Translation},
		{54, 92, BondMobility::Mobility::Translation},
		{54, 120, BondMobility::Mobility::Translation},
		{55, 92, BondMobility::Mobility::Translation},
		{55, 93, BondMobility::Mobility::Translation},
		{55, 122, BondMobility::Mobility::Translation},
		{56, 93, BondMobility::Mobility::Translation},
		{56, 96, BondMobility::Mobility::Translation},
		{56, 124, BondMobility::Mobility::Translation},
		{57, 96, BondMobility::Mobility::Translation},
		{57, 97, BondMobility::Mobility::Translation},
		{57, 126, BondMobility::Mobility::Translation},
		{58, 100, BondMobility::Mobility::Translation},
		{58, 128, BondMobility::Mobility::Translation},
		{58, 153, BondMobility::Mobility::Translation},
		{59, 102, BondMobility::Mobility::Translation},
		{59, 128, BondMobility::Mobility::Translation},
		{59, 129, BondMobility::Mobility::Translation},
		{60, 104, BondMobility::Mobility::Translation},
		{60, 129, BondMobility::Mobility::Translation},
		{60, 132, BondMobility::Mobility::Translation},
		{61, 106, BondMobility::Mobility::Translation},
		{61, 132, BondMobility::Mobility::Translation},
		{61, 133, BondMobility::Mobility::Translation},
		{62, 108, BondMobility::Mobility::Translation},
		{62, 133, BondMobility::Mobility::Translation},
		{62, 136, BondMobility::Mobility::Translation},
		{63, 110, BondMobility::Mobility::Translation},
		{63, 136, BondMobility::Mobility::Translation},
		{63, 137, BondMobility::Mobility::Translation},
		{64, 112, BondMobility::Mobility::Translation},
		{64, 137, BondMobility::Mobility::Translation},
		{64, 140, BondMobility::Mobility::Translation},
		{65, 114, BondMobility::Mobility::Translation},
		{65, 140, BondMobility::Mobility::Translation},
		{65, 141, BondMobility::Mobility::Translation},
		{66, 116, BondMobility::Mobility::Translation},
		{66, 141, BondMobility::Mobility::Translation},
		{66, 144, BondMobility::Mobility::Translation},
		{67, 118, BondMobility::Mobility::Translation},
		{67, 144, BondMobility::Mobility::Translation},
		{67, 145, BondMobility::Mobility::Translation},
		{68, 120, BondMobility::Mobility::Translation},
		{68, 145, BondMobility::Mobility::Translation},
		{68, 148, BondMobility::Mobility::Translation},
		{69, 122, BondMobility::Mobility::Translation},
		{69, 148, BondMobility::Mobility::Translation},
		{69, 149, BondMobility::Mobility::Translation},
		{70, 124, BondMobility::Mobility::Translation},
		{70, 149, BondMobility::Mobility::Translation},
		{70, 152, BondMobility::Mobility::Translation},
		{71, 126, BondMobility::Mobility::Translation},
		{71, 152, BondMobility::Mobility::Translation},
		{71, 153, BondMobility::Mobility::Translation},
		{73, 74, BondMobility::Mobility::Translation},
		{73, 75, BondMobility::Mobility::Translation},
		{77, 78, BondMobility::Mobility::Translation},
		{77, 79, BondMobility::Mobility::Translation},
		{81, 82, BondMobility::Mobility::Translation},
		{81, 83, BondMobility::Mobility::Translation},
		{85, 86, BondMobility::Mobility::Translation},
		{85, 87, BondMobility::Mobility::Translation},
		{89, 90, BondMobility::Mobility::Translation},
		{89, 91, BondMobility::Mobility::Translation},
		{93, 94, BondMobility::Mobility::Translation},
		{93, 95, BondMobility::Mobility::Translation},
		{97, 98, BondMobility::Mobility::Translation},
		{97, 99, BondMobility::Mobility::Translation},
		{100, 101, BondMobility::Mobility::Translation},
		{100, 102, BondMobility::Mobility::Translation},
		{102, 103, BondMobility::Mobility::Translation},
		{104, 105, BondMobility::Mobility::Translation},
		{104, 106, BondMobility::Mobility::Translation},
		{106, 107, BondMobility::Mobility::Translation},
		{108, 109, BondMobility::Mobility::Translation},
		{108, 110, BondMobility::Mobility::Translation},
		{110, 111, BondMobility::Mobility::Translation},
		{112, 113, BondMobility::Mobility::Translation},
		{112, 114, BondMobility::Mobility::Translation},
		{114, 115, BondMobility::Mobility::Translation},
		{116, 117, BondMobility::Mobility::Translation},
		{116, 118, BondMobility::Mobility::Translation},
		{118, 119, BondMobility::Mobility::Translation},
		{120, 121, BondMobility::Mobility::Translation},
		{120, 122, BondMobility::Mobility::Translation},
		{122, 123, BondMobility::Mobility::Translation},
		{124, 125, BondMobility::Mobility::Translation},
		{124, 126, BondMobility::Mobility::Translation},
		{126, 127, BondMobility::Mobility::Translation},
		{129, 130, BondMobility::Mobility::Translation},
		{129, 131, BondMobility::Mobility::Translation},
		{133, 134, BondMobility::Mobility::Translation},
		{133, 135, BondMobility::Mobility::Translation},
		{137, 138, BondMobility::Mobility::Translation},
		{137, 139, BondMobility::Mobility::Translation},
		{141, 142, BondMobility::Mobility::Translation},
		{141, 143, BondMobility::Mobility::Translation},
		{145, 146, BondMobility::Mobility::Translation},
		{145, 147, BondMobility::Mobility::Translation},
		{149, 150, BondMobility::Mobility::Translation},
		{149, 151, BondMobility::Mobility::Translation},
		{153, 154, BondMobility::Mobility::Translation},
		{153, 155, BondMobility::Mobility::Translation}
	};
	c.addWorld(false, 1, ROOT_MOBILITY::WELD, flexibilities_w0);

	// World 1
	std::vector<BOND_FLEXIBILITY> flexibilities_w1 = {
		{2, 10, BondMobility::Mobility::Torsion},
		{3, 8, BondMobility::Mobility::Torsion},
		{8, 9, BondMobility::Mobility::Torsion},
		{10, 11, BondMobility::Mobility::Torsion},
	};
	c.addWorld(true, 1, ROOT_MOBILITY::WELD, flexibilities_w1);

	// // World 2
	// std::vector<BOND_FLEXIBILITY> flexibilities_w2 = {
	// };
	// c.addWorld(true, 1, ROOT_MOBILITY::WELD, flexibilities_w2);

	// Does OMMVV have to be first?
	// Add samplers
	c.getWorld(0).addSampler(SamplerName::HMC, IntegratorName::OMMVV, ThermostatName::ANDERSEN, false);

	// // Does not work if I set OMMVV instead of VV. How do I check if it is working?
	c.getWorld(1).addSampler(SamplerName::HMC, IntegratorName::Verlet, ThermostatName::ANDERSEN, true);

	// // // Does not work if I set OMMVV instead of VV. How do I check if it is working?
	// c.getWorld(2).addSampler(SamplerName::HMC, IntegratorName::Verlet, ThermostatName::ANDERSEN, true);

	int nofReplicas = 1;
	SimTK::Real temperature = 300;
	std::vector<SimTK::Real> temperatures, boostTemperatures;
	for (int i = 0; i < nofReplicas; i++) {
		temperatures.push_back(temperature + (i * 100));
		boostTemperatures.push_back(temperature + (i * 100)); // used for openmm velocities
	}

	std::vector<AcceptRejectMode> acceptRejectModes = { AcceptRejectMode::MetropolisHastings, AcceptRejectMode::MetropolisHastings, AcceptRejectMode::MetropolisHastings };
	std::vector<SimTK::Real> timesteps = { 0.001, 0.006, 0.008 };
	std::vector<int> worldIndexes = { 0, 1, 2 };
	std::vector<int> mdsteps = { 50, 10, 20 };
	std::vector<int> boostMDSteps = { 50, 10, 20 };
	std::vector<int> samplesPerRound = { 1, 1, 1 };

	std::vector<int> distortOptions = { 0, 0, 0 };
	std::vector<std::string> distortArgs = { "0", "0", "0" };
	std::vector<int> flow = { 0, 0, 0 };
	std::vector<int> work = { 0, 0, 0 };

	acceptRejectModes.pop_back();
	timesteps.pop_back();
	worldIndexes.pop_back();
	mdsteps.pop_back();
	boostMDSteps.pop_back();
	samplesPerRound.pop_back();
	distortOptions.pop_back();
	distortArgs.pop_back();
	flow.pop_back();
	work.pop_back();

	for (int i = 0; i < nofReplicas; i++) {
		c.addReplica(i);
		c.addThermodynamicState(i,
			temperatures[i],
			boostTemperatures[i],
			acceptRejectModes,
			distortOptions,
			distortArgs,
			flow,
			work,
			worldIndexes,
			timesteps,
			mdsteps,
			boostMDSteps
			);
	}

	c.initialize();

	// pas how many rounds to run for here
	c.RunREXNew(0, 1);

	return 0;
}