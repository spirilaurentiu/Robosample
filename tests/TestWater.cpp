#include "Robo.hpp"
#include "Context.hpp"
#include <fstream>
#include <array>

constexpr auto PERCENT = 0.25;
constexpr auto BONDLEN_EQ = 9.73000000E-01; // water bond lenght equilibrium in angstroms
constexpr auto BONDLEN_BEGIN = BONDLEN_EQ - (BONDLEN_EQ * PERCENT);
constexpr auto BONDLEN_END = BONDLEN_EQ + (BONDLEN_EQ * PERCENT);
constexpr auto BOND_FORCE_CONSTANT = 5.63500000E+02; // kcal / (mol * angstroms^2)

constexpr auto ANGLE_EQ = 1.85860192;
constexpr auto ANGLE_BEGIN = ANGLE_EQ - (ANGLE_EQ * PERCENT);
constexpr auto ANGLE_END = ANGLE_EQ + (ANGLE_EQ * PERCENT);
constexpr auto ANGLE_FORCE_CONSTANT = 4.22000000E+01;

constexpr auto TEMPERATURE = 300.0;    
constexpr auto NPOINTS = 1'000;
constexpr auto DECIMAL_PLACES = 2;

constexpr auto KB = 0.0019872041; // Boltzmann's constant in kcal / (mol * kelvin)

bool areEqual(SimTK::Real a, SimTK::Real b) {
    SimTK::Real epsilon = std::pow(10, -DECIMAL_PLACES);
    return std::abs(a - b) < epsilon;
}

void CreateXRange(SimTK::Real start, SimTK::Real stop, std::size_t size, std::vector<SimTK::Real>& range, SimTK::Real& dx) {
    range.reserve(size);

    dx = (stop - start) / static_cast<SimTK::Real>(size - 1);

    for (std::size_t i = 0; i < size; ++i) {
        range.push_back(start + dx * static_cast<SimTK::Real>(i));
    }
}

std::vector<SimTK::Real> CreateMaxwellBoltzmannDistribution(SimTK::Real k, SimTK::Real x0, SimTK::Real beta, SimTK::Real dx, const std::vector<SimTK::Real>& range) {
    std::vector<SimTK::Real> p;
    p.reserve(range.size());

    SimTK::Real Z = sqrt(SimTK::Pi / (beta * k));

    for (std::size_t i = 0; i < range.size(); i++) {
        p.push_back(exp(-beta * k * (range[i] - x0) * (range[i] - x0)) * dx / Z);
    }

    return p;
}

// can be used only for sample szie <= 25, so there went 12 hours of my life
SimTK::Real AndersonDarling(const std::vector<SimTK::Real>& sample, const std::vector<SimTK::Real>& reference) {
    auto n = static_cast<int64_t>(sample.size());
    SimTK::Real sampleCDF = 0.0,
        referenceCDF = 0.0,
        A2 = 0.0,
        epsilon = 1e-16;

    for(int64_t i = 0; i < n; i++) {
        sampleCDF += sample[i];
        referenceCDF += reference[i];

        auto term_0 = log(sampleCDF);
        if (std::isnan(term_0) || std::isinf(term_0)) term_0 = epsilon;

        auto term_1 = log(1.0 - sampleCDF);
        if (std::isnan(term_1) || std::isinf(term_1)) term_1 = epsilon;

        auto term_2 = log(referenceCDF);
        if (std::isnan(term_2) || std::isinf(term_2)) term_2 = epsilon;

        auto term_3 = log(1.0 - referenceCDF);
        if (std::isnan(term_3) || std::isinf(term_3)) term_3 = epsilon;

        A2 -= (2 * (i - 1) - 1) * (term_0 + term_1 - term_2 - term_3);
        std::cout << "i=" << i << ", A2=" << A2 << ", sampleCDF=" << sampleCDF << ", referenceCDF" << referenceCDF << std::endl;
        std::cout << "\t";
        std::cout << "log(sampleCDF)=" << term_0 << ", log(1.0 - sampleCDF)=" << term_1 << 
            ", log(referenceCDF)=" << term_2 << ", log(1.0 - referenceCDF)=" << term_3 << std::endl;
    }

    A2 /= n;
    A2 += 1.0 / n;
    return A2;
}

bool TestMaxwellBoltzmannDistribution(
    const std::vector<SimTK::Real>& sample,
    SimTK::Real k,
    SimTK::Real mu,
    SimTK::Real begin,
    SimTK::Real end,
    SimTK::Real npoints,
    SimTK::Real temperature,
    const std::string& name)
{
    // set constants
    SimTK::Real beta = 1 / (KB * temperature); // 1 / kcal
    
    // set points where we want to compute our target distribution
    std::vector<SimTK::Real> xrange;
    SimTK::Real dx;
    CreateXRange(begin, end, npoints, xrange, dx);

    // this is the theoretical maxwell boltzmann distribution and we will test against it
    const auto theoretical_pdf = CreateMaxwellBoltzmannDistribution(k, mu, beta, dx, xrange);

    // calculate theoretical distribution statitics
    // probably there is an analytical way to do this, but i have not been bestowed upon with this knowledge
    SimTK::Real th_expectance = 0.0,
        th_expectance_sq = 0.0,
        th_variance = 0.0,
        th_stdev = 0.0;
    for(std::size_t i = 0; i < npoints; i++) {
        th_expectance += xrange[i] * theoretical_pdf[i];
        th_expectance_sq += xrange[i] * xrange[i] * theoretical_pdf[i];
    }

    th_variance = th_expectance_sq - th_expectance * th_expectance;
    th_stdev = sqrt(th_variance);

    // now compare values to actual distribution
    SimTK::Real ac_mean = 0.0;
    for (std::size_t i = 0; i < sample.size(); i++) {
        ac_mean += sample[i];
    }
    ac_mean /= static_cast<SimTK::Real>(sample.size());

    SimTK::Real ac_variance = 0.0,
        ac_stdev = 0.0;
    for (std::size_t i = 0; i < sample.size(); i++) {
        ac_variance += (sample[i] - ac_mean) * (sample[i] - ac_mean);
    }
    ac_variance /= static_cast<SimTK::Real>(sample.size());
    ac_stdev = sqrt(ac_variance);

    std::cout << "th_expectance=" << th_expectance << ", th_stdev=" << th_stdev << std::endl;
    std::cout << "ac_mean=" << ac_mean << ", ac_stdev=" << ac_stdev << std::endl;

    if (!areEqual(th_expectance, ac_mean)) {
        std::cout << "Failed Boltzmann distribution test for " << name << std::endl;
        std::cout << "   expected mean = " << th_expectance << ", actual mean = " << ac_mean << std::endl;
        return false;
    }

    if (!areEqual(th_stdev, ac_stdev)) {
        std::cout << "Failed Boltzmann distribution test for " << name << std::endl;
        std::cout << "   expected stdev = " << th_stdev << ", actual stdev = " << ac_stdev << std::endl;
        return false;
    }

    return true;
}

bool TestSlider01() {
    Context c;
    if (!c.initializeFromFile("inp.wat.slider01")) {
        std::cout << "Failed to initialize from inp.wat.slider01" << std::endl;
        return false;
    }

	const auto Rounds = c.getRequiredNofRounds();
    std::vector<SimTK::Real> bl01, bl02; // bond lengths between atom 0-1 and 0-2
    std::vector<SimTK::Real> angles;
    std::vector<SimTK::Real> energies;

	bl01.reserve(Rounds);
	bl02.reserve(Rounds);
	angles.reserve(Rounds);
	energies.reserve(Rounds);

    for (int i = 0; i < Rounds; i++) {
        c.RunOneRound();

        const auto b01 = c.Distance(0, 0, 0, 0, 1);
        bl01.push_back(b01 * 10);

		const auto b02 = c.Distance(0, 0, 0, 0, 2);
		bl02.push_back(b02 * 10);

        const auto a = c.Roboangle(0, 0, 0, 0, 1, 2);
        angles.push_back(a);

        const auto e = c.getPotentialEnergy(0, 0);
        energies.push_back(e);
    }

    // test bond lengths for the the mobile atom pair
    if (!TestMaxwellBoltzmannDistribution(bl01, BOND_FORCE_CONSTANT, BONDLEN_EQ, BONDLEN_BEGIN, BONDLEN_END, NPOINTS, TEMPERATURE, "slider01")) {
        return false;
    }

    // test if the other atoms did not move relative to eachother
    for (auto b : bl02) {
        if (!areEqual(b, bl02[0])) {
            std::cout << "Atom distance changed for inp.wat.slider01" << std::endl;
            std::cout << "   expected = " << bl02[0] << ", actual = " << b << std::endl;
            return false;
        }
    }

    for (auto a : angles) {
        if (!areEqual(a, angles[0])) {
            std::cout << "Atom angle changed for inp.wat.slider01" << std::endl;
            std::cout << "   expected = " << angles[0] << ", actual = " << a << std::endl;
            return false;
        }
    }

    return true;
}

bool TestSlider01BS02() {
    Context c;
    if (!c.initializeFromFile("inp.wat.slider01.bs02")) {
        std::cout << "Failed to initialize from inp.wat.slider01.bs02" << std::endl;
        return false;
    }

	const auto Rounds = c.getRequiredNofRounds();
    std::vector<SimTK::Real> bl01, bl02; // bond lengths between atom 0-1 and 0-2
    std::vector<SimTK::Real> angles;
    std::vector<SimTK::Real> energies;

	bl01.reserve(Rounds);
	bl02.reserve(Rounds);
	angles.reserve(Rounds);
	energies.reserve(Rounds);

    for (int i = 0; i < Rounds; i++) {
        c.RunOneRound();

        const auto b01 = c.Distance(0, 0, 0, 0, 1);
        bl01.push_back(b01 * 10);

		const auto b02 = c.Distance(0, 0, 0, 0, 2);
		bl02.push_back(b02 * 10);

        const auto a = c.Roboangle(0, 0, 0, 0, 1, 2);
        angles.push_back(a);

        const auto e = c.getPotentialEnergy(0, 0);
        energies.push_back(e);

        const SimTK::State& pdbState = c.getWorld(0)->integ->updAdvancedState();
        c.getWorld(0)->updateAtomListsFromCompound(pdbState);
        for(int mol_i = 0; mol_i < c.getNofMolecules(); mol_i++) {
            c.getWorld(0)->getTopology(mol_i).writeAtomListPdb(c.getOutputDir(),
                "/pdbs/sb." +
                c.getPdbPrefix() + "." + std::to_string(mol_i) + ".",
                ".pdb", 10, i);
        }
    }

    std::ofstream fdistances("distances01.txt");
    for (const auto d : bl01) {
        fdistances << d << std::endl;
    }

	std::ofstream fdistances2("distances02.txt");
	for (const auto d : bl02) {
		fdistances2 << d << std::endl;
	}

    std::ofstream fangles("angles.txt");
    for (const auto a : angles) {
        fangles << a << std::endl;
    }

    // test bond lengths for the the mobile atom pair
    if (!TestMaxwellBoltzmannDistribution(bl01, BOND_FORCE_CONSTANT, BONDLEN_EQ, BONDLEN_BEGIN, BONDLEN_END, NPOINTS, TEMPERATURE, "slider01.bs02")) {
        return false;
    }

    if (!TestMaxwellBoltzmannDistribution(bl02, BOND_FORCE_CONSTANT, BONDLEN_EQ, BONDLEN_BEGIN, BONDLEN_END, NPOINTS, TEMPERATURE, "slider01.bs02")) {
        return false;
    }

    if (!TestMaxwellBoltzmannDistribution(angles, ANGLE_FORCE_CONSTANT, ANGLE_EQ, ANGLE_BEGIN, ANGLE_END, NPOINTS, TEMPERATURE, "slider01.bs02")) {
        return false;
    }

    return true;
}

bool TestCartesian() {
    Context c;
    if (!c.initializeFromFile("inp.wat.cartesian")) {
        std::cout << "Failed to initialize from inp.wat.cartesian" << std::endl;
        return false;
    }

	const auto Rounds = c.getRequiredNofRounds();
    std::vector<SimTK::Real> bl01, bl02; // bond lengths between atom 0-1 and 0-2
    std::vector<SimTK::Real> angles;
    std::vector<SimTK::Real> energies;

	bl01.reserve(Rounds);
	bl02.reserve(Rounds);
	angles.reserve(Rounds);
	energies.reserve(Rounds);

    for (int i = 0; i < Rounds; i++) {
        c.RunOneRound();

        const auto b01 = c.Distance(0, 0, 0, 0, 1);
        bl01.push_back(b01 * 10);

		const auto b02 = c.Distance(0, 0, 0, 0, 2);
		bl02.push_back(b02 * 10);

        const auto a = c.Roboangle(0, 0, 0, 0, 1, 2);
        angles.push_back(a);

        const auto e = c.getPotentialEnergy(0, 0);
        energies.push_back(e);
    }

    // std::ofstream fdistances("distances01.txt");
    // for (const auto d : bl01) {
    //     fdistances << d << std::endl;
    // }

	// std::ofstream fdistances2("distances02.txt");
	// for (const auto d : bl02) {
	// 	fdistances2 << d << std::endl;
	// }

    // std::ofstream fangles("angles.txt");
    // for (const auto a : angles) {
    //     fangles << a << std::endl;
    // }

    // test bond lengths for the the mobile atom pair
    if (!TestMaxwellBoltzmannDistribution(bl01, BOND_FORCE_CONSTANT, BONDLEN_EQ, BONDLEN_BEGIN, BONDLEN_END, NPOINTS, TEMPERATURE, "cartesian")) {
        return false;
    }

    if (!TestMaxwellBoltzmannDistribution(bl02, BOND_FORCE_CONSTANT, BONDLEN_EQ, BONDLEN_BEGIN, BONDLEN_END, NPOINTS, TEMPERATURE, "cartesian")) {
        return false;
    }

    if (!TestMaxwellBoltzmannDistribution(angles, ANGLE_FORCE_CONSTANT, ANGLE_EQ, ANGLE_BEGIN, ANGLE_END, NPOINTS, TEMPERATURE, "cartesian")) {
        return false;
    }

    return true;
}

// return 0 for success
// return 1 for fail
int main() {
    // if (!TestSlider01()) {
    //     return 1;
    // }

    // TOOO what about skewness? KolmogorovSmirnov? energy? dcd?

    if (!TestSlider01BS02()) {
        return 1;
    }

    // if (!TestCartesian()) {
    //     return 1;
    // }

    return 0;
}