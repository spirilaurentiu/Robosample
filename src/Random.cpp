#include "Random.hpp"

Random32 buildRandom32(uint32_t seed) {
	std::seed_seq seedSeq { seed, seed + 1, seed + 2 };
	std::array<uint32_t, 3> seeds;
	seedSeq.generate(seeds.begin(), seeds.end());

	const uint64_t seed_bits = (uint64_t(seeds[0]) << 32) | seeds[1];

	return Random32 { seed_bits };
}

Random64 buildRandom64(Random32& seeder) {
	const auto high = (uint64_t(seeder()) << 32) | (seeder() & 0xFFFFFFFF);
	const auto low = (uint64_t(seeder()) << 32) | (seeder() & 0xFFFFFFFF);

	return Random64 { PCG_128BIT_CONSTANT(high, low) };
}

void RANDOM_CACHE::initialize(int argNU) {
    nu = argNU;
    V.resize(nu);

    FillWithGaussian();
}

int RANDOM_CACHE::getNU() const {
    return nu;
}

void RANDOM_CACHE::setSeed(Random32& seeder) {
    Random64 = buildRandom64(seeder);
}

void RANDOM_CACHE::generateGaussianVelocities() {
    task = std::async(std::launch::async, FillWithGaussian);
}

void RANDOM_CACHE::wait() {
    task.wait();
}

const SimTK::Vector& RANDOM_CACHE::getV() const {
    return V;
}
