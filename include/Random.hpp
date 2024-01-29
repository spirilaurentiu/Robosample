#pragma once

#include "Robo.hpp"

using Random32 = pcg32_k1024;
using Random64 = pcg64_k1024;

Random32 buildRandom32(uint32_t seed);
Random64 buildRandom64(Random32& seeder);

class RANDOM_CACHE {
public:
	void initialize(int argNU);
    int getNU() const;
    void setSeed(Random32& seeder);
    void generateGaussianVelocities();
    void wait();
    const SimTK::Vector& getV() const;

private:
	std::normal_distribution<> Gaussian = std::normal_distribution<>(0.0, 1.0);
	pcg64_k1024 Random64;

	std::function<SimTK::Real()> GenerateGaussian = [this]() mutable {
		return this->Gaussian(this->Random64);
	};

	std::function<void()> FillWithGaussian = [this]() mutable {
		std::generate(V.begin(), V.end(), GenerateGaussian);
	};

	std::future<void> task;

	SimTK::Vector V;
	int nu = -1;
};
