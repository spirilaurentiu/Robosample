#include "ElementCache.hpp"

bool ELEMENT_CACHE::addElement(int atomicNumber, SimTK::mdunits::Mass atomicMass) {
	if (atomicNumber < 1 || atomicNumber > 118) {
		std::cerr << "Invalid atomic number " << atomicNumber << ". It must be between 1 and 118." << std::endl;
		return false;
	}

	cache.insert(std::make_pair(key(atomicNumber, atomicMass),
		SimTK::Element(
			atomicNumber,
			elements[atomicNumber - 1].first.c_str(),
			elements[atomicNumber - 1].second.c_str(),
			atomicMass
	)));

	return true;
}

const SimTK::Element::Symbol& ELEMENT_CACHE::getSymbolByAtomicNumber(int atomicNumber) const {
	if (atomicNumber < 1 || atomicNumber > 118) {
		std::cerr << "Invalid atomic number " << atomicNumber << ". It must be between 1 and 118." << std::endl;
		return null;
	}

	return elements[atomicNumber - 1].second;
}

const SimTK::Element& ELEMENT_CACHE::getElement(int atomicNumber, SimTK::mdunits::Mass atomicMass) const {
	return cache.at(key(atomicNumber, atomicMass));
}

void ELEMENT_CACHE::print() const {
	for (const auto& e : cache) {
		std::cout << e.first << " " << e.second << std::endl;
	}
}

inline std::size_t ELEMENT_CACHE::key(int i, int j) const {
	return (std::size_t) i << 32 | (unsigned int) j;
}