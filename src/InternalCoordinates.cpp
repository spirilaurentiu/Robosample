#include "InternalCoordinates.hpp"

AmberAtom::AmberAtom(const bSpecificAtom& b) {
	mass = b.mass;
	bonds = static_cast<int>(b.bondsInvolved.size());
	amberId = b.number;
}

AmberAtom::AmberAtom(const bSpecificAtom* b) {
	mass = b->mass;
	bonds = static_cast<int>(b->bondsInvolved.size());
	amberId = b->number;
}

bool AmberAtom::operator==(const AmberAtom& rhs) {
	return amberId == rhs.amberId;
}

bool AmberAtom::operator!=(const AmberAtom& rhs) {
	return !operator==(rhs);
}

void InternalCoordinates::compute(const std::vector<bSpecificAtom>& bAtomList) {
    computeRoot(bAtomList);

	selectedAtoms.reserve(bAtomList.size());
	indexMap = std::vector<int>(bAtomList.size(), -1);

	selectAtom(getRoot().first);
	selectAtom(getRoot().second);
	selectAtom(getRoot().third);

	while (selectedAtoms.size() < bAtomList.size()) {
		for (int i = 0; i < selectedAtoms.size(); i++) {
			const auto a1 = selectedAtoms[i];

			a0_list.clear();
			a0_list.reserve(bAtomList[a1].neighbors.size());
			for (const auto& b : bAtomList[a1].neighbors) {
				AmberAtom a0(b);
				if (!isSelected(a0))
					a0_list.push_back(a0);
			}
			sortByMass(a0_list, false);

			for(const auto& a0 : a0_list) {

				a2_list.clear();
				a2_list.reserve(bAtomList[a1].neighbors.size() - 1);
				for (const auto& b : bAtomList[a1].neighbors) {
					AmberAtom a2(b);
					if (isSelected(a2) && !isTerminal(a2) && a2 != a0) {
						a2_list.push_back(a2);
					}
				}
				sortByMass(a2_list, false);
				if(a2_list.empty()) continue;
				const auto& a2 = a2_list[0];

				a3_list.clear();
				a3_list.reserve(bAtomList[a2.amberId].neighbors.size());
				for (const auto& b : bAtomList[a2.amberId].neighbors) {
					AmberAtom a3(b);
					if (isSelected(a3) && a3.amberId != a1) {
						a3_list.push_back(a3);
					}
				}
				sortByMass(a3_list, false);
				if(a3_list.empty()) continue;
				const auto& a3 = a3_list[0];

				selectAtom(a0);
				bonds.push_back({ a0.amberId, a1 });
				angles.push_back({ a0.amberId, a1, a2.amberId });
				torsions.push_back({ a0.amberId, a1, a2.amberId, a3.amberId });
			}
		}
	}
}

const TORSION& InternalCoordinates::getRoot() const {
    return root;
}

const std::vector<BOND>& InternalCoordinates::getBonds() const {
	return bonds;
}

const std::vector<ANGLE>& InternalCoordinates::getAngles() const {
	return angles;
}

const std::vector<TORSION>& InternalCoordinates::getTorsions() const {
	return torsions;
}

int InternalCoordinates::amber2BAT(int amberIx) const {
	return indexMap[amberIx];
}

int InternalCoordinates::BAT2amber(int bat) const {
	return selectedAtoms[bat];
}

void InternalCoordinates::computeRoot(const std::vector<bSpecificAtom>& bAtomList) {
	// get terminal atoms
	std::vector<AmberAtom> terminalAtoms;
	terminalAtoms.reserve(bAtomList.size());

	for (int i=0; i < bAtomList.size(); i++) {
		AmberAtom a(bAtomList[i]);
		if (isTerminal(a)) {
			terminalAtoms.push_back(a);
		}
	}

	// prioritize by mass and amber id
	sortByMass(terminalAtoms, true);

	// find root
	root.first = terminalAtoms[0].amberId;
	root.second = bAtomList[root.first].neighbors[0]->number;

	if (bAtomList[root.second].neighbors.size() == 1) {
		root.third = bAtomList[root.second].neighbors[0]->number;
	} else {
		std::vector<AmberAtom> candidates;
		for (const auto& b : bAtomList[root.second].neighbors) {
			AmberAtom a(b);
			if (!isTerminal(a)) {
				candidates.push_back(a);
			}
		}
		sortByMass(candidates, true);
		root.third = candidates[0].amberId;
	}
}

void InternalCoordinates::sortByMass(std::vector<AmberAtom>& v, bool reverse) {
	std::sort(v.begin(), v.end(), [reverse](const AmberAtom& lhs, const AmberAtom& rhs) {
		if (lhs.mass == rhs.mass) {
			if (reverse) return lhs.amberId > rhs.amberId;
			return lhs.amberId < rhs.amberId;
		}

		if (reverse) return lhs.mass > rhs.mass;
		return lhs.mass < rhs.mass;
	});
};

bool InternalCoordinates::isSelected(const AmberAtom& a) const {
	return indexMap[a.amberId] != -1;
}

void InternalCoordinates::selectAtom(const AmberAtom& a) {
	selectedAtoms.push_back(a.amberId);
	indexMap[a.amberId] = selectedAtoms.size() - 1;
}

void InternalCoordinates::selectAtom(int a) {
	selectedAtoms.push_back(a);
	indexMap[a] = selectedAtoms.size() - 1;
}

bool InternalCoordinates::isTerminal(const AmberAtom& a) const {
	return a.bonds == 1;
}