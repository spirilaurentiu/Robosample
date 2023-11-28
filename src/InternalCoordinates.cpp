#include "InternalCoordinates.hpp"

bool BOND::operator==(const BOND& rhs) const {
	return first == rhs.first && second == rhs.second;
}

bool BOND::operator!=(const BOND& rhs) const {
	return !operator==(rhs);
}

bool ANGLE::operator==(const ANGLE& rhs) const {
	return first == rhs.first && second == rhs.second && third == rhs.third;
}

bool ANGLE::operator!=(const ANGLE& rhs) const {
	return !operator==(rhs);
}

bool TORSION::operator==(const TORSION& rhs) const {
	return first == rhs.first && second == rhs.second && third == rhs.third && fourth == rhs.fourth;
}

bool TORSION::operator!=(const TORSION& rhs) const {
	return !operator==(rhs);
}

AmberAtom::AmberAtom(const bSpecificAtom& b) {
	mass = b.mass;
	bonds = static_cast<int>(b.bondsInvolved.size());
	amberId = b.getNumber();
}

AmberAtom::AmberAtom(const bSpecificAtom* b) {
	mass = b->mass;
	bonds = static_cast<int>(b->bondsInvolved.size());
	amberId = b->getNumber();
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

	selectAtom(root.first);
	selectAtom(root.second);
	selectAtom(root.third);

	// add root bonds
	bonds.push_back({ root.first, root.second });
	bonds.push_back({ root.second, root.third });

	while (selectedAtoms.size() < bAtomList.size()) {
		for (int i = 0; i < selectedAtoms.size(); i++) {
			// a1 is the base atom
			const auto a1 = selectedAtoms[i];

			// a0 is a new atom connected to a1
			a0_list.clear();
			a0_list.reserve(bAtomList[a1].neighbors.size());
			for (const auto& b : bAtomList[a1].neighbors) {
				AmberAtom a0(b);
				if (!isSelected(a0))
					a0_list.push_back(a0);
			}
			sortByMass(a0_list, false);

			// iterate all a0 atoms
			for(const auto& a0 : a0_list) {

				// a2 is an atom connected to a1, is not a terminal atom and has been selected
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

				// a3 is connected to a2, has been selected and is not a1
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

				// created a new torsion
				selectAtom(a0);
				bonds.push_back({ a0.amberId, a1 });
				angles.push_back({ a0.amberId, a1, a2.amberId });
				torsions.push_back({ a0.amberId, a1, a2.amberId, a3.amberId });
			}
		}
	}

	computeLevelsAndOffsets(bAtomList);
}

struct DEPTH_NODE {
	int i = -1;
	int depth = -1;
	int parent = -1;
};

bool equal_bond(const BOND& lhs, const BOND& rhs) {
	return (lhs.first == rhs.first && lhs.second == rhs.second) ||
		(lhs.first == rhs.second && lhs.second == rhs.first);
}

void InternalCoordinates::computeLevelsAndOffsets(const std::vector<bSpecificAtom>& bAtomList) {
	// adj list is bonds
	// node values are i&j from adj list
	// create vector<vector<int>> v; v[level][offset] is amber id (bAtomList indexing)
	// getOrderedAtoms() -> bAtomList indexing
	// getOrderedBonds() -> order of bonds as in orderedAtoms; i and j are from bAtomList

	// iterate original bonds
	for (int i = 0; i < bAtomList.size(); i++) {
		for (const auto& b : bAtomList[i].neighbors) {
			const int j = b->getNumber();

			if (j >= i) continue;

			BOND b1 = { i, j };
			bool found = false;

			for (const auto& b0 : bonds) {
				if (equal_bond(b0, b1)) {
					found = true;
					// cout << "bond " << b0.first << " " << b0.second << endl;
					break;
				}
			}

			if (!found) {
				ringClosingBonds.push_back(b1);
				// missing bond
				// std::cout << "missing bond " << b1.first << " " << b1.second << std::endl;
			}	
		}
	}

	// build adjacency lists for unordered graph
	std::vector<std::vector<int>> adjacency(bAtomList.size());
	adjacency[root.first].push_back(root.second);
	adjacency[root.second].push_back(root.first);
	adjacency[root.second].push_back(root.third);
	adjacency[root.third].push_back(root.second);

	for (const auto& b : bonds) {
		adjacency[b.first].push_back(b.second);
		adjacency[b.second].push_back(b.first);
	}

	// compute depth of the molecule graph using dfs
	// we start from the first atom
	std::stack<DEPTH_NODE> nodesToVisit;
	nodesToVisit.push({ root.first, 0, -1 });

	// store levels and offsets
	std::vector<std::vector<AmberAtom>> depthGraph;
	std::vector<bool> v(bAtomList.size(), false);

	// non recursive bfs
	while (!nodesToVisit.empty()) {
		const auto atom = nodesToVisit.top();
		nodesToVisit.pop();

		if (v[atom.i])
			continue;
		v[atom.i] = true;

		int depth = atom.depth;
		if (depthGraph.size() < depth + 1)
			depthGraph.push_back({});

		AmberAtom a(bAtomList[atom.i]);
		a.parent = atom.parent;
		depthGraph[depth].emplace_back(a);

		for (const int i : adjacency[atom.i])
			nodesToVisit.push({ i, depth + 1, atom.i });
	}

	// sort by mass (ligher first). if equal, sort by amber index (lower first)
	for (auto& level : depthGraph)
		sortByMass(level, false);

	// copy the graph to a more managable structure
	for (int level = 0; level < depthGraph.size(); level++) {

		// allocate space for the level
		levelGraph.push_back({});
		graphedAtoms.push_back({});

		if (level > 0) {
			graphedBonds.push_back({});
		}

		for (int off = 0; off < depthGraph[level].size(); off++) {
			levelGraph[level].push_back({depthGraph[level][off].amberId, depthGraph[level][off].parent});
			graphedAtoms[level].push_back(levelGraph[level][off].amberIndex);

			const auto& a = levelGraph[level][off];
			if (a.parent != -1) {
				graphedBonds[level - 1].push_back({ a.amberIndex, a.parent });
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

const std::vector<BOND>& InternalCoordinates::getRingClosingBonds() const {
	return ringClosingBonds;
}

const std::vector<ANGLE>& InternalCoordinates::getAngles() const {
	return angles;
}

const std::vector<TORSION>& InternalCoordinates::getTorsions() const {
	return torsions;
}

const std::vector<std::vector<BAT_ATOM>>& InternalCoordinates::getLevelGraph() const {
	return levelGraph;
}

const std::vector<std::vector<int>>& InternalCoordinates::getGraphedAtoms() const {
	return graphedAtoms;
}

const std::vector<std::vector<BOND>>& InternalCoordinates::getGraphedBonds() const {
	return graphedBonds;
}

int InternalCoordinates::amber2BAT(int amberIx) const {
	return indexMap[amberIx];
}

int InternalCoordinates::BAT2amber(int bat) const {
	return selectedAtoms[bat];
}

void InternalCoordinates::computeRoot(const std::vector<bSpecificAtom>& bAtomList) {
	// the root consists of three atoms
	// first is the heaviest terminal atom with the largest amber id
	// second is bonded to the first atom (there is only one bonded since first is terminal)
	// third is the heaviest atom bonded to second. if there are more than three, it cannot be terminal

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
	root.second = bAtomList[root.first].neighbors[0]->getNumber();

	if (bAtomList[root.second].neighbors.size() == 1) {
		root.third = bAtomList[root.second].neighbors[0]->getNumber();
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