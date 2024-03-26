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
	bonds = b.getNumBonds();
	amberId = b.getNumber();
}

AmberAtom::AmberAtom(const bSpecificAtom* b) {
	mass = b->mass;
	bonds = b->getNumBonds();
	amberId = b->getNumber();
}

bool AmberAtom::operator==(const AmberAtom& rhs) {
	return amberId == rhs.amberId;
}

bool AmberAtom::operator!=(const AmberAtom& rhs) {
	return !operator==(rhs);
}



/*!
 * <!-- This function computes the BAT graph and calls computeLevelsAndOffsets
 * to calculate it into a level / offset format -->
*/
void InternalCoordinates::compute(std::vector<bSpecificAtom>& bAtomList) {

    computeRoot(bAtomList);

	computeBAT( bAtomList );

	// TODO the results of this function are not used anywhere
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

/*!
 * <!-- This function computes the root -->
*/
bool InternalCoordinates::computeRoot(
	const std::vector<bSpecificAtom>& bAtomList)
{
	// the root consists of three atoms
	// first is the heaviest terminal atom with the largest amber id
	// second is bonded to the first atom (there is only one bonded since first is terminal)
	// third is the heaviest atom bonded to second.
	// if there are more than three, it cannot be terminal

	// get terminal atoms
	std::vector<AmberAtom> terminalAtoms;
	terminalAtoms.reserve(bAtomList.size());

	for (const auto& atom : bAtomList) {
		if (atom.getVisited()) { // Laurentiu
			continue;
		}

		AmberAtom a(atom);
		if (isTerminal(a)) {
			terminalAtoms.push_back(a);
		}
	}

	//  && bAtomList.size() > 1
	if (terminalAtoms.size() == 0) {
		return false;
	}

	// Create this root
	roots.push_back(root); // Laurentiu

	// First root atom
	sortByMass(terminalAtoms, true);
	roots.back().first = terminalAtoms[0].amberId;
	const auto& firstAtom = bAtomList[roots.back().first];

	// Second root atom
	if (firstAtom.getNumBonds() == 0) {
		return true;
	}
	roots.back().second = bAtomList[firstAtom.neighborsIndex[0]].getNumber();
	const auto& secondAtom = bAtomList[roots.back().second];

	// Third root atom
	if (secondAtom.getNumBonds() == 0) {
		return true;
	}
	else if (secondAtom.getNumBonds() == 1) {
		roots.back().third = bAtomList[secondAtom.neighborsIndex[0]].getNumber();
	} else {
		std::vector<AmberAtom> candidates;
		for (auto i : secondAtom.neighborsIndex) {
			AmberAtom a(bAtomList[i]);
			if (!isTerminal(a)) {
				candidates.push_back(a);
			}
		}
		sortByMass(candidates, true);
		roots.back().third = candidates[0].amberId;
	}

	return true;
}


/*!
 * <!-- This function computes the BAT graph and calls computeLevelsAndOffsets
 * to calculate it into a level / offset format -->
*/
void InternalCoordinates::computeBAT(std::vector<bSpecificAtom>& bAtomList) {

    selectedAtoms.reserve(bAtomList.size());
    indexMap = std::vector<int>(bAtomList.size(), -1);

	perMolBonds.push_back(std::vector<BOND>());     // Laurentiu
    std::vector<BOND>& bonds = perMolBonds.back();  // Laurentiu
    perMolAngles.push_back(std::vector<ANGLE>());       // Laurentiu
    std::vector<ANGLE>& angles = perMolAngles.back();   // Laurentiu
    perMolTorsions.push_back(std::vector<TORSION>());       // Laurentiu
    std::vector<TORSION>& torsions = perMolTorsions.back(); // Laurentiu

	if (roots.back().first != -1) {
		selectAtom(roots.back().first);
		bAtomList[roots.back().first].setParentNumber(-1);
	}
    
	if (roots.back().second != -1) {
		selectAtom(roots.back().second);
		bAtomList[roots.back().second].setParentNumber(roots.back().first);
		bonds.push_back({ roots.back().second, roots.back().first });
	}

	if (roots.back().third != -1) {
		selectAtom(roots.back().third);
		bAtomList[roots.back().third].setParentNumber(roots.back().second);
		bonds.push_back({ roots.back().third, roots.back().second });
	}

    std::size_t lastSelectedAtomsSize = 0;                      // Laurentiu

	while (selectedAtoms.size() < bAtomList.size() && bAtomList.size() >= 3) {

		if(!(selectedAtoms.size() > lastSelectedAtomsSize)){	// Laurentiu
			break;												// Laurentiu
		}														// Laurentiu
		lastSelectedAtomsSize = selectedAtoms.size(); 			// Laurentiu

		// a1 is the base atom
		for (auto a1 : selectedAtoms) {

			// a0 is a new atom connected to a1
			a0_list.clear();
			a0_list.reserve(bAtomList[a1].neighborsIndex.size());
			for (auto a0_index : bAtomList[a1].neighborsIndex) {
				AmberAtom a0(bAtomList[a0_index]);
				if (!isSelected(a0))
					a0_list.push_back(a0);
			}
			sortByMass(a0_list, false);

			// iterate all a0 atoms
			for(const auto& a0 : a0_list) {

				// a2 is an atom connected to a1, is not a terminal atom and has been selected
				a2_list.clear();
				a2_list.reserve(bAtomList[a1].neighborsIndex.size() - 1);
				for (auto a1_index : bAtomList[a1].neighborsIndex) {
					AmberAtom a2(bAtomList[a1_index]);
					if (isSelected(a2) && !isTerminal(a2) && a2 != a0) {
						a2_list.push_back(a2);
					}
				}
				sortByMass(a2_list, false);
				if(a2_list.empty()) continue;
				const auto& a2 = a2_list[0];

				// a3 is connected to a2, has been selected and is not a1
				a3_list.clear();
				a3_list.reserve(bAtomList[a2.amberId].neighborsIndex.size());
				for (auto a3_index : bAtomList[a2.amberId].neighborsIndex) {
					AmberAtom a3(bAtomList[a3_index]);
					if (isSelected(a3) && a3.amberId != a1) {
						a3_list.push_back(a3);
					}
				}
				sortByMass(a3_list, false);
				if(a3_list.empty()) continue;
				const auto& a3 = a3_list[0];

				// created a new torsion
				selectAtom(a0);
				BOND bond;
				bond.first = a0.amberId;
				bond.second = a1;
				bonds.emplace_back(bond);

				bAtomList[a0.amberId].setParentNumber(a1);

				//angles.push_back({ a0.amberId, a1, a2.amberId });
				//torsions.push_back({ a0.amberId, a1, a2.amberId, a3.amberId });
			}
		}
	}

	for(size_t cnt = 0; cnt < bonds.size(); cnt++){

		BOND& bond  = bonds[cnt];

		int child    = bond.first;
		int parent   = bond.second;
		int gparent  = bAtomList[parent].getParentNumber();
		int ggparent = -2;
		if(gparent >= 0){
			ggparent = bAtomList[gparent].getParentNumber();
		}

		angles.push_back({gparent, parent, child});
		torsions.push_back({ggparent, gparent, parent, child});

		// std::cout << "BAT pushed "
		// 	<< ggparent << " " << gparent << " " << parent << " "   << child 
		// << eol;

	}

}


/*!
 * <!-- This function computes the BAT graph and calls computeLevelsAndOffsets
 * to calculate it into a level / offset format -->
*/
void InternalCoordinates::computeBATOld(std::vector<bSpecificAtom>& bAtomList) {

	// selectedAtoms.reserve(bAtomList.size());
	// indexMap = std::vector<int>(bAtomList.size(), -1);

	// selectAtom(root.first);
	// selectAtom(root.second);
	// selectAtom(root.third);


	// perMolBonds.push_back(std::vector<BOND>());		// Laurentiu
	// std::vector<BOND>& bonds = perMolBonds.back();	// Laurentiu
	// perMolAngles.push_back(std::vector<ANGLE>());		// Laurentiu
	// std::vector<ANGLE>& angles = perMolAngles.back();	// Laurentiu
	// perMolTorsions.push_back(std::vector<TORSION>());		// Laurentiu
	// std::vector<TORSION>& torsions = perMolTorsions.back();	// Laurentiu

	// // add root bonds
	// bonds.push_back({ root.second, root.first });
	// bonds.push_back({ root.third, root.second });

	// std::size_t lastSelectedAtomsSize = 0;						// Laurentiu

	// while (selectedAtoms.size() < bAtomList.size()) {

	// 	if(!(selectedAtoms.size() > lastSelectedAtomsSize)){	// Laurentiu
	// 		break;												// Laurentiu
	// 	}														// Laurentiu
	// 	lastSelectedAtomsSize = selectedAtoms.size(); 			// Laurentiu

	// 	for (int i = 0; i < selectedAtoms.size(); i++) {
	// 		// a1 is the base atom
	// 		const auto a1 = selectedAtoms[i];

	// 		// a0 is a new atom connected to a1
	// 		a0_list.clear();
	// 		a0_list.reserve(bAtomList[a1].neighbors.size());
	// 		for (const auto& b : bAtomList[a1].neighbors) {
	// 			AmberAtom a0(b);
	// 			if (!isSelected(a0))
	// 				a0_list.push_back(a0);
	// 		}
	// 		sortByMass(a0_list, false);

	// 		// iterate all a0 atoms
	// 		for(const auto& a0 : a0_list) {

	// 			// a2 is an atom connected to a1, is not a terminal atom and has been selected
	// 			a2_list.clear();
	// 			a2_list.reserve(bAtomList[a1].neighbors.size() - 1);
	// 			for (const auto& b : bAtomList[a1].neighbors) {
	// 				AmberAtom a2(b);
	// 				if (isSelected(a2) && !isTerminal(a2) && a2 != a0) {
	// 					a2_list.push_back(a2);
	// 				}
	// 			}
	// 			sortByMass(a2_list, false);
	// 			if(a2_list.empty()) continue;
	// 			const auto& a2 = a2_list[0];

	// 			// a3 is connected to a2, has been selected and is not a1
	// 			a3_list.clear();
	// 			a3_list.reserve(bAtomList[a2.amberId].neighbors.size());
	// 			for (const auto& b : bAtomList[a2.amberId].neighbors) {
	// 				AmberAtom a3(b);
	// 				if (isSelected(a3) && a3.amberId != a1) {
	// 					a3_list.push_back(a3);
	// 				}
	// 			}
	// 			sortByMass(a3_list, false);
	// 			if(a3_list.empty()) continue;
	// 			const auto& a3 = a3_list[0];

	// 			// created a new torsion
	// 			selectAtom(a0);
	// 			bonds.push_back({ a0.amberId, a1 });
	// 			angles.push_back({ a0.amberId, a1, a2.amberId });
	// 			torsions.push_back({ a0.amberId, a1, a2.amberId, a3.amberId });
	// 		}
	// 	}
	// }

}

void InternalCoordinates::updateVisited(std::vector<bSpecificAtom>& bAtomList)
{
	std::size_t cnt = 0;
	std::vector<int>::iterator imIt;
	//for(imIt = indexMap.begin(); imIt != indexMap.end(); imIt++, cnt++)
	for(cnt = 0; cnt < bAtomList.size(); cnt++)
	{
		if(indexMap[cnt] != -1){
			bAtomList[cnt].setVisited(1);
		}
	}
}

/*!
 * <!-- Calculate level / offset format from the atom list read from Amber 
 * prmtop and deposit in levelGraph graphedAtoms
 * @todo what is levelGraph and graphedAtoms -->
*/
void InternalCoordinates::computeLevelsAndOffsets(const std::vector<bSpecificAtom>& bAtomList) {
	// adj list is bonds
	// node values are i&j from adj list
	// create vector<vector<int>> v; v[level][offset] is amber id (bAtomList indexing)
	// getOrderedAtoms() -> bAtomList indexing
	// getOrderedBonds() -> order of bonds as in orderedAtoms; i and j are from bAtomList

	if (bAtomList.size() == 1) {
		std::cerr << "[ERROR]: No atoms in the molecule" << std::endl;
		return;
	}

	std::vector<BOND> allBonds;
	for(auto && v : perMolBonds){
		allBonds.insert(allBonds.end(), v.begin(), v.end());
	}

	// iterate original bonds
	for (int i = 0; i < bAtomList.size(); i++) {
		for (auto smth : bAtomList[i].neighborsIndex) {
			const int j = bAtomList[smth].getNumber();

			if (j >= i) continue;

			BOND b1 = { i, j };
			bool found = false;

			for (const auto& b0 : allBonds) {
				// std::cout << "comparing " << b0.first << " " << b0.second << " to " << b1.first << " " << b1.second << std::endl;
				if (equal_bond(b0, b1)) {
					found = true;
					// cout << "found bond " << b0.first << " " << b0.second << endl;
					break;
				}
			}

			if (!found) {
				ringClosingBonds.push_back(b1);
				// // missing bond
				// std::cout << "missing bond " << b1.first << " " << b1.second << std::endl;
			}	
		}
	}

	// for (const auto& b : ringClosingBonds) {
	// 	std::cout << "ringClosingBonds " << b.first << " " << b.second << std::endl;
	// }

	// build adjacency lists for unordered graph
	std::vector<std::vector<int>> adjacency(bAtomList.size());

	// TODO
	// for (const auto& root : roots) {
	// }

	const auto first = roots.front().first;
	const auto second = roots.front().second;
	const auto third = roots.front().third;

	if (first != -1 && second != -1) {
		adjacency[first].push_back(second);
		adjacency[second].push_back(first);
	}

	if (second != -1 && third != -1) {
		adjacency[second].push_back(third);
		adjacency[third].push_back(second);
	}

	// adjacency[roots.front().first].push_back(roots.front().second);
	// adjacency[roots.front().second].push_back(roots.front().first);
	// adjacency[roots.front().second].push_back(roots.front().third);
	// adjacency[roots.front().third].push_back(roots.front().second);

	for (const auto& b : perMolBonds.back()) {
		adjacency[b.first].push_back(b.second);
		adjacency[b.second].push_back(b.first);
	}

	// compute depth of the molecule graph using dfs
	// we start from the first atom
	std::stack<DEPTH_NODE> nodesToVisit;
	nodesToVisit.push({ roots.front().first, 0, -1 });

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
			levelGraph[level].push_back({int(depthGraph[level][off].amberId), int(depthGraph[level][off].parent)});
			graphedAtoms[level].push_back(levelGraph[level][off].amberIndex);

			const auto& a = levelGraph[level][off];
			if (a.parent != -1) {
				graphedBonds[level - 1].push_back({ a.amberIndex, a.parent });
			}
		}
	}
}

/*!
 * <!-- Get the last root added -->
*/
const TORSION& InternalCoordinates::getLastRoot() const {
    return root;
}

/*!
 * <!-- Return the root of the whichMoleculeth molecule -->
*/
const TORSION& InternalCoordinates::getRoot( int which ) const
{

	return roots[which];

}

/*!
 * <!-- Get all the roots -->
*/
const std::vector<TORSION>& InternalCoordinates::getRoots() const
{
	return roots;
}

/*!
 * <!-- Get the last molecule's bond vector added -->
*/
const std::vector<BOND>& InternalCoordinates::getLastMolsBonds() const {
	return perMolBonds.back();
}

/*!
 * <!-- Return a particular molecule's bond list -->
*/
const std::vector<BOND>& InternalCoordinates::getMoleculeBonds( int which ) const
{
	return perMolBonds[which];
}

/*!
 * <!-- Return the BOND 2D vector -->
*/	
const std::vector<std::vector<BOND>>& InternalCoordinates::getBonds() const
{
	return perMolBonds;
}

const std::vector<BOND>& InternalCoordinates::getRingClosingBonds() const {
	return ringClosingBonds;
}

const std::vector<ANGLE>& InternalCoordinates::getAngles() const {
	return perMolAngles.back();
}

const std::vector<TORSION>& InternalCoordinates::getTorsions() const {
	return perMolTorsions.back();
}

/*!
 * <!-- Function to search for a BOND with a specific 'second' value -->
*/
int
InternalCoordinates::findBondByFirst(
	int wIx,
	int targetFirstValue)
{

	std::vector<BOND>::iterator it =
		std::find_if(
			perMolBonds[wIx].begin(),
			perMolBonds[wIx].end(),
			[targetFirstValue](const BOND& someBond)
			{
				return someBond.first == targetFirstValue;
			}
		);

	if(it == perMolBonds[wIx].end()){
		return -1;
	}else{
		return std::distance(perMolBonds[wIx].begin(), it);
	}
}

/*!
 * <!--  -->
*/
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

/*!
 * <!-- Sort atoms in place first by mass and then by indeces -->
*/
void InternalCoordinates::sortByMass(std::vector<AmberAtom>& v, bool reverse)
{
	// Sort atoms first by mass and then by indeces
	std::sort(v.begin(), v.end(),
		[reverse](const AmberAtom& lhs, const AmberAtom& rhs)
	{
		// If masses are equal sort by indeces
		if (lhs.mass == rhs.mass) {
			if (reverse){
				return lhs.amberId > rhs.amberId;
			}else{
				return lhs.amberId < rhs.amberId;
			}
		}

		// Sort mass
		if (reverse) {
			return lhs.mass > rhs.mass;
		}else{
			return lhs.mass < rhs.mass;
		}

	});
};

/*!
 * <!-- -->
*/
bool InternalCoordinates::isSelected(const AmberAtom& a) const {
	return indexMap[a.amberId] != -1;
}

/*!
 * <!-- -->
*/
void InternalCoordinates::selectAtom(const AmberAtom& a) {
	selectedAtoms.push_back(a.amberId);
	indexMap[a.amberId] = selectedAtoms.size() - 1;
}

/*!
 * <!-- -->
*/
void InternalCoordinates::selectAtom(int a) {
	selectedAtoms.push_back(a);
	indexMap[a] = selectedAtoms.size() - 1;
}

/*!
 * <!-- -->
*/
bool InternalCoordinates::isTerminal(const AmberAtom& a) const {
	return a.bonds == 0 || a.bonds == 1;
}