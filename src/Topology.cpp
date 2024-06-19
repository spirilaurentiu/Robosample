/**@file
Implementation of Topology class. **/

#include "Topology.hpp"

using namespace SimTK;

#ifndef TRACE_GRAPH
#define TRACE_GRAPH true
#endif

/** Default constructor.Sets the name of this molecule to 'no_name '.
The name has no particular function and is not guaranteed to be unique **/
Topology::Topology(){
	this->name = std::string("no_name");
	this->setCompoundName((this->name));
}

/** Constructor that sets the name of the molecule. The name has no particular
function and is not guaranteed to be unique **/
Topology::Topology(std::string nameOfThisMolecule){
	this->name = nameOfThisMolecule;
	this->setCompoundName((this->name));
}

/** Default destructor. It deallocates bAtomType of every atom in the bAtomList
because we want to allow the valence to change during the simulation
e.g. semi-grand canonical ensemble. **/
Topology::~Topology() {
	// for (auto& atom : bAtomList) {
	// 	atom.destroy();
	// }
}

/** Print atom and bonds list with details**/
void Topology::PrintAtomList(int whichWorld)
{
	// Atoms
	std::cout<<"Topology::PrintAtomList\n";
	for(unsigned int i = 0; i < subAtomList.size(); i++){
		subAtomList[i].Print(whichWorld);
	}

	// Bonds
	for(unsigned int i = 0; i < subBondList.size(); i++){
		subBondList[i].Print(whichWorld);
	}
}

/** The following functions are used to build the molecular graph using bonding
information from bonds list and bondsInvolved list of each atom in bAtomList.
**/
void Topology::generateAIx2TopXMaps(void)
{
	for (unsigned int aix = 0; aix < getNumAtoms(); ++aix) {
		aIx2TopTransform.insert(std::make_pair(
			(subAtomList[aix]).getCompoundAtomIndex(), SimTK::Transform()));
	}
}

/** Print Molmodel specific types as introduced in Gmolmodel **/
void Topology::PrintMolmodelAndDuMMTypes(
	SimTK::DuMMForceFieldSubsystem& dumm) const
{
	scout("Print Molmodel And DuMM Types:"); ceolf; 
	for (size_t sAIx = 0; sAIx < subAtomList.size(); ++sAIx){

		std::cout << " list ix " << sAIx
			<< " CompoundAtomIndex " << (subAtomList[sAIx]).getCompoundAtomIndex()
			//<< " DuMMAtomIndex " << getDuMMAtomIndex((subAtomList[i]).getCompoundAtomIndex())
			<< " biotypename " << (subAtomList[sAIx]).getBiotype()
			<< " name " << (subAtomList[sAIx]).getName()
			<< " BiotypeIndex " << (subAtomList[sAIx]).getBiotypeIndex()
			<< " ChargedAtomTypeIndex "<< (subAtomList[sAIx]).getChargedAtomTypeIndex()
			<< " AtomClassIx " << (subAtomList[sAIx]).getDummAtomClassIndex()
			<< " partialChargeInE " << (subAtomList[sAIx]).charge
			<< " chargedAtomTypeIndex "
			<< (subAtomList[sAIx]).getChargedAtomTypeIndex()
			<< " DuMM VdW Radius "
			<< dumm.getVdwRadius((subAtomList[sAIx]).getDummAtomClassIndex())
			<< " DuMM VdW Well Depth "
			<< dumm.getVdwWellDepth((subAtomList[sAIx]).getDummAtomClassIndex())
			<< std::endl << std::flush;
	}
}

bool Topology::checkIfTripleUnorderedAreEqual(
		std::vector<SimTK::Compound::AtomIndex> &first,
		std::vector<SimTK::Compound::AtomIndex> &second)
{
	assert(!"Deprecated function.");

	// if(first == second){
	// 	return true;
	// }
	// else if(
	// 		(first[0] == second[2]) &&
	// 		(first[1] == second[1]) &&
	// 		(first[2] == second[0])
	// 		){
	// 	return true;
	// }
	// else{
	// 	return false;
	// }

}

// Helper function for calcLogDetMBATAnglesContribution
// Finds all triple runs - TODO VERY INEFFICIENT
void Topology::loadTriples_SP_NEW()
{
	// // Assign Compound coordinates by matching bAtomList coordinates
	// std::map<AtomIndex, Vec3> atomTargets;
	// for(int ix = 0; ix < getNumAtoms(); ++ix){
	// 	Vec3 vec(subAtomList[ix].getX(),
	// 			 subAtomList[ix].getY(),
	// 			 subAtomList[ix].getZ());
	// 	atomTargets.insert(pair<AtomIndex, Vec3> (
	// 		subAtomList[ix].getCompoundAtomIndex(), vec));
	// }
	// std::vector< std::vector<Compound::AtomIndex> > bondedAtomRuns =
	// getBondedAtomRuns(3, atomTargets);
	// // Find root bAtomList index
	// int ix = -1;
	// for(const auto& atom: subAtomList){
	// 	ix++;
	// }
	// // Find neighbour with maximum atomIndex
	// int maxAIx = -1;
	// Compound::AtomIndex aIx;
	// for(auto atom: subAtomList[ix].neighbors){
	// 	aIx = atom->getCompoundAtomIndex();
	// 	if(aIx > maxAIx){
	// 		maxAIx = aIx;
	// 	}
	// }
	// int flag;
	// int bIx = -1;
	// for(auto bAR: bondedAtomRuns){ // Iterate bondedAtomRuns
	// 	bIx++;
	// 	flag = 0;
	// 	for(auto tripleEntry: triples){ // Iterate triples gathered so far
	// 		if(checkIfTripleUnorderedAreEqual(bAR, tripleEntry)){
	// 			flag = 1;
	// 			break;
	// 		}
	// 	} // END Iterate triples gathered so far
	// 	if(!flag){ // Not found in gathered triples
	// 		if((bAR[0] < bAR[1]) || (bAR[2] < bAR[1]) // Only level changing branches
	// 		|| ((bAR[1] == 0) && (bAR[2] == maxAIx)) // except for the root atom
	// 		){
	// 			triples.push_back(bAR);
	// 		}
	// 	}
	// } // END Iterate bondedAtomRuns


}

// Numerically unstable around -pi, 0 and pi due to the log(0)
SimTK::Real Topology::calcLogSineSqrGamma2(const SimTK::State &quatState)
{
	SimTK::Compound::AtomIndex aIx = subAtomList[rootAtomIx].getCompoundAtomIndex();
	SimTK::Transform X = calcAtomFrameInGroundFrame(quatState, aIx);
	SimTK::Quaternion quat = (X.R()).convertRotationToQuaternion();

	SimTK::Real w = quat[0];
	SimTK::Real x = quat[1];
	SimTK::Real y = quat[2];
	SimTK::Real z = quat[3];
	SimTK::Real sinPitch = 2 * ((w * y) - (z * x));

	//std::cout << std::setprecision(20) << std::fixed;
	SimTK::Real pitch = std::asin(sinPitch);
	//std::cout << "Topology pitch " << pitch << std::endl;

	SimTK::Real result = std::log(sinPitch * sinPitch);

	// Quick and dirty
	if(result < -14.0){ // Around double precision log(0)
		result = -14.0;
	}

	// More elaborate fix
	//TODO: Given a difference compute the limit of this log

	return result;
}


SimTK::Real Topology::calcLogDetMBATGamma2Contribution(const SimTK::State& quatState){
	//State& eulerState;
	//matter.convertToEulerAngles(quatState, eulerState);
	//std::cout << "calcLogDetMBATGamma2Contribution quaternionState " << quatState << std::endl;
	//std::cout << "calcLogDetMBATGamma2Contribution	  eulerState " << eulerState << std::endl;

	bSpecificAtom *root = &(subAtomList[bSpecificAtomRootIndex]);
	SimTK::Compound::AtomIndex aIx = root->getCompoundAtomIndex();
	SimTK::Transform X = calcAtomFrameInGroundFrame(quatState, aIx);
	SimTK::Quaternion quat = (X.R()).convertRotationToQuaternion();
	//std::cout << "calcLogDetMBATGamma2Contribution quaternion " << quat << std::endl;

	SimTK::Real w = quat[0];
	SimTK::Real x = quat[1];
	SimTK::Real y = quat[2];
	SimTK::Real z = quat[3];
	SimTK::Real sinPitch = 2 * (w * y - z * x);

	std::cout << std::setprecision(20) << std::fixed;
	std::cout << "sinpitch " << sinPitch << std::endl;
	std::cout << "sinpitchsq" << sinPitch * sinPitch << std::endl;

	SimTK::Real pitch = std::asin(sinPitch);
	std::cout << "pitch " << pitch << std::endl;
	//if(pitch < 0){
	//	pitch = pitch + (2*SimTK_PI);
	//	std::cout << "sin converted pitch " << std::sin(pitch) << std::endl;
	//}

	if(sinPitch < SimTK::Eps){ // consider using SimTK::Eps
		return -SimTK::Infinity;
	}
	SimTK::Real result = std::log(sinPitch * sinPitch);
	return result;
}


SimTK::Real Topology::calcLogDetMBATDistsContribution(const SimTK::State&){
	// function args were const SimTK::State& someState

	// Assign Compound coordinates by matching bAtomList coordinates
	std::map<AtomIndex, Vec3> atomTargets;
	for(int ix = 0; ix < getNumAtoms(); ++ix){
		Vec3 vec(subAtomList[ix].getX(), subAtomList[ix].getY(), subAtomList[ix].getZ());
		atomTargets.insert(pair<AtomIndex, Vec3> (subAtomList[ix].getCompoundAtomIndex(), vec));
	}

	//std::cout << "Topology::calcLogDetMBATDistsContribution dists: " ;
	SimTK::Real result = 0.0;

	for(auto bond: subBondList){

		int relativeBond_I = bond.i - subBondList.get_offset();
		int relativeBond_J = bond.j - subBondList.get_offset();

		// scout("Topology::calcLogDetMBATDistsContribution_SP_NEW ") << bond.i <<" " << bond.j <<" " << subBondList.get_offset() <<" " << relativeBond_I <<" " << relativeBond_J <<" " << eol;

		SimTK::Vec3 atom1pos = SimTK::Vec3(subAtomList[relativeBond_I].getX(), subAtomList[relativeBond_I].getY(), subAtomList[relativeBond_I].getZ());
		SimTK::Vec3 atom2pos = SimTK::Vec3(subAtomList[relativeBond_J].getX(), subAtomList[relativeBond_J].getY(), subAtomList[relativeBond_J].getZ());

		//SimTK::Real distSqr = (atom2pos - atom1pos).normSqr(); // funny results ?
		SimTK::Real dist = std::sqrt(
				std::pow(atom2pos[0] - atom1pos[0], 2) +
				std::pow(atom2pos[1] - atom1pos[1], 2) +
				std::pow(atom2pos[2] - atom1pos[2], 2));

		//std::cout << "atom " << bond.j << " 2*logDistSqr " << std::log(distSqr) + std::log(distSqr) << " " ;
		//std::cout << "dist " << bond.j << " = " << dist << " " ;

		result = result + (4.0 * std::log(dist));

	}
	//std::cout << std::endl;

	return result;
}


SimTK::Real Topology::calcLogDetMBATAnglesContribution(const SimTK::State&){
	// function args were const SimTK::State& someState

	// Assign Compound coordinates by matching bAtomList coordinates
	std::map<AtomIndex, Vec3> atomTargets;
	for(int ix = 0; ix < getNumAtoms(); ++ix){
			Vec3 vec(subAtomList[ix].getX(), subAtomList[ix].getY(), subAtomList[ix].getZ());
			atomTargets.insert(pair<AtomIndex, Vec3> (subAtomList[ix].getCompoundAtomIndex(), vec));
	}

	//std::cout << "Topology::calcLogDetMBATAnglesContribution angles: " ;

	SimTK::Real result = 0.0;
	for(auto triple: triples){

		//spacecout("TRIPLE:", triple[0], triple[1], triple[2]);
		//SimTK::Vec3 vec0 = calcAtomLocationInGroundFrame(someState, triple[0]);
		//SimTK::Vec3 vec1 = calcAtomLocationInGroundFrame(someState, triple[1]);
		//SimTK::Vec3 vec2 = calcAtomLocationInGroundFrame(someState, triple[2]);

		SimTK::UnitVec3 v1(atomTargets.find(triple[0])->second - atomTargets.find(triple[1])->second);
		SimTK::UnitVec3 v2(atomTargets.find(triple[2])->second - atomTargets.find(triple[1])->second);

		SimTK::Real dotProduct = dot(v1, v2);
		assert(dotProduct < 1.1);
		assert(dotProduct > -1.1);
		if (dotProduct > 1.0) dotProduct = 1.0;
		if (dotProduct < -1.0) dotProduct = -1.0;
		// SimTK::Real angle = std::acos(dotProduct);
		//std::cout << SimTK_RADIAN_TO_DEGREE * angle << " " ;

		SimTK::Real sinSqAngle = 1 - (dotProduct * dotProduct);

		result = result + std::log(sinSqAngle);

	}

	//std::cout << std::endl;

	return result;
}


SimTK::Real Topology::calcLogDetMBATMassesContribution(const SimTK::State&)
{
	// function args were const SimTK::State& someState

	//std::cout << "Topology::calcLogDetMBATMassesContribution masses: " ;
	SimTK::Real result = 0.0;
	for(const auto& atom: subAtomList){
		//std::cout << 3.0 * std::log(atom.mass) << " " ;
		//std::cout << atom.mass << " " ;
		result += 3.0 * std::log(atom.mass);
	}
	//std::cout << std::endl;

	return result;
}


SimTK::Real Topology::calcLogDetMBATInternal(const SimTK::State& someState)
{
	SimTK::Real distsContribution = calcLogDetMBATDistsContribution(someState);
	SimTK::Real anglesContribution = calcLogDetMBATAnglesContribution(someState);
	SimTK::Real massesContribution = calcLogDetMBATMassesContribution(someState);

	// std::cout << std::setprecision(20) << std::fixed;
	// std::cout << "MBAT dists masses angles contributions: "
	//          << distsContribution << " "
	//          << massesContribution << " "
	//          << anglesContribution << std::endl;

	return distsContribution + anglesContribution + massesContribution;
}


/**
 ** Interface **
 **/

/** Get the number of atoms in the molecule **/
int Topology::getNAtoms() const{
	return getNumAtoms();
}

/*!
 * <!-- Get the number of bonds in the molecule -->
*/
int Topology::getNBonds() const{
	return subBondList.size();
}

/** Get a pointer to an atom object in the atom list inquiring
by its Molmodel assigned atom index (SimTK::Compound::AtomIndex) .**/
// TODO: Optimize use CompoundAtomIx2GmolAtomIx instead
bSpecificAtom * Topology::updAtomByAtomIx(int cAIx) {
	for (int aix = 0; aix < natoms; aix++){
		if(subAtomList[aix].getCompoundAtomIndex() == cAIx){
			return &subAtomList[aix];
		}
	}

	return nullptr;
}

/** Get a pointer to an atom object in the atom list inquiring
by atom name. **/
bSpecificAtom * Topology::getAtomByName(std::string) const {
	// function args were std::string name
	assert(!"Not implemented."); throw std::exception();

	return nullptr;
}

/** Get the neighbours in the graph. **/
std::vector<bSpecificAtom *> Topology::getNeighbours(int) const {
	assert(!"Not implemented."); throw std::exception();

	return {};
}

/* Check if a1 and a2 are bonded */
bool Topology::checkBond(int a1, int a2)
{
	for(int i = 0; i < nbonds; i++){
		if( (bonds[i]).isThisMe(a1, a2) )
		{
			return true;
		}
	}
	return false;
}

/** **/
const bBond& Topology::getBond(int a1, int a2) const
{
	// TODO is this fast enough?
	const auto bond = std::find_if(bonds.begin(), bonds.end(), [a1, a2](bBond b) { return b.isThisMe(a1, a2); });

	#ifndef NDEBUG
		std::string assert_string("No bond with these atom indeces found: " + to_string(a1) + " " + to_string(a2));
		assert(bond != bonds.end() && assert_string.size());
	#endif

	return *bond;
}

/** Create MobilizedBodyIndex vs Compound::AtomIndex maps **/
void Topology::loadAIx2MbxMap()
{

	// If the map is empty fill with empty vectors first
	if(aIx2mbx.empty()){
		
		// Iterate through atoms and get their MobilizedBodyIndeces
		for (unsigned int aCnt = 0; aCnt < getNumAtoms(); ++aCnt) {

			// Get atomIndex from atomList
			SimTK::Compound::AtomIndex aIx = (subAtomList[aCnt]).getCompoundAtomIndex();

			// Insert
			aIx2mbx.insert(
				std::pair< SimTK::Compound::AtomIndex, std::vector<SimTK::MobilizedBodyIndex> >
					(aIx, std::vector<SimTK::MobilizedBodyIndex>())
			);
		}

	}

	// Iterate through atoms and get their MobilizedBodyIndeces
	for (unsigned int aCnt = 0; aCnt < getNumAtoms(); ++aCnt) {

		// Get atomIndex from atomList
		SimTK::Compound::AtomIndex aIx = (subAtomList[aCnt]).getCompoundAtomIndex();

		// Get MobilizedBodyIndex from CompoundAtom
		SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndex(aIx);

		// Insert
		//aIx2mbx.insert(
		//		std::pair<SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex>
		//		(aIx, mbx));
		aIx2mbx[aIx].emplace_back(mbx);
	}
}

/** Compound AtomIndex to bAtomList number **/
void Topology::loadCompoundAtomIx2GmolAtomIx()
{
	for (unsigned int i = 0; i < getNumAtoms(); ++i) {
		SimTK::Compound::AtomIndex aIx = (subAtomList[i]).getCompoundAtomIndex();
		int gmolIx = (subAtomList[i]).getNumber();

		CompoundAtomIx2GmolAtomIx.insert(
			std::pair<SimTK::Compound::AtomIndex, int>
			(aIx, gmolIx));
	}
}

/**  **/
int Topology::getNumber(SimTK::Compound::AtomIndex cAIx)
{
	return CompoundAtomIx2GmolAtomIx[cAIx];
}

/*!
 * <!-- Calculate all atom frames in top frame. It avoids calling
 * calcDefaultAtomFrameInCompoundFrame multiple times. This has
 * to be called every time the coordinates change though. -->
*/
void Topology::calcAtomsTopTransforms()
{
	for (unsigned int i = 0; i < getNumAtoms(); ++i) {
		SimTK::Compound::AtomIndex aIx = (subAtomList[i]).getCompoundAtomIndex();
		aIx2TopTransform[aIx] = calcDefaultAtomFrameInCompoundFrame(aIx);
	}
}

/*!
 * <!--  -->
*/
void Topology::printTopTransforms()
{
	std::cout << "Topology TopTransforms " << std::endl;
	for (unsigned int i = 0; i < getNumAtoms(); ++i) {
		SimTK::Compound::AtomIndex aIx = (subAtomList[i]).getCompoundAtomIndex();
		std::cout << aIx << " " << aIx2TopTransform[aIx] << std::endl;
	}
}

/*!
 * <!-- Get atom Top level transform from the existing Topology map -->
*/
SimTK::Transform Topology::getTopTransform_FromMap(SimTK::Compound::AtomIndex aIx)
{
	return aIx2TopTransform[aIx];
}

// Return mbx by calling DuMM functions
SimTK::MobilizedBodyIndex Topology::getAtomMobilizedBodyIndexThroughDumm(
	SimTK::Compound::AtomIndex aIx,
	SimTK::DuMMForceFieldSubsystem& dumm)
{
	SimTK::DuMM::AtomIndex dAIx = getDuMMAtomIndex(aIx);
	return dumm.getAtomBody(dAIx);
}

// Get atom location on mobod through DuMM functions
SimTK::Vec3 Topology::getAtomLocationInMobilizedBodyFrameThroughDumm(
	SimTK::Compound::AtomIndex aIx,
	SimTK::DuMMForceFieldSubsystem& dumm)
{
	SimTK::DuMM::AtomIndex dAIx = getDuMMAtomIndex(aIx);
	return dumm.getAtomStationOnBody(dAIx);
}

SimTK::Vec3 Topology::calcAtomLocationInGroundFrameThroughSimbody(
	SimTK::Compound::AtomIndex aIx,
	SimTK::DuMMForceFieldSubsystem& dumm,
	SimTK::SimbodyMatterSubsystem& matter,
	const SimTK::State& someState)
{
	const SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndexThroughDumm(aIx, dumm);
	const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(mbx);

	const Transform&    X_GB = mobod.getBodyTransform(someState);
	const Rotation&     R_GB = X_GB.R();
	const Vec3&         p_GB = X_GB.p();

	SimTK::Vec3 station = getAtomLocationInMobilizedBodyFrameThroughDumm(aIx, dumm);

	const Vec3 p_BS_G = R_GB * station;
	return p_GB + p_BS_G;

}

/** Print maps **/
void Topology::printMaps()
{
/*
	std::cout << "Topology " << name << " maps " << std::endl;
	std::cout << "mbx2aIx:" << std::endl;
	map<SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex>::const_iterator mbx2aIxIt;
	for(mbx2aIxIt = mbx2aIx.begin();
	   mbx2aIxIt != mbx2aIx.end(); ++mbx2aIxIt)
	{
		std::cout << "mbx " << mbx2aIxIt->first
			<< " atomIndex " << mbx2aIxIt->second << std::endl;
	}
*/
	std::cout << "Topology map aIx2mbx:" << std::endl;
	map<SimTK::Compound::AtomIndex, std::vector<SimTK::MobilizedBodyIndex>>::const_iterator aIx2mbxIt;
	for(aIx2mbxIt = aIx2mbx.begin();
	   aIx2mbxIt != aIx2mbx.end(); ++aIx2mbxIt)
	{
		std::cout << "atomIndex " << aIx2mbxIt->first << " mbxs:";
		for (auto val : (aIx2mbxIt->second)){
			std::cout << " " << val ;
		}
		std::cout << std::endl << std::flush;
	}

	std::cout << "Topology map CompoundAtomIx2GmolAtomIx:" << std::endl;
	std::map< SimTK::Compound::AtomIndex, int >::const_iterator aIx2gmolaIxIt;
	for(aIx2gmolaIxIt = CompoundAtomIx2GmolAtomIx.begin();
	   aIx2gmolaIxIt != CompoundAtomIx2GmolAtomIx.end(); ++aIx2gmolaIxIt)
	{
		std::cout << "atomIndex " << aIx2gmolaIxIt->first
			<< " gmolaIx " << aIx2gmolaIxIt->second
			<< std::endl << std::flush;
	}
/*
	map<SimTK::Compound::BondIndex, int>::const_iterator bondIx2GmolBondIt;
	for(bondIx2GmolBondIt = bondIx2GmolBond.begin();
	   bondIx2GmolBondIt != bondIx2GmolBond.end(); ++bondIx2GmolBondIt)
	{
		std::cout << "Compound bondIndex " << bondIx2GmolBondIt->first
			<< " bBond index " << bondIx2GmolBondIt->second << std::endl;
	}

	map<int, SimTK::Compound::BondIndex>::const_iterator GmolBond2bondIxIt;
	for(GmolBond2bondIxIt = GmolBond2bondIx.begin();
		GmolBond2bondIxIt != GmolBond2bondIx.end(); ++GmolBond2bondIxIt)
	{
		std::cout << "bBond index " << GmolBond2bondIxIt->first
			<< " Compound index " << GmolBond2bondIxIt->second << std::endl;
	}
*/
}

/** Write a pdb with bAtomList coordinates and inNames **/
void Topology::writeAtomListPdb(std::string dirname,
	std::string prefix,
	std::string sufix,
	int maxNofDigits,
	int index) const
{
	// Using floor here is no buneo because the index can be zero
	std::string zeros("");
	int nofDigits = static_cast<int>(std::to_string(index).size());
	if(maxNofDigits > nofDigits){
		zeros = std::string(maxNofDigits - nofDigits, '0');
	}

	std::stringstream sstream;
	sstream << dirname << "/"
		<< prefix << zeros << std::to_string(index) << sufix;
	string ofilename = sstream.str();

	FILE *oF = fopen (ofilename.c_str(),"w");
	if (oF) {
		// Pdb lines
		for(int i = 0; i < getNumAtoms(); i++){
			fprintf(oF, "%-6s%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f  %4.2f%6.2f          %2s\n"
				, "ATOM"                 // record
				, i                      // index
				, subAtomList[i].getInName().c_str()  // name
				, "UNK"                  // residue name
				, 'A'                    // chain
				, 1                      // residue index
				, 10.0*subAtomList[i].getX()    // x in A
				, 10.0*subAtomList[i].getY()    // y in A
				, 10.0*subAtomList[i].getZ()    // z in A
				, 1.0                    // occupancy
				, 0.0                    // beta factor
				, "  ");                 // element
		}

		fclose(oF);
		//std::cout << "\tTopology written to '" << ofilename << "'\n";
	} else {
		std::cout << "FAILED TO OPEN '" << ofilename << "' TO WRIE!\n";
	}
}

/** Get bAtomList coordinates **/
void Topology::getCoordinates(
		std::vector<SimTK::Real>& Xs,
		std::vector<SimTK::Real>& Ys,
		std::vector<SimTK::Real>& Zs)
{
	assert(Xs.size() == static_cast<size_t>(getNumAtoms()));
	assert(Ys.size() == static_cast<size_t>(getNumAtoms()));
	assert(Zs.size() == static_cast<size_t>(getNumAtoms()));
	for(int ix = 0; ix < getNumAtoms(); ++ix){
		Xs[ix] = subAtomList[ix].getX();
		Ys[ix] = subAtomList[ix].getY();
		Zs[ix] = subAtomList[ix].getZ();
	}
}

/*!
 * <!--  -->
*/
void Topology::setSubAtomList(
	std::vector<bSpecificAtom>::iterator beginArg,
	std::vector<bSpecificAtom>::iterator endArg,
	ELEMENT_CACHE& elementCacheArg)
{		
	subAtomList.set_view( beginArg, endArg );
	rootAtomIx = bSpecificAtomRootIndex - beginArg->getNumber();

	natoms = subAtomList.size();

	atomFrameCache.resize(natoms);

	elementCache = elementCacheArg;
}

/*!
 * <!--  -->
*/
void Topology::setSubBondList(
	std::vector<bBond>::iterator beginArg,
	std::vector<bBond>::iterator endArg)
{		

	subBondList.set_view( beginArg, endArg );

	nbonds = (subBondList).size();
}

/** Get own CompoundIndex in CompoundSystem **/
const CompoundSystem::CompoundIndex &Topology::getCompoundIndex() const
{
	return compoundIndex;
}

/** Set the compoundIndex which is the position in the vector of Compounds
 * of the CompoundSystem **/
void Topology::setCompoundIndex(
	const CompoundSystem::CompoundIndex &compoundIndex)
{
	//Topology::compoundIndex = compoundIndex;
	this->compoundIndex = compoundIndex;
}

/** Get the neighbor atom bonded to aIx atom in the parent mobilized body.
TODO: No chemical parent for satelite atoms or first atom. **/
SimTK::Compound::AtomIndex
Topology::getChemicalParent_IfIAmRoot(
	SimTK::SimbodyMatterSubsystem *matter,
	//std::unique_ptr<SimTK::SimbodyMatterSubsystem> matter,
	SimTK::Compound::AtomIndex aIx,
	SimTK::DuMMForceFieldSubsystem& dumm)
{

	SimTK::Compound::AtomIndex chemParentAIx;
	int gmolAtomIndex = -111111;

	if(getAtomLocationInMobilizedBodyFrameThroughDumm(aIx, dumm) == 0){

		// Get body, parentBody, parentAtom
		SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndexThroughDumm(aIx, dumm);
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
		SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

		if(parentMobod.getMobilizedBodyIndex() != 0){ // parent not Ground
			// Find the true bSpecificAtom (CHEMICAL) parent
			bSpecificAtom *originSpecAtom = updAtomByAtomIx(aIx); //TODO: optimize

			// TODO: Check is neighbors and bondsInvolved are redundant
			// Loop through neighbor atoms (bSpecificAtom)
			for(auto k : originSpecAtom->neighborsIndex) {

				// Loop through bonds that this atom is involved in (bBond);
				for (auto bondIndex : originSpecAtom->bondsInvolvedIndex) {

					// Check if this neighbor is involved in this bond
					if( bonds[bondIndex].isThisMe(originSpecAtom->getNumber(), subAtomList[k].getNumber()
						) ){

						Compound::AtomIndex candidateChemParentAIx = subAtomList[k].getCompoundAtomIndex();

						// Check if neighbor atom's mobod is a parent mobod
						if(getAtomMobilizedBodyIndexThroughDumm(candidateChemParentAIx, dumm) == parentMbx){

							if(!bonds[bondIndex].isRingClosing()){ // No ring atoms are allowed
								chemParentAIx = candidateChemParentAIx;
								gmolAtomIndex = subAtomList[k].getNumber();
								break;
							}
						}
					}
				}
			}
		}
		
	}else{
		std::cout << "Warning: requiring chemical parent for non-root atom\n";
		bSpecificAtom *originSpecAtom = updAtomByAtomIx(aIx); //TODO: optimize
		for(auto k : originSpecAtom->neighborsIndex) {
			Compound::AtomIndex candidateChemParentAIx = subAtomList[k].getCompoundAtomIndex();
			if(getAtomLocationInMobilizedBodyFrameThroughDumm(candidateChemParentAIx, dumm) == 0){ // atom is at body's origin // DANGER
				chemParentAIx = candidateChemParentAIx;
				gmolAtomIndex = subAtomList[k].getNumber();
				break;
			}
		}
	}

	return chemParentAIx;
}



