#include "TrivalentAtomTetra.hpp"

/*
#ifndef DEBUG_LEVEL01
#define DEBUG_LEVEL01
#endif

#ifndef DEBUG_LEVEL02
#define DEBUG_LEVEL02
#endif
*/

//using namespace std;
using namespace SimTK;

/****
 *  Trivalent Atom Class with tetrahedral geometry.
 *  Bond centers are named "bond1", "bond2", and "bond3"
 ****/
TrivalentAtomTetra::TrivalentAtomTetra(
        const Compound::AtomName& atomName,   ///< name for new atom
        const Element& element   ///< chemical element for new atom
) : Compound::SingleAtom(atomName, element)
{
    static const Angle TetrahedralAngle = 109.47 * Deg2Rad;
  
    // BondCenter1 dihedral will be relative to BondCenter2
    addFirstBondCenter( "bond1", atomName );
  
    // bond centers 2 and 3 dihedrals relative to bond center 1
    addSecondBondCenter( "bond2", atomName,  TetrahedralAngle);
    //addLeftHandedBondCenter( "bond2", atomName, TetrahedralAngle, TetrahedralAngle );
    addRightHandedBondCenter( "bond3", atomName, TetrahedralAngle, TetrahedralAngle );
  
    // Choice of inboard bond may differ from bond priority - user may change this
    setInboardBondCenter("bond1");
    setCompoundName("TrivalentAtomTetra"); // overridden
}

