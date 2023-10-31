#ifndef __BSPECIFICATOM__
#define __BSPECIFICATOM__

/* -------------------------------------------------------------------------- *
 *                                gMolModel                                   *
 * -------------------------------------------------------------------------- *
 *  This is an extension of SimTK Molmodel package.                           *
 * -------------------------------------------------------------------------- */
/** @file
 * This defines the bMoleculeReader class and additional heloer classes
 **/

#include "bgeneral.hpp"
#include "Robo.hpp"
#include "Simbody.h"
#include "Molmodel.h"

#include "bBond.hpp"
//class bBond;

//==============================================================================
//                           CLASS SpecificAtom
//==============================================================================
/** 
 * gMolmodel Specific Atom Type Class.
 * This incorporates additional Amber forcefield data.
**/
class bSpecificAtom{ /*Information like in sdf*/
public:
    void setAtomCompoundType(const SimTK::Compound::AtomName &atomName,
        int atomicNumber,
        SimTK::Element::Name elementName,
        SimTK::Element::Symbol elementSymbol,
        SimTK::mdunits::Mass atomicMass);

    void destroy();

    std::string biotype; // moleculeName + force field type + fftype (MOL0AAAAoh, MOL0AAAOho)
    std::string residueName; // Residue and chain
    std::string elem; // H, C, N, O, S, P, F, Cl, Br, I etc
    std::string fftype; // atom name from the force field ("O1", "C1", "C2", "H1", "H10")
    std::string name; // combination of four letters that depends on the amber index
    std::string inName; // original name in the sdf file ("O1", "C1", "C2", "H1", "H10")

    // Graph useful vars
    std::vector<bSpecificAtom *> neighbors;
    std::vector<bBond *> bondsInvolved;

    SimTK::Real charge = SimTK::NaN;
    int residueIndex = std::numeric_limits<int>::min();

    SimTK::Real mass = SimTK::NaN;
    SimTK::Real vdwRadius = SimTK::NaN;
    SimTK::Real LJWellDepth = SimTK::NaN;

    SimTK::Real x = SimTK::NaN;
    SimTK::Real y = SimTK::NaN;
    SimTK::Real z = SimTK::NaN;

    int nbonds = std::numeric_limits<int>::min();
    int freebonds = std::numeric_limits<int>::min();
    int number = std::numeric_limits<int>::min(); // amber index
    int atomicNumber = std::numeric_limits<int>::min();
    int visited = std::numeric_limits<int>::min();
    int moleculeIndex = std::numeric_limits<int>::min();
    
    SimTK::DuMM::AtomClassIndex dummAtomClassIndex;
    SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;
    SimTK::BiotypeIndex biotypeIndex;
    SimTK::Compound::AtomIndex compoundAtomIndex;
  
public:
    void Print(int whichWorld);

    // Interface
    int getNBonds() const;
    int getFreebonds() const;
    std::string getName() const;
    std::string getInName() const;
    int getNumber() const;
    std::string getElem() const;
    SimTK::mdunits::Mass getMass() const;
    void setMass(SimTK::Real);
    SimTK::Real getX() const;
    SimTK::Real getY() const;
    SimTK::Real getZ() const;
    std::string getBiotype() const;

    const SimTK::Compound::SingleAtom& getSingleAtom() const;
    void setCompoundName(const SimTK::Compound::Name& name);
    void addFirstBondCenter(const SimTK::Compound::BondCenterName& centerName,
        const SimTK::Compound::AtomPathName& atomName);
    void addFirstTwoBondCenters(const SimTK::Compound::BondCenterName& centerName1,
        const SimTK::Compound::BondCenterName& centerName2,
        const SimTK::Compound::AtomPathName& atomName,
        SimTK::UnitVec3 dir1,
        SimTK::UnitVec3 dir2);
    void addLeftHandedBondCenter(const SimTK::Compound::BondCenterName& centerName,
        const SimTK::Compound::AtomName& atomName,
        SimTK::Angle bondAngle1,
        SimTK::Angle bondAngle2);
    void addRightHandedBondCenter(const SimTK::Compound::BondCenterName& centerName,
        const SimTK::Compound::AtomName& atomName,
        SimTK::Angle bondAngle1,
        SimTK::Angle bondAngle2);
    void setInboardBondCenter(const SimTK::Compound::BondCenterName& centerName);
    void setDefaultInboardBondLength(SimTK::mdunits::Length length);

    SimTK::Compound::AtomIndex getCompoundAtomIndex() const;
    SimTK::Real getCharge() const;
    int getIsVisited() const;
    std::string getFftype() const;

    SimTK::DuMM::AtomClassIndex getDummAtomClassIndex() const;
    void setDummAtomClassIndex(SimTK::DuMM::AtomClassIndex);

    int getAtomicNumber() const;
    void setAtomicNumber(int);

    SimTK::Real getVdwRadius() const;
    void setVdwRadius(SimTK::Real);

    SimTK::Real getLJWellDepth() const;
    void setLJWellDepth(SimTK::Real);

    SimTK::DuMM::ChargedAtomTypeIndex getChargedAtomTypeIndex() const;
    void setChargedAtomTypeIndex(const SimTK::DuMM::ChargedAtomTypeIndex);

    SimTK::BiotypeIndex getBiotypeIndex() const;
    void setBiotypeIndex(SimTK::BiotypeIndex);

    void setNbonds(int);
    void setFreebonds(int);
    void generateName(int index);
    void setInName(const std::string& name);
    void setNumber(int);
    void setElem(const std::string& elem);

    void setX(SimTK::Real);
    void setY(SimTK::Real);
    void setZ(SimTK::Real);
    void setFfType(const std::string&);
    void setBiotype(const std::string&);
    void setCompoundAtomIndex(SimTK::Compound::AtomIndex);
    void setCharge(SimTK::Real);
    void setVisited(int);

    void addNeighbor(bSpecificAtom *);
    void addBond(bBond *);

public:
    // wasted 6 hours trying to make this unique_ptr or allocated on the stack
    SimTK::Compound::SingleAtom* compoundSingleAtom = nullptr;
};

// Update Molmodel MolAtom dest with Gmolmodel bSpecificAtom src values
int bAtomAssign(MolAtom *dest, const bSpecificAtom *src);

#endif  //__BSPECIFICATOM__


