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
    std::string biotype; // NEW
    std::string residueName; // Residue and chain
    std::string chain; // Residue and chain

    // Graph useful vars
    std::vector<bSpecificAtom *> neighbors;
    std::vector<bBond *> bondsInvolved;

    //char biotype[20];

    // used to be initialized with bZeroCharArray
    // as per c++11
    // it should be initialized to default
    char fftype[20] {};

    // We need this to be a pointer because we don't know it's size at
    // initializetion. It could be Univalent, Bivalent, TrivalentAtom.
    // The size of this variable will be known  after Topology loads info
    // from an inputReader
    SimTK::Compound::SingleAtom *compoundSingleAtom = nullptr;

    double charge = 0.0;
    long int residueIndex = std::numeric_limits<long int>::min(); // Residue and chain

    SimTK::Real mass = std::numeric_limits<SimTK::Real>::min();
    SimTK::Real vdwRadius = std::numeric_limits<SimTK::Real>::min();
    SimTK::Real LJWellDepth = std::numeric_limits<SimTK::Real>::min(); // Lennard-Jones well depth

    SimTK::Real x = std::numeric_limits<SimTK::Real>::min();
    SimTK::Real y = std::numeric_limits<SimTK::Real>::min();
    SimTK::Real z = std::numeric_limits<SimTK::Real>::min();

    char name[5] {};
    char inName[5] {};

    int nbonds = 0;
    int freebonds = 0;
    int number = 0;
    int atomicNumber = 0;
    int mobile = 0;
    int visited = 0;
    int moleculeIndex = std::numeric_limits<int>::min(); // Residue and chain
    
    SimTK::DuMM::AtomClassIndex atomClassIndex;
    SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;
    SimTK::BiotypeIndex biotypeIndex;
    SimTK::Compound::AtomIndex compoundAtomIndex;

    std::string elem;
  
public:
    void Print(int whichWorld);
    void Zero();

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
    SimTK::Compound::SingleAtom * getBAtomType() const;
    SimTK::Compound::AtomIndex getCompoundAtomIndex() const;
    SimTK::Real getCharge() const;
    int getIsMobile() const;
    int getIsVisited() const;
    std::string getFftype() const;

    SimTK::DuMM::AtomClassIndex getAtomClassIndex() const;
    void setAtomClassIndex(SimTK::DuMM::AtomClassIndex);

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
    void setName(std::string);
    void setInName(std::string);
    void setNumber(int);
    void setElem(std::string);

    void setX(SimTK::Real);
    void setY(SimTK::Real);
    void setZ(SimTK::Real);
    void setFfType(std::string);
    void setBiotype(std::string);
    void setBiotype(const char *);
    void setBAtomType(SimTK::Compound::SingleAtom *);
    void setCompoundAtomIndex(SimTK::Compound::AtomIndex);
    void setCharge(SimTK::Real);
    void setIsMobile(int);
    void setVisited(int);

    void addNeighbor(bSpecificAtom *);
    void addBond(bBond *);
};

// Update Molmodel MolAtom dest with Gmolmodel bSpecificAtom src values
int bAtomAssign(MolAtom *dest, const bSpecificAtom *src);

#endif  //__BSPECIFICATOM__


