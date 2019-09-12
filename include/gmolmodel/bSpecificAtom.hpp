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
    int nbonds;
    int freebonds;
    char name[5];
    char inName[5];
    int number;
    char elem;
    int atomicNumber;
    SimTK::Real mass;
    SimTK::Real vdwRadius;
    SimTK::Real LJWellDepth; // Lennard-Jones well depth
    double charge;
    char fftype[20];
    SimTK::DuMM::AtomClassIndex atomClassIndex;



    float x;
    float y;
    float z;
    //char biotype[20];
    std::string biotype; // NEW
    SimTK::BiotypeIndex biotypeIndex;

    // We need this to be a pointer because we don't know it's size at
    // initializetion. It could be Univalent, Bivalent, TrivalentAtom.
    // The size of this variable will be known  after Topology loads info
    // from an inputReader
    SimTK::Compound::SingleAtom *bAtomType;
    SimTK::Compound::AtomIndex atomIndex;
    int mobile;
    int visited;
    
    // Graph useful vars
    std::vector<bSpecificAtom *> neighbors;
    std::vector<bBond *> bondsInvolved;

    // Residue and chain
    std::string residueName;
    long int residueIndex;
    std::string chain;
    int moleculeIndex;
  
    bSpecificAtom();
    ~bSpecificAtom();
    void Print(void);
    void Zero(void);

    // Interface
    int getNBonds(void);
    int getFreebonds(void);
    std::string getName(void);
    std::string getInName(void);
    int getNumber(void);
    char getElem(void);
    SimTK::mdunits::Mass getMass(void);
    void setMass(SimTK::Real);
    SimTK::Real getX(void) const;
    SimTK::Real getY(void) const;
    SimTK::Real getZ(void) const;
    std::string getFftype(void);
    SimTK::DuMM::AtomClassIndex getAtomClassIndex(void);
    void setAtomClassIndex(SimTK::DuMM::AtomClassIndex);
    std::string getBiotype(void);
    SimTK::Compound::SingleAtom * getBAtomType(void);
    SimTK::Compound::AtomIndex getCompoundAtomIndex(void);
    SimTK::Real getCharge(void);
    int getIsMobile(void);
    int getIsVisited(void);

    int getAtomicNumber(void);
    void setAtomicNumber(int);
    SimTK::Real getVdwRadius(void);
    void setVdwRadius(SimTK::Real);
    SimTK::Real getLJWellDepth(void);
    void setLJWellDepth(SimTK::Real);

    const SimTK::DuMM::ChargedAtomTypeIndex getChargedAtomTypeIndex() const;
    void setChargedAtomTypeIndex(const SimTK::DuMM::ChargedAtomTypeIndex);

    SimTK::BiotypeIndex getBiotypeIndex(void);
    void setBiotypeIndex(SimTK::BiotypeIndex);

    void setNbonds(int);
    void setFreebonds(int);
    void setName(std::string);
    void setInName(std::string);
    void setNumber(int);
    void setElem(char);


    void setX(SimTK::Real);
    void setY(SimTK::Real);
    void setZ(SimTK::Real);
    void setFftype(std::string);
    void setBiotype(std::string);
    void setBiotype(const char *);
    void setBAtomType(SimTK::Compound::SingleAtom *);
    void setCompoundAtomIndex(SimTK::Compound::AtomIndex);
    void setCharge(SimTK::Real);
    void setIsMobile(int);
    void setVisited(int);

    void addNeighbor(bSpecificAtom *);
    void addBond(bBond *);

private:
    SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;

};

// Update Molmodel MolAtom dest with Gmolmodel bSpecificAtom src values
int bAtomAssign(MolAtom *dest, const bSpecificAtom *src);

#endif  //__BSPECIFICATOM__


