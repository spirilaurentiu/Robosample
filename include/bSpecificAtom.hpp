#pragma once

/* -------------------------------------------------------------------------- *
 *                                gMolModel                                   *
 * -------------------------------------------------------------------------- *
 *  This is an extension of SimTK Molmodel package.                           *
 * -------------------------------------------------------------------------- */
/** @file
 * This defines the SpecificAtom class
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
 * gMolmodel Specific Atom Class which incorporates everything below.
 * This incorporates additional Amber forcefield data.
**/
class bSpecificAtom{ /*Information like in sdf*/

public:
    ///////////////////////
    //     INTERFACE     //
    ///////////////////////

    //------------------------------------------------------------------------------
    /** @name Interface - Interface. **/

    /**@{**/

    /** TODO comment. **/
    void setAtomCompound(const SimTK::Element &element);

    /** TODO comment. **/
    void deleteAtomCompound();

public:
    void Print(int whichWorld) const;

    // Interface
    int getNBonds() const;
    int getFreebonds() const;
    std::string getName() const;
    void setName(const std::string& name);

    std::string getInName() const;
    int getNumber() const;
    std::string getElem() const;
    SimTK::mdunits::Mass getMass() const;
    void setMass(SimTK::Real);
    SimTK::Real getX() const;
    SimTK::Real getY() const;
    SimTK::Real getZ() const;
    const double* getCartesians() const;

    std::string getBiotype() const;

    const SimTK::Compound::SingleAtom& getSingleAtom() const;
    void setCompoundName(const SimTK::Compound::Name& name);

    SimTK::Compound::AtomIndex getCompoundAtomIndex() const;
    SimTK::Real getCharge() const;
    bool wasVisited() const;
    int getVisited() const;

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

    void setNbonds(std::size_t);
    void setFreebonds(std::size_t);
    void incrFreebonds(void);
    void decrFreebonds(void);

    void generateName(int index);
    void setInName(const std::string& name);
    void setNumber(int);
    void setElem(const std::string& elem);

    void setX(SimTK::Real);
    void setY(SimTK::Real);
    void setZ(SimTK::Real);

    void setCartesians(SimTK::Real, SimTK::Real, SimTK::Real);
    // void setCartesians(double *);

    void setCartesians(SimTK::Vec3);
    void setCartesians(OpenMM::Vec3 );

    void setFfType(const std::string&);
    void setBiotype(const std::string&);
    void setCompoundAtomIndex(SimTK::Compound::AtomIndex);
    void setCharge(SimTK::Real);
    void setVisited(bool);

    // void addNeighbor(bSpecificAtom *);
    // void addBond(bBond *);

    void addNeighborIndex(int index);
    void addBondIndex(int index);
    std::size_t getNumBonds() const;

    // Getter and setter for the residueName property
    std::string getResidueName() const;
    void setResidueName(const std::string& residueName);

    const int getMoleculeIndex() const;
    void setMoleculeIndex(int);

    void setParentNumber(int n);
    int getParentNumber() const;

    void setIsRoot(bool argIsRoot){this->isRoot = argIsRoot;}
    const bool getIsRoot() const {return this->isRoot;}

    void setDuMMAtomIndex(SimTK::DuMM::AtomIndex dAIx_){this->dAIx = dAIx_;}
    const SimTK::DuMM::AtomIndex& getDuMMAtomIndex(void){return this->dAIx;}
    SimTK::DuMM::AtomIndex updDuMMAtomIndex(void) const {return this->dAIx;}


    // End of interface
    /**@}**/

    // Graph useful vars
    // std::vector<bSpecificAtom *> neighbors;
    // std::vector<bBond *> bondsInvolved;

    std::vector<int> neighborsIndex;
    std::vector<int> bondsInvolvedIndex;

    SimTK::Real charge = SimTK::NaN;
    int residueIndex = std::numeric_limits<int>::min();

    SimTK::Real mass = SimTK::NaN;
    SimTK::Real vdwRadius = SimTK::NaN;
    SimTK::Real LJWellDepth = SimTK::NaN;

    SimTK::Real x = SimTK::NaN;
    SimTK::Real y = SimTK::NaN;
    SimTK::Real z = SimTK::NaN;

    double Cartesians[3];
    OpenMM::Vec3 OpenMMCoords;

private:
    //int nbonds = std::numeric_limits<int>::min(); // SP_OLD
    std::size_t nbonds = 0; // SP_NEW
    std::size_t freebonds = std::numeric_limits<std::size_t>::min();
    int number = std::numeric_limits<int>::min(); // amber index
    int parentNumber = std::numeric_limits<int>::min(); // amber index of parent atom
    int atomicNumber = std::numeric_limits<int>::min(); // atomic number
    //int visited = std::numeric_limits<int>::min(); // SP_OLD
    int visited = 0; // SP_NEW
    int moleculeIndex = -111111;
    SimTK::DuMM::AtomClassIndex dummAtomClassIndex;
    SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;
    SimTK::BiotypeIndex biotypeIndex;
    SimTK::DuMM::AtomIndex dAIx;

public:
    
    SimTK::Compound::AtomIndex compoundAtomIndex;
  

public:
    // wasted 6 hours trying to make this unique_ptr or allocated on the stack
    // after more hours wasted, i read this article: https://www.cppstories.com/2014/05/vector-of-objects-vs-vector-of-pointers/
    // as expected, it should be faster to have this allocated on the stack
    // however, taking into account that this is called only a few times, it is not worth the trouble
    // there are two use cases for this: addition of bond centers and building the molecule graph
    // in the first case, the calls are sequential, so we waste some time
    // in the second case, the calls are random and the article shows that it is actually better to have the objects allocated on the heap
    SimTK::Compound::SingleAtom* compoundSingleAtom = nullptr;

private:

    std::string biotype; // moleculeName + force field type + fftype (MOL0AAAAoh, MOL0AAAOho)
    std::string residueName; // Residue and chain
    std::string elem; // H, C, N, O, S, P, F, Cl, Br, I etc
    std::string fftype; // atom name from the force field ("O1", "C1", "C2", "H1", "H10")
    std::string name; // combination of four letters that depends on the amber index
    std::string inName; // original name in the sdf file ("O1", "C1", "C2", "H1", "H10")

    bool isRoot = false; // specifies if this atom is a base atom


};

// // Update Molmodel MolAtom dest with Gmolmodel bSpecificAtom src values
// int bAtomAssign(MolAtom *dest, const bSpecificAtom *src);
