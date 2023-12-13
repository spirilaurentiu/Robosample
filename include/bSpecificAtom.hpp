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
    bool setAtomCompoundType(const SimTK::Element &element);

    /** TODO comment. **/
    void destroy();

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

    void setNbonds(const int&);
    void setFreebonds(const int&);
    void incrFreebonds(void);
    void decrFreebonds(void);

    void setName(const std::string& name);
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
    void setVisited(const int&);

    void addNeighbor(bSpecificAtom *);
    void addBond(bBond *);

    // Getter and setter for the residueName property
    std::string getResidueName() const;
    void setResidueName(const std::string& residueName);

    // End of interface
    /**@}**/

    // Graph useful vars
    std::vector<bSpecificAtom *> neighbors;
    std::vector<bBond *> bondsInvolved;

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
    int nbonds = std::numeric_limits<int>::min();
    int freebonds = std::numeric_limits<int>::min();
    int number = std::numeric_limits<int>::min(); // amber index
    int atomicNumber = std::numeric_limits<int>::min();
    int visited = std::numeric_limits<int>::min();
    int moleculeIndex = std::numeric_limits<int>::min();
    SimTK::DuMM::AtomClassIndex dummAtomClassIndex;
    SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;
    SimTK::BiotypeIndex biotypeIndex;

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

};

// // Update Molmodel MolAtom dest with Gmolmodel bSpecificAtom src values
// int bAtomAssign(MolAtom *dest, const bSpecificAtom *src);
