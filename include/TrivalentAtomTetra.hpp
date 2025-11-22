#ifndef __TRIVALENTATOMTETRA__
#define __TRIVALENTATOMTETRA__

/* -------------------------------------------------------------------------- *
 *                                gMolModel                                   *
 * -------------------------------------------------------------------------- *
 *  This is an extension of SimTK Molmodel package.                           *
 * -------------------------------------------------------------------------- */
/** @file
 * This defines the TrivalentAtomTetra class
 **/

#include "Robo.hpp"
#include "Simbody.h"
#include "Molmodel.h"

/*
#ifndef DEBUG_LEVEL01
#define DEBUG_LEVEL01
#endif

#ifndef DEBUG_LEVEL02
#define DEBUG_LEVEL02
#endif
*/

//==============================================================================
//                           CLASS TrivalentAtomTetra
//==============================================================================
/** 
 * Trivalent Atom Class with tetrahedral geometry (ex. positive N)
 * Bond centers are named "bond1", "bond2", and "bond3"
**/

class  TrivalentAtomTetra : public SimTK::Compound::SingleAtom {
 public:
  TrivalentAtomTetra(
    const SimTK::Compound::AtomName& atomName,   ///< name for new atom
    const SimTK::Element& element              /// element for new atom
  );
};

#endif  //__TRIVALENTATOMTETRA__

