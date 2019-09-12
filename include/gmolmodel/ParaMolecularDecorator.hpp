#ifndef __PARAMOLECULARDECORATOR_HPP__
#define __PARAMOLECULARDECORATOR_HPP__

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

#include "Robo.hpp"
#include "Topology.hpp"

using namespace SimTK;

class ParaMolecularDecorator : public DecorationGenerator {
public:
    ParaMolecularDecorator(SimTK::CompoundSystem *argCompoundSystem,
        SimTK::SimbodyMatterSubsystem *argMatter,
        //Topology *argResidue, // RE
        SimTK::DuMMForceFieldSubsystem *argDumm,
        SimTK::GeneralForceSubsystem *argForces);

    void AddMolecule(Topology *argMolecule);

    void loadPoint(const Vec3 point);

    void loadLine(const Vec3 p1, const Vec3 p2);

    void clearPoints(void);

    void clearLines(void);

    void generateDecorations(const State& state,
        Array_<DecorativeGeometry>& geometry);

    ~ParaMolecularDecorator(void);

    // Gmolmodel specific
    void setAtomTargets(std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>> residueAtomLocations);

private:
    SimTK::CompoundSystem *compoundSystem;
    SimTK::SimbodyMatterSubsystem *matter;
    SimTK::DuMMForceFieldSubsystem *dumm;
    SimTK::GeneralForceSubsystem *forces;

    // Topology *molecule; // RE
    std::vector<Topology *> molecules; // NEW

    Array_< Vec3 >  points;
    Array_< std::pair< Vec3, Vec3 > > lines;

    // Gmolmodel specific
    std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;

};

#endif // __PARAMOLECULARDECORATOR_HPP__
