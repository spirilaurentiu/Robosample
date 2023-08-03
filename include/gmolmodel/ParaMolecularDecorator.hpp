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
		SimTK::DuMMForceFieldSubsystem *argDumm,
		SimTK::GeneralForceSubsystem *argForces
	);

	void AddMolecule(Topology *argMolecule);

	void loadPoint(const Vec3 point);

	void loadLine(const Vec3 p1, const Vec3 p2);

	void loadArrow(const Vec3 p1, const Vec3 p2);

	void updateArrow(int which, const Vec3 p1, const Vec3 p2);

	void clearPoints(void);

	void clearLines(void);

	// Draw a frame
	void drawFrame(
		Array_<DecorativeGeometry>& geometry, SimTK::Transform G_X_F,
		SimTK::Real scaleFactor, SimTK::Real lineThickness, SimTK::Vec3 color,
		std::string text, SimTK::Real textScaleFactor, SimTK::Vec3 textColor,
		SimTK::Vec3 textOffset
	);

	// Draw a line
	void drawLine(Array_<DecorativeGeometry>& geometry,
		SimTK::Transform G_X_B, SimTK::Transform G_X_M,
		int numOfDofs, SimTK::Real lineThickness
	);

	// Draw a line
	void drawArrow(Array_<DecorativeGeometry>& geometry,
		SimTK::Transform G_X_B, SimTK::Transform G_X_M,
		int numOfDofs, SimTK::Real lineThickness
	);

	void drawLoadedArrows(Array_<DecorativeGeometry>& geometry);

	// Draw DuMM based geometry
	void drawDummBasedGeometry(Array_<DecorativeGeometry>& geometry,
		const State& someState);

	void generateDecorations(const State& state,
		Array_<DecorativeGeometry>& geometry);

	//~ParaMolecularDecorator(void);

	// 
	void setAtomTargets(std::vector<std::pair<bSpecificAtom *,
		SimTK::Vec3>> residueAtomLocations);

	void updPCommVars(SimTK::Real argCommVar);
	void updFCommVars(SimTK::Real argCommVar);
	void updMCommVars(SimTK::Real argCommVar);
	void updBCommVars(SimTK::Real argCommVar);
	

private:
	SimTK::CompoundSystem *compoundSystem;
	SimTK::SimbodyMatterSubsystem *matter;
	SimTK::DuMMForceFieldSubsystem *dumm;
	SimTK::GeneralForceSubsystem *forces;

	std::vector<Topology *> molecules;

	Array_< Vec3 >  points;
	Array_< std::pair< Vec3, Vec3 > > lines;
	Array_< std::pair< Vec3, Vec3 > > arrows;

	// Gmolmodel specific
	std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;

	SimTK::Real PCommVar, FCommVar, MCommVar, BCommVar;

};

#endif // __PARAMOLECULARDECORATOR_HPP__
