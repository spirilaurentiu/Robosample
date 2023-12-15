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

//==============================================================================
//                   CLASS ForceArrowGenerator
//==============================================================================
/* 
*	Added by Teodor
*	Shamelessly copied from the ExampleContactPlayground CPP file
*	(from Simbody) and adapted to work in Robosample
*/

class ForceArrowGenerator : public DecorationGenerator {
public:
    ForceArrowGenerator(const MultibodySystem& system,
                        const CompliantContactSubsystem& complCont) 
    :   m_system(system), m_compliant(complCont) {}

    virtual void generateDecorations(const State& state, Array_<DecorativeGeometry>& geometry) override {
		static const Real ForceScale = .25;
        const Vec3 frcColors[] = {Red,Orange,Cyan};
        const Vec3 momColors[] = {Blue,Green,Purple};
        m_system.realize(state, Stage::Velocity);

        const int ncont = m_compliant.getNumContactForces(state);
        for (int i=0; i < ncont; ++i) {
            const ContactForce& force = m_compliant.getContactForce(state,i);
            const ContactId     id    = force.getContactId();
            const Vec3& frc = force.getForceOnSurface2()[1];
            const Vec3& mom = force.getForceOnSurface2()[0];
            Real  frcMag = frc.norm(), momMag=mom.norm();
            int frcThickness = 1, momThickness = 1;
            Real frcScale = ForceScale, momScale = ForceScale;
            while (frcMag > 10)
                frcThickness++, frcScale /= 10, frcMag /= 10;
            while (momMag > 10)
                momThickness++, momScale /= 10, momMag /= 10;
            DecorativeLine frcLine(force.getContactPoint(),
                force.getContactPoint() + frcScale*frc);
            DecorativeLine momLine(force.getContactPoint(),
                force.getContactPoint() + momScale*mom);
            frcLine.setColor(frcColors[id%3]);
            momLine.setColor(momColors[id%3]);
            frcLine.setLineThickness(2*frcThickness);
            momLine.setLineThickness(2*momThickness);
            geometry.push_back(frcLine);
            geometry.push_back(momLine);

            ContactPatch patch;
            const bool found = m_compliant.calcContactPatchDetailsById(state,id,patch);
            //cout << "patch for id" << id << " found=" << found << endl;
            //cout << "resultant=" << patch.getContactForce() << endl;
            //cout << "num details=" << patch.getNumDetails() << endl;
            for (int i=0; i < patch.getNumDetails(); ++i) {
                const ContactDetail& detail = patch.getContactDetail(i);
                const Real peakPressure = detail.getPeakPressure();
                // Make a black line from the element's contact point in the normal
                // direction, with length proportional to log(peak pressure)
                // on that element. 
                DecorativeLine normal(detail.getContactPoint(),
                    detail.getContactPoint()+ std::log10(peakPressure)
                                                * detail.getContactNormal());
                normal.setColor(Black);
                geometry.push_back(normal);
                // Make a red line that extends from the contact
                // point in the direction of the slip velocity, of length 3*slipvel.
                DecorativeLine slip(detail.getContactPoint(),
                    detail.getContactPoint()+3*detail.getSlipVelocity());
                slip.setColor(Red);
                geometry.push_back(slip);
            }
        }
    }
private:
    const MultibodySystem&              m_system;
    const CompliantContactSubsystem&    m_compliant;
};


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
