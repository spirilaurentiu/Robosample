/**@file
Implementation of HMCSampler class. **/

#include "Robo.hpp"
#include "ParaMolecularDecorator.hpp"

using namespace SimTK;


//=============================================================================
//                   CLASS ForceArrowGenerator
//=============================================================================
/* 
*	Added by Teodor
*	Shamelessly copied from the ExampleContactPlayground CPP file
*	(from Simbody) and adapted to work in Robosample
*/

void ForceArrowGenerator::generateDecorations(const State& state, Array_<DecorativeGeometry>& geometry) {
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




ParaMolecularDecorator::ParaMolecularDecorator(
	SimTK::CompoundSystem *argCompoundSystem,
	SimTK::SimbodyMatterSubsystem *argMatter,
	//Topology *argResidue, RE
	SimTK::DuMMForceFieldSubsystem *argDumm,
	SimTK::GeneralForceSubsystem *argForces
)
{
	this->compoundSystem = argCompoundSystem;
	this->matter = argMatter;
	// this->molecule = argResidue; //RE
	this->dumm = argDumm;
	this->forces = argForces;

	MCommVar = SimTK::NaN;
	
}

void ParaMolecularDecorator::AddMolecule(Topology *argMolecule)
{
	molecules.push_back(argMolecule);
}

void ParaMolecularDecorator::loadPoint(const Vec3 point)
{
	points.push_back(point);
}

void ParaMolecularDecorator::loadLine(const Vec3 p1, const Vec3 p2)
{
	lines.push_back(std::pair<Vec3, Vec3>(p1, p2));
}

void ParaMolecularDecorator::loadArrow(const Vec3 p1, const Vec3 p2)
{
	arrows.push_back(std::pair<Vec3, Vec3>(p1, p2));
}

void ParaMolecularDecorator::updateArrow(int which, const Vec3 p1, const Vec3 p2)
{
	arrows[which] = std::pair<Vec3, Vec3>(p1, p2);
	std::cout << "ParaMolecularDecorator::updateArrow " << p1 << " " << p2 
		<< " " << arrows[which].first << " " << arrows[which].second << "\n";
}

void ParaMolecularDecorator::clearPoints(void)
{
	points.clear();
}

void ParaMolecularDecorator::clearLines(void)
{
	lines.clear();
}

// Gmolmodel specific
void ParaMolecularDecorator::setAtomTargets(
	std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>>
	residueAtomLocations)
{
	atomTargets.clear();
	for(unsigned int j = 0; j < residueAtomLocations.size(); j++){
		SimTK::Compound::AtomIndex atomIndex =
		((residueAtomLocations[j]).first)->compoundAtomIndex;
		SimTK::Vec3 location = ((residueAtomLocations[j]).second);
		atomTargets.insert(std::pair<
		SimTK::Compound::AtomIndex, SimTK::Vec3>(atomIndex, location));
	}
}

void ParaMolecularDecorator::updPCommVars(SimTK::Real argCommVar)
{
	this->PCommVar = argCommVar;
}


void ParaMolecularDecorator::updFCommVars(SimTK::Real argCommVar)
{
	this->FCommVar = argCommVar;
}
void ParaMolecularDecorator::updMCommVars(SimTK::Real argCommVar)
{
	this->MCommVar = argCommVar;
}

void ParaMolecularDecorator::updBCommVars(SimTK::Real argCommVar)
{
	this->BCommVar = argCommVar;
}


void ParaMolecularDecorator::drawLine(
	Array_<DecorativeGeometry>& geometry,
	SimTK::Transform G_X_B, SimTK::Transform G_X_M,
	int numOfDofs, SimTK::Real lineThickness
)
{
	DecorativeLine decorativeLineBM(G_X_B.p(), G_X_M.p());

	decorativeLineBM.setLineThickness(lineThickness);
	if(numOfDofs == 3){
		decorativeLineBM.setColor(SimTK::Vec3(SimTK::Orange));
	}else if(numOfDofs == 5){
		decorativeLineBM.setColor(SimTK::Vec3(SimTK::Magenta));
	}else if (numOfDofs == 2){
		decorativeLineBM.setColor(SimTK::Vec3(SimTK::Cyan));
	}else if (numOfDofs == 1){
		decorativeLineBM.setColor(SimTK::Vec3(SimTK::Green));
	}else{
		decorativeLineBM.setColor(SimTK::Vec3(SimTK::Red));
	}
	geometry.push_back(decorativeLineBM);
}

void ParaMolecularDecorator::drawArrow(
	Array_<DecorativeGeometry>& geometry,
	SimTK::Transform G_X_B, SimTK::Transform G_X_M,
	int numOfDofs, SimTK::Real lineThickness
)
{
	DecorativeArrow decorativeLineBM(G_X_B.p(), G_X_M.p());
	decorativeLineBM.setLineThickness(lineThickness);
	if(numOfDofs == 3){
		decorativeLineBM.setColor(SimTK::Vec3(SimTK::Orange));
	}else if(numOfDofs == 5){
		decorativeLineBM.setColor(SimTK::Vec3(SimTK::Magenta));
	}else if (numOfDofs == 2){
		decorativeLineBM.setColor(SimTK::Vec3(SimTK::Cyan));
	}else if (numOfDofs == 1){
		decorativeLineBM.setColor(SimTK::Vec3(SimTK::Green));
	}else{
		decorativeLineBM.setColor(SimTK::Vec3(SimTK::Red));
	}
	geometry.push_back(decorativeLineBM);
}

void ParaMolecularDecorator::drawLoadedArrows(Array_<DecorativeGeometry>& geometry)
{
	for(size_t tz = 0; tz < arrows.size(); tz++){
		DecorativeLine decorativeArrow(arrows[tz].first, arrows[tz].second);
		decorativeArrow.setLineThickness(4);
		decorativeArrow.setColor(SimTK::Vec3(SimTK::Orange));
		geometry.push_back(decorativeArrow);

		/* DecorativeSphere tip(0.2);
		tip.setColor(SimTK::Vec3(1, 0, 0));
		tip.setOpacity(0.5);
		tip.setTransform(SimTK::Transform(SimTK::Rotation(), arrows[tz].second));
		geometry.push_back(tip); */ 
	}

}

void ParaMolecularDecorator::drawFrame(
	Array_<DecorativeGeometry>& geometry, SimTK::Transform G_X_F,
	SimTK::Real scaleFactor, SimTK::Real lineThickness, SimTK::Vec3 color,
	std::string text, SimTK::Real textScaleFactor, SimTK::Vec3 textColor,
	SimTK::Vec3 textOffset
)
{
		DecorativeFrame decorativeFrameF;
		decorativeFrameF.setTransform(G_X_F);
		decorativeFrameF.setScaleFactors(SimTK::Vec3(
			scaleFactor, scaleFactor, scaleFactor));
		decorativeFrameF.setLineThickness(lineThickness);
		decorativeFrameF.setColor(color);
		geometry.push_back( decorativeFrameF );

		DecorativeText decorativeTextF(text);
		SimTK::Transform textOffsetF(SimTK::Rotation(), textOffset);
		decorativeTextF.setTransform(G_X_F * textOffsetF);
		decorativeTextF.setScaleFactors(SimTK::Vec3(textScaleFactor,
			textScaleFactor, textScaleFactor));
		decorativeTextF.setColor(textColor);
		geometry.push_back(decorativeTextF);
}



// Draw DuMM based geometry
void ParaMolecularDecorator::drawDummBasedGeometry(Array_<DecorativeGeometry>& geometry,
const State& someState)
{
	for (SimTK::DuMM::BondIndex bIx(0); bIx < dumm->getNumBonds(); ++bIx) {
		const SimTK::DuMM::AtomIndex dAIx1 = dumm->getBondAtom(bIx, 0);
		const SimTK::DuMM::AtomIndex dAIx2 = dumm->getBondAtom(bIx, 1);
		const SimTK::MobilizedBodyIndex mbx1 = dumm->getAtomBody(dAIx1);
		const SimTK::MobilizedBodyIndex mbx2 = dumm->getAtomBody(dAIx2);
		const SimTK::MobilizedBody mobod1 = matter->getMobilizedBody(mbx1);
		const SimTK::MobilizedBody mobod2 = matter->getMobilizedBody(mbx2);

		SimTK::Transform X_GB1 = mobod1.getBodyTransform(someState);
		SimTK::Transform X_GB2 = mobod2.getBodyTransform(someState);

		SimTK::Vec3 p_BS1 = dumm->getAtomStationOnBody(dAIx1);
		SimTK::Vec3 p_GS1 = X_GB1 * p_BS1;

		SimTK::Vec3 p_BS2 = dumm->getAtomStationOnBody(dAIx2);
		SimTK::Vec3 p_GS2 = X_GB2 * p_BS2;


		if( mbx1 == mbx2 ){
			geometry.push_back(DecorativeLine( p_GS1, p_GS2 ));
			(geometry.back()).setLineThickness(3);
			(geometry.back()).setColor( SimTK::Gray );
		}else{
			geometry.push_back(DecorativeLine( p_GS1, p_GS2 ));
			(geometry.back()).setLineThickness(2);
			(geometry.back()).setColor( SimTK::Black );
		}
	}
}

void ParaMolecularDecorator::generateDecorations(const State& someState,
		Array_<DecorativeGeometry>& geometry) 
{

	// Draw loaded arrows
	drawLoadedArrows(geometry);

	// Draw DuMM based geometry
	drawDummBasedGeometry(geometry, someState);

	//   
	/*
	// DuMM
	for (DuMM::AtomIndex daIx(0); daIx < dumm->getNumAtoms(); ++daIx) {
		const SimTK::MobilizedBodyIndex mbx = dumm->getAtomBody(daIx);
		const SimTK::MobilizedBody mobod = matter->getMobilizedBody(mbx);

		SimTK::Transform X_GB = mobod.getBodyTransform(someState);
		SimTK::Vec3 p_BS = dumm->getAtomStationOnBody(daIx);
		SimTK::Vec3 p_GS = X_GB * p_BS;
		SimTK::Transform X_BD(Rotation(), p_GS);

		//Real shrink = 0.3;
		//Real opacity = dumm->getAtomElement(daIx)==1?0.5:1;
		Real opacity = 0.5;
		Real r = dumm->getAtomRadius(daIx);
		if (r < 0.01){
			r = 0.1; //nm
		}

		geometry.push_back( DecorativeSphere(shrink * r) );
		(geometry.back()).setColor(dumm->getAtomDefaultColor(daIx));
		(geometry.back()).setOpacity(opacity);
		(geometry.back()).setResolution(3);
		(geometry.back()).setTransform(X_BD);
		

		// Text
		//std::ostringstream streamObj;
		//streamObj << std::fixed;
		//streamObj << std::setprecision(3);
		//streamObj //<< X_BD.p()[0] << ' ' << X_BD.p()[1] << ' ' 
		//    << X_BD.p()[2];

		//std::string text1 = streamObj.str();
		//DecorativeText decorativeText1(text1);
		//decorativeText1.setTransform(X_BD);
		//decorativeText1.setScaleFactors(SimTK::Vec3(0.02, 0.02, 0.02));
		//decorativeText1.setColor(SimTK::Vec3(1, 0, 1));
		//geometry.push_back(decorativeText1);

	}
// */

	SimTK::Transform G_X_T = molecules[0]->getTopLevelTransform();

    /*  DecorativeFrame decorativeFrameT;
			decorativeFrameT.setTransform(G_X_T);
			decorativeFrameT.setScaleFactors(SimTK::Vec3(0.06, 0.06, 0.06));
			decorativeFrameT.setLineThickness(5);
			decorativeFrameT.setColor(SimTK::Vec3(0.5, 0.5, 0));
			geometry.push_back( decorativeFrameT ); */

	// Draw Compound transforms for root atoms NEW WAY
	/* 		
			SimTK::Transform P_X_F[matter->getNumBodies()]; // related to X_PFs
			SimTK::Transform T_X_root[matter->getNumBodies()]; // related to CompoundAtom.frameInMobilizedBodyFrame s
			SimTK::Transform T_X_Proot;
			SimTK::Transform root_X_M0[matter->getNumBodies()];
			SimTK::Angle inboardBondDihedralAngles[matter->getNumBodies()]; */ // related to X_FMs
			//SimTK::Real inboardBondLengths[matter->getNumBodies()]; // related to X_FMs
			//SimTK::Vec3 locs[molecules[0]->getNumAtoms()];

 	/*

			// Iterate through atoms - get T_X_roots for all the bodies
			for (SimTK::Compound::AtomIndex aIx(0); aIx < molecules[0]->getNumAtoms(); ++aIx){
				if(molecules[0]->getAtomLocationInMobilizedBodyFrame(aIx) == 0){ // atom is at body's origin
					// Get body
					SimTK::MobilizedBodyIndex mbx = molecules[0]->getAtomMobilizedBodyIndex(aIx);
					T_X_root[int(mbx)] = molecules[0]->calcDefaultAtomFrameInCompoundFrame(aIx);
				}
			}

			P_X_F[1] = G_X_T * T_X_root[1]; // NEW
			// Iterate through atoms - get P_X_F for all the bodies
			for (SimTK::Compound::AtomIndex aIx(1); aIx < molecules[0]->getNumAtoms(); ++aIx){
				if(molecules[0]->getAtomLocationInMobilizedBodyFrame(aIx) == 0){ // atom is at body's origin

					// Get body, parentBody, parentAtom
					SimTK::MobilizedBodyIndex mbx = molecules[0]->getAtomMobilizedBodyIndex(aIx);
					const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
					const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
					SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

					if(parentMobod.getMobilizedBodyIndex() != 0){ // parent not Ground
					    T_X_Proot = T_X_root[parentMbx];

					    // Get inboard dihedral angle and put in root_X_M0
					    inboardBondDihedralAngles[int(mbx)] = molecules[0]->bgetDefaultInboardDihedralAngle(aIx);
					    //inboardBondLengths[int(mbx)] = molecules[0]->bgetDefaultInboardBondLength(aIx);
					    root_X_M0[int(mbx)] = SimTK::Transform(
					        SimTK::Rotation(inboardBondDihedralAngles[int(mbx)], SimTK::XAxis)
					    );

					    // Get Proot_X_M0
					    SimTK::Transform T_X_M0 = T_X_root[int(mbx)] * root_X_M0[int(mbx)];
					    SimTK::Transform Proot_X_T = ~T_X_Proot;
					    SimTK::Transform Proot_X_M0 = Proot_X_T * T_X_M0;
					    P_X_F[int(mbx)] = Proot_X_M0 * M_X_pin;

					    // Draw root 
					    SimTK::Transform G_X_root = G_X_T * T_X_root[int(mbx)];
					    //std::ostringstream streamObj_r;
					    //streamObj_r << std::string("r") + std::to_string(int(mbx));
					    //std::string text_r = streamObj_r.str();
					    //DecorativeText decorativeText_r(text_r);
					    //SimTK::Transform textOffset_r(SimTK::Rotation(), SimTK::Vec3(0.01, 0.0, 0.0));
					    //decorativeText_r.setTransform(G_X_root * textOffset_r);
					    //decorativeText_r.setScaleFactors(SimTK::Vec3(0.008, 0.008, 0.008));
					    //decorativeText_r.setColor(SimTK::Vec3(0, 1, 0));
					    //geometry.push_back(decorativeText_r);
			
					    DecorativeFrame decorativeFrame_r;
					    //decorativeFrame_r.setTransform(G_X_root * textOffset_r);
					    decorativeFrame_r.setTransform(G_X_root);
					    decorativeFrame_r.setScaleFactors(SimTK::Vec3(0.04, 0.04, 0.04));
					    decorativeFrame_r.setLineThickness(4);
					    decorativeFrame_r.setColor(SimTK::Vec3(0, 1, 0));
					    geometry.push_back( decorativeFrame_r );

					    // Draw F (Proot_X_Mr in this case) in Ground recovered from P_X_F
					    SimTK::Transform G_X_Proot = G_X_T * T_X_Proot;
					    SimTK::Transform G_X_F = G_X_Proot * P_X_F[int(mbx)];

					    std::ostringstream streamObjf;
					    streamObjf << std::string("a") + std::to_string(int(aIx));
					    std::string textf = streamObjf.str();
					    DecorativeText decorativeTextf(textf);
					    SimTK::Transform textOffsetf(SimTK::Rotation(), SimTK::Vec3(0.01, 0.0, 0.0));
					    decorativeTextf.setTransform(G_X_F * textOffsetf);
					    decorativeTextf.setScaleFactors(SimTK::Vec3(0.008, 0.008, 0.008));
					    decorativeTextf.setColor(SimTK::Vec3(0, 0, 1));
					    geometry.push_back(decorativeTextf);
			
					    DecorativeFrame decorativeFramef;
					    //decorativeFramef.setTransform(G_X_F * textOffsetf);
					    decorativeFramef.setTransform(G_X_F);
					    decorativeFramef.setScaleFactors(SimTK::Vec3(0.04, 0.04, 0.04));
					    decorativeFramef.setLineThickness(4);
					    decorativeFramef.setColor(SimTK::Vec3(0, 0, 1));
					    geometry.push_back( decorativeFramef );

					    // Line between P and F in Ground
					    //DecorativeLine decorativeLinepf(G_X_Proot.p() + textOffsetf.p(), G_X_F.p() + textOffsetf.p());
//                        DecorativeLine decorativeLinepf(G_X_Proot.p(), G_X_F.p());
//                        decorativeLinepf.setLineThickness(3);
//                        geometry.push_back(decorativeLinepf);

					} //END if parent not Ground
				}
			}
// */

	// Test Transforms operations

	/*    DecorativeFrame decorativeFrame_G;
	decorativeFrame_G.setScaleFactors(SimTK::Vec3(1, 1, 1));
	decorativeFrame_G.setColor(SimTK::Vec3(0, 0, 0));
	geometry.push_back( decorativeFrame_G );

	SimTK::Vec3 v1(1,0,0);

	DecorativeLine decorativeLine_v1(Vec3(0,0,0), v1);
	decorativeLine_v1.setLineThickness(2);
	decorativeLine_v1.setColor(SimTK::Vec3(0, 0, 0));
	geometry.push_back(decorativeLine_v1);

	SimTK::Angle bond12Angle = 120 * Deg2Rad;
	SimTK::Vec3 v2 = SimTK::Rotation(120 * Deg2Rad, ZAxis) * v1;

	DecorativeLine decorativeLine_v2(Vec3(0,0,0), v2);
	decorativeLine_v2.setLineThickness(2);
	decorativeLine_v2.setColor(SimTK::Vec3(1, 0, 1));
	geometry.push_back(decorativeLine_v2);

	SimTK::Angle dihedral3 = 120 * Deg2Rad;
	SimTK::Angle dihedral4 = -120 * Deg2Rad;
	SimTK::Vec3 v3 = SimTK::Rotation(dihedral3, YAxis) * SimTK::Rotation(bond12Angle, ZAxis) * v1;
	SimTK::Vec3 v4 = SimTK::Rotation(dihedral4, YAxis) * SimTK::Rotation(bond12Angle, ZAxis) * v1;

	DecorativeLine decorativeLine_v3(Vec3(0,0,0), v3);
	decorativeLine_v3.setLineThickness(2);
	decorativeLine_v3.setColor(SimTK::Vec3(1, 0, 0));
	geometry.push_back(decorativeLine_v3);

	DecorativeLine decorativeLine_v4(Vec3(0,0,0), v4);
	decorativeLine_v4.setLineThickness(2);
	decorativeLine_v4.setColor(SimTK::Vec3(0, 0, 1));
	geometry.push_back(decorativeLine_v4);*/
	
	/*
	 // Ground
	 SimTK::Transform offset(SimTK::Rotation(), SimTK::Vec3(0.0, 0.01, 0.0));
	 SimTK::Transform G_X_G;
	 DecorativeFrame decorativeFrame_G;
	 decorativeFrame_G.setTransform(G_X_G);
	 decorativeFrame_G.setScaleFactors(SimTK::Vec3(1, 1, 1));
	 decorativeFrame_G.setColor(SimTK::Vec3(0, 0, 0));
	 geometry.push_back( decorativeFrame_G );
 
	 // v1
	 SimTK::Vec3 G_v1(0.5, -1.31, 1.1);
	 DecorativeLine decorativeLine_v1(G_X_G.p() , G_v1);
	 decorativeLine_v1.setLineThickness(2);
	 decorativeLine_v1.setColor(SimTK::Vec3(0, 0, 0));
	 geometry.push_back(decorativeLine_v1);
 
	 // F1
	 SimTK::Transform G_X_F1( SimTK::Rotation(1.2, SimTK::UnitVec3(2.3, 4.1, -0.9)), SimTK::Vec3(1, 1, 1) );
	 DecorativeFrame decorativeFrame_F1;
	 decorativeFrame_F1.setTransform(G_X_F1);
	 decorativeFrame_F1.setScaleFactors(SimTK::Vec3(1, 1, 1));
	 decorativeFrame_F1.setColor(SimTK::Vec3(0, 0, 1));
	 geometry.push_back( decorativeFrame_F1 );
   
	 // Check v1 in F1
	 SimTK::Vec3 F1_v1 = ~(G_X_F1.R()) * G_v1;
	 SimTK::Vec3 checkG_v1 = G_X_F1.R() * F1_v1;
	 DecorativeLine decorativeLine_checkG_v1(G_X_G.p() + offset.p(), checkG_v1 + offset.p());
	 decorativeLine_checkG_v1.setLineThickness(2);
	 decorativeLine_checkG_v1.setColor(SimTK::Vec3(0, 1, 0));
	 geometry.push_back(decorativeLine_checkG_v1);
   
	 // v1 F1 Bond
	 DecorativeLine decorativeLine_v1F1(checkG_v1, G_X_F1.p());
	 decorativeLine_v1F1.setLineThickness(2);
	 decorativeLine_v1F1.setColor(SimTK::Vec3(0, 0, 1));
	 geometry.push_back(decorativeLine_v1F1);
   
	 // Transform 3
	 SimTK::Transform F1_X_F3 = alignFlipAndTranslateFrameAlongXAxis(G_X_F1, G_v1);
	 SimTK::Transform G_X_F3 = G_X_F1 * F1_X_F3;
	 DecorativeFrame decorativeFrame_F3;
	 decorativeFrame_F3.setTransform(G_X_F3);
	 decorativeFrame_F3.setScaleFactors(SimTK::Vec3(1., 1., 1.));
	 decorativeFrame_F3.setColor(SimTK::Vec3(1, 0, 0));
	 geometry.push_back( decorativeFrame_F3 );
 
 // */

	///*
	// Draw Compound transforms for periferic atoms NEW WAY
	// Set transforms inside the bodies = root_X_atom.p; Set locations for everyone
	/*
	for (SimTK::Compound::AtomIndex aIx(1); aIx < molecules[0]->getNumAtoms(); ++aIx){
		SimTK::MobilizedBodyIndex mbx = molecules[0]->getAtomMobilizedBodyIndex(aIx);
		if(molecules[0]->getAtomLocationInMobilizedBodyFrame(aIx) != 0){ // atom is not at body's origin
			SimTK::Transform G_X_root = G_X_T * T_X_root[int(mbx)]; // NEW
			assert(atomTargets[aIx][0]); // We have to set it in World
			SimTK::Vec3 G_vchild = atomTargets[aIx]; //NEW
			SimTK::Vec3 G_vroot = G_X_root.p(); // NEW
			SimTK::Vec3 G_v = G_vchild - G_vroot; // NEW
			SimTK::Vec3 root_v = ~(G_X_root.R()) * G_v; //NEW
			SimTK::Vec3 mroot_v = -1 * root_v;
			SimTK::UnitVec3 G_XAxis(1,0,0);
			SimTK::UnitVec3 root_XAxis = ~(G_X_root.R()) * G_XAxis;

			SimTK::Transform root_X_child = alignFlipAndTranslateFrameAlongXAxis(G_X_root, G_vchild);
			SimTK::Transform G_X_child = G_X_root * root_X_child;

			// Draw F (Proot_X_Mr in this case) in Ground recovered from P_X_F
			std::ostringstream streamObj_c;
			streamObj_c << std::string("c") + std::to_string(int(aIx));
			std::string text_c = streamObj_c.str();
			DecorativeText decorativeText_c(text_c);
			SimTK::Transform textOffset_c(SimTK::Rotation(), SimTK::Vec3(0.0, 0.01, 0.0));
			//decorativeText_c.setTransform(G_X_child * textOffset_c);
			decorativeText_c.setTransform(G_X_child);
			decorativeText_c.setScaleFactors(SimTK::Vec3(0.008, 0.008, 0.008));
			decorativeText_c.setColor(SimTK::Vec3(1, 0, 1));
			geometry.push_back(decorativeText_c);

			// Text
			//std::ostringstream streamObj;
			//streamObj << std::fixed;
			//streamObj << std::setprecision(3);
			//streamObj << G_X_child.p()[2];
			//std::string text1 = streamObj.str();
			//DecorativeText decorativeText1(text1);
			//decorativeText1.setTransform(G_X_child);
			//decorativeText1.setScaleFactors(SimTK::Vec3(0.02, 0.02, 0.02));
			//decorativeText1.setColor(SimTK::Vec3(0, 0, 0));
			//geometry.push_back(decorativeText1);

			DecorativeFrame decorativeFrame_c;
			//decorativeFrame_c.setTransform(G_X_child * textOffset_c);
			decorativeFrame_c.setTransform(G_X_child);
			decorativeFrame_c.setScaleFactors(SimTK::Vec3(0.04, 0.04, 0.04));
			decorativeFrame_c.setLineThickness(4);
			decorativeFrame_c.setColor( SimTK::Gray );
			geometry.push_back( decorativeFrame_c );

			// Line between root and child in Ground
			//DecorativeLine decorativeLine_rc2(G_vroot, G_vchild);
			//decorativeLine_rc2.setLineThickness(2);
			//decorativeLine_rc2.setColor(SimTK::Vec3(0, 1, 1));
			//geometry.push_back(decorativeLine_rc2);

		}
	}
// */

	// /*
	// Draw Default Compound 
	//DecorativeBrick topDecorativeBrick;
	//topDecorativeBrick.setTransform(G_X_T);
	//topDecorativeBrick.setScaleFactors(SimTK::Vec3(0.03, 0.03, 0.03));
	//topDecorativeBrick.setColor(SimTK::Vec3(0, 0, 0));
	//geometry.push_back( topDecorativeBrick );
	//for (SimTK::Compound::AtomIndex aIx(0); aIx < molecules[0]->getNumAtoms(); ++aIx){
	//    SimTK::Transform T_X_atom =  molecules[0]->calcDefaultAtomFrameInCompoundFrame(SimTK::Compound::AtomIndex(aIx));
	//    SimTK::Transform G_X_atom = G_X_T * T_X_atom;
		// Bricks
		//DecorativeBrick decorativeBrick(SimTK::Vec3(0.03, 0.03, 0.03));
		//decorativeBrick.setTransform(G_X_atom);
		//decorativeBrick.setOpacity(0.5);
		//geometry.push_back( decorativeBrick );
	
	// Text
		//std::ostringstream streamObj;
		//streamObj << std::fixed;
		//streamObj << std::setprecision(3);
		//streamObj << int(aIx)
			//<< G_X_atom.p()[0] << ' ' << G_X_atom.p()[1] << ' ' 
			//<< G_X_atom.p()[2]
	//;

	//std::string text1 = streamObj.str();
		//DecorativeText decorativeText1(text1);
		//decorativeText1.setTransform(G_X_atom);
		//decorativeText1.setScaleFactors(SimTK::Vec3(0.02, 0.02, 0.02));
		//decorativeText1.setColor(SimTK::Vec3(0, 0, 0));
		//geometry.push_back(decorativeText1);
		
	//}
// */

	// Ground frame
	// drawFrame(geometry, Transform(),
	// 	0.05, 4, SimTK::Vec3(0.5, 0.5, 0.5),
	// 	"G", 0.009, SimTK::Vec3(0.5, 0.5, 0.5), SimTK::Vec3(0.02, 0.0, 0.0));

	// Top frame	
	// drawFrame(geometry, G_X_T,
	// 	0.05, 4, SimTK::Vec3(0.5, 0.5, 0.5),
	// 	"Top", 0.009, SimTK::Vec3(0.5, 0.0, 0.5), SimTK::Vec3(0.02, 0.0, 0.0));

	for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();

		SimTK::Transform G_X_B = mobod.getBodyTransform(someState);
		
//        // Bricks
//        DecorativeBrick decorativeBrick(SimTK::Vec3(0.03, 0.03, 0.03));
//        decorativeBrick.setTransform(G_X_B);
//        //decorativeBrick.setColor(SimTK::Vec3(10, 0, 0));
//        decorativeBrick.setOpacity(0.5);
//        geometry.push_back( decorativeBrick );

		if((mbx > 0)){
		//if(int(mbx) == 3){
			SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
			const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();

			// Get all possible transforms
			SimTK::Transform B_X_M = mobod.getOutboardFrame(someState);
			SimTK::Transform G_X_P = parentMobod.getBodyTransform(someState);
			SimTK::Transform G_X_B = mobod.getBodyTransform(someState);
			SimTK::Transform G_X_M = G_X_B * B_X_M;

			SimTK::Transform P_X_F = mobod.getInboardFrame(someState);
			SimTK::Transform G_X_F = G_X_P * P_X_F;
			SimTK::Transform F_X_M = mobod.getMobilizerTransform(someState);

			// F
			//DecorativeSphere decorativeSphereF(0.02);
			//decorativeSphereF.setColor(SimTK::Vec3(0, 0, 1));
			//decorativeSphereF.setOpacity(0.5);
			//decorativeSphereF.setTransform(G_X_F);
			//geometry.push_back(decorativeSphereF);

			// M
			//DecorativeSphere decorativeSphereM(0.02);
			//decorativeSphereM.setColor(SimTK::Vec3(1, 0, 0));
			//decorativeSphereM.setOpacity(0.5);
			//decorativeSphereM.setTransform(G_X_M);
			//geometry.push_back(decorativeSphereM);

			// Draw frames
			std::vector<int> chosenBodies = {3, 6};
			/* {        1,  2,  3,  4,  5, 6,  7,  8,  9
			, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
			, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29}; */
			/* if(chosenBodies.size() >= matter->getNumBodies()){
				std::cerr << "Visualizer chosenBody greater than allowed.\n"
					<< std::flush;
				exit(1);
			} */

			bool found = false;
			found = 
				(std::find(std::begin(chosenBodies), std::end(chosenBodies), int(mbx))
				!= std::end(chosenBodies));

			//if(found) {
			if( (int(mbx) >= 1) ) {
			//if(false){

				// // Frame P
				// std::ostringstream streamObjP;
				// streamObjP << std::string("P") + std::to_string(int(mbx));
				// std::string textP = streamObjP.str();
				// drawFrame(geometry, G_X_P,
				// 	0.05, 4, SimTK::Vec3(0, 0, 0),
				// 	streamObjP.str(), 0.009, SimTK::Vec3(0, 0, 0), SimTK::Vec3(-0.02, 0.0, 0.0));
				// Frame F
				std::ostringstream streamObjF;
				streamObjF << std::string("F") + std::to_string(int(mbx)); //+ " " + std::to_string(this->FCommVar)
				std::string textF = streamObjF.str();
				drawFrame(geometry, G_X_F,
					0.04, 4, SimTK::Vec3(0, 0, 1),
					streamObjF.str(), 0.008, SimTK::Vec3(0, 0, 1), SimTK::Vec3(-0.02, 0.0, 0.0));
				// // Frame M
				// std::ostringstream streamObjM;
				// std::setprecision(2);
				// streamObjM << std::string("M") + std::to_string(int(mbx));
				// drawFrame(geometry, G_X_M,
				// 	0.05, 4, SimTK::Vec3(1, 0, 0),
				// 	streamObjM.str(), 0.008, SimTK::Vec3(1, 0, 0), SimTK::Vec3(-0.03, 0.0, 0.0));
				// // Frame B
				// std::ostringstream streamObjB;
				// streamObjB << std::string("B") + std::to_string(int(mbx)); //+ " " + std::to_string(this->BCommVar);
				// std::string textB = streamObjB.str();
				// drawFrame(geometry, G_X_B,
				// 	0.04, 4, SimTK::Vec3(0, 0, 0),
				// 	streamObjB.str(), 0.008, SimTK::Vec3(0, 0, 0), SimTK::Vec3(-0.04, 0.0, 0.0));

			}

			// Draw lines
			//if(found) {
			if( (int(mbx) >= 1) ) {
				// BM expressed in Ground
				drawLine(geometry, G_X_B, G_X_M,
					mobod.getNumU(someState), 4);

			}

		}

	}


}


/*
ParaMolecularDecorator::~ParaMolecularDecorator()
{
	points.clear();
	lines.clear();
}
*/












