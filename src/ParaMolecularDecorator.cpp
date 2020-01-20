/**@file
Implementation of HamiltonianMonteCarloSampler class. **/

#include "Robo.hpp"
#include "ParaMolecularDecorator.hpp"

using namespace SimTK;

ParaMolecularDecorator::ParaMolecularDecorator(SimTK::CompoundSystem *argCompoundSystem,
                 SimTK::SimbodyMatterSubsystem *argMatter,
                 //Topology *argResidue, RE
                 SimTK::DuMMForceFieldSubsystem *argDumm,
                 SimTK::GeneralForceSubsystem *argForces)
{
    this->compoundSystem = argCompoundSystem;
    this->matter = argMatter;
    // this->molecule = argResidue; //RE
    this->dumm = argDumm;
    this->forces = argForces;
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

void ParaMolecularDecorator::clearPoints(void)
{
    points.clear();
}

void ParaMolecularDecorator::clearLines(void)
{
    lines.clear();
}

// Gmolmodel specific
void ParaMolecularDecorator::setAtomTargets(std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>> residueAtomLocations)
{
    atomTargets.clear();
    for(unsigned int j = 0; j < residueAtomLocations.size(); j++){
        SimTK::Compound::AtomIndex atomIndex = ((residueAtomLocations[j]).first)->atomIndex;
        SimTK::Vec3 location = ((residueAtomLocations[j]).second);
        atomTargets.insert(std::pair<SimTK::Compound::AtomIndex, SimTK::Vec3>(atomIndex, location));
    }
}

void ParaMolecularDecorator::generateDecorations(const State& someState,
        Array_<DecorativeGeometry>& geometry) 
{
    /*
    for (auto p = points.begin(); p != points.end(); p++) {
        geometry.push_back(DecorativePoint(*p));
    }
    for (auto p = lines.begin(); p != lines.end(); p++) {
        geometry.push_back(DecorativeLine( p->first, p->second ));
        (geometry.back()).setLineThickness(5);
        (geometry.back()).setColor( Vec3(1, 0, 1) );
    }
    */

    // Draw DuMM based geometry
 //   /*
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
//    */
 /*
    // DuMM
    for (DuMM::AtomIndex daIx(0); daIx < dumm->getNumAtoms(); ++daIx) {
        const SimTK::MobilizedBodyIndex mbx = dumm->getAtomBody(daIx);
        const SimTK::MobilizedBody mobod = matter->getMobilizedBody(mbx);

        SimTK::Transform X_GB = mobod.getBodyTransform(someState);
        SimTK::Vec3 p_BS = dumm->getAtomStationOnBody(daIx);
        SimTK::Vec3 p_GS = X_GB * p_BS;
        SimTK::Transform X_BD(Rotation(), p_GS);

        Real shrink = 0.3;
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

            DecorativeFrame decorativeFrameT;
            decorativeFrameT.setTransform(G_X_T);
            decorativeFrameT.setScaleFactors(SimTK::Vec3(0.06, 0.06, 0.06));
            decorativeFrameT.setLineThickness(5);
            decorativeFrameT.setColor(SimTK::Vec3(0.5, 0, 0));
            geometry.push_back( decorativeFrameT );


    // Draw Compound transforms for root atoms NEW WAY
            SimTK::Transform M_X_pin = SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis); // Moves rotation from X to Z
            SimTK::Transform P_X_F[matter->getNumBodies()]; // related to X_PFs
            SimTK::Transform T_X_root[matter->getNumBodies()]; // related to CompoundAtom.frameInMobilizedBodyFrame s
            SimTK::Transform T_X_Proot; // NEW
            SimTK::Transform root_X_M0[matter->getNumBodies()];
            SimTK::Angle inboardBondDihedralAngles[matter->getNumBodies()]; // related to X_FMs
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

 ///*
    // Draw Rigid bodies
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        //const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        //const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();

        //SimTK::Transform G_X_B = mobod.getBodyTransform(someState);
        
        // Bricks
/*        DecorativeBrick decorativeBrick(SimTK::Vec3(0.03, 0.03, 0.03));
        decorativeBrick.setTransform(G_X_B);
        //decorativeBrick.setColor(SimTK::Vec3(10, 0, 0));
        decorativeBrick.setOpacity(0.5);
        geometry.push_back( decorativeBrick );

        std::ostringstream streamObjB;
        streamObjB << std::string("B") + std::to_string(int(mbx)); 
        std::string textB = streamObjB.str();
        DecorativeText decorativeTextB(textB);
        SimTK::Transform textOffsetB(SimTK::Rotation(), SimTK::Vec3(0.0, 0.0, -0.01));
        decorativeTextB.setTransform(G_X_B * textOffsetB);
        decorativeTextB.setScaleFactors(SimTK::Vec3(0.008, 0.008, 0.008));
        decorativeTextB.setColor(SimTK::Vec3(0, 0, 0));
        geometry.push_back(decorativeTextB);

        DecorativeFrame decorativeFrame;
        decorativeFrame.setTransform(G_X_B);
        decorativeFrame.setScaleFactors(SimTK::Vec3(0.03, 0.03, 0.03));
        decorativeFrame.setLineThickness(4); 
        decorativeFrame.setColor(SimTK::Vec3(0, 0, 0));
        decorativeFrame.setRepresentation(SimTK::DecorativeGeometry::Representation::DrawPoints); 
        geometry.push_back( decorativeFrame );*/

        if((mbx > 0)){
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
            //SimTK::Transform G_X_M = G_X_F * F_X_M; // = G_X_B * B_X_M;

            DecorativeFrame decorativeFrameB;
            decorativeFrameB.setTransform(G_X_B);
            decorativeFrameB.setScaleFactors(SimTK::Vec3(0.04, 0.04, 0.04));
            decorativeFrameB.setLineThickness(4);
            decorativeFrameB.setColor(SimTK::Vec3(0, 0, 0));
            geometry.push_back( decorativeFrameB );

            // F
            //DecorativeSphere decorativeSphereF(0.02);
            //decorativeSphereF.setColor(SimTK::Vec3(0, 0, 1));
            //decorativeSphereF.setOpacity(0.5);
            //decorativeSphereF.setTransform(G_X_F);
            //geometry.push_back(decorativeSphereF);

            //X_PF, X_BM frames and text
	    /*
            DecorativeFrame decorativeFrameF;
            decorativeFrameF.setTransform(G_X_F);
            decorativeFrameF.setScaleFactors(SimTK::Vec3(0.04, 0.04, 0.04));
            decorativeFrameF.setLineThickness(4);
            decorativeFrameF.setColor(SimTK::Vec3(0, 0, 1));
            geometry.push_back( decorativeFrameF );

            std::ostringstream streamObjF;
            streamObjF << std::string("F") + std::to_string(int(mbx)); 
            std::string textF = streamObjF.str();
            DecorativeText decorativeTextF(textF);
            SimTK::Transform textOffsetF(SimTK::Rotation(), SimTK::Vec3(-0.01, 0.0, 0.0));
            decorativeTextF.setTransform(G_X_F * textOffsetF);
            decorativeTextF.setScaleFactors(SimTK::Vec3(0.008, 0.008, 0.008));
            decorativeTextF.setColor(SimTK::Vec3(0, 0, 1));
            geometry.push_back(decorativeTextF);
	    */

            // Line between P and F
            //DecorativeLine decorativeLinePF(G_X_P.p(), G_X_F.p());
            //decorativeLinePF.setLineThickness(3);
            //decorativeLinePF.setColor(SimTK::Vec3(1, 0, 0));
            //geometry.push_back(decorativeLinePF);

            // M
            //DecorativeSphere decorativeSphereM(0.02);
            //decorativeSphereM.setColor(SimTK::Vec3(1, 0, 0));
            //decorativeSphereM.setOpacity(0.5);
            //decorativeSphereM.setTransform(G_X_M);
            //geometry.push_back(decorativeSphereM);

	    /*
            DecorativeFrame decorativeFrameM;
            decorativeFrameM.setTransform(G_X_M);
            decorativeFrameM.setScaleFactors(SimTK::Vec3(0.05, 0.05, 0.05));
            decorativeFrameM.setLineThickness(4);
            decorativeFrameM.setColor(SimTK::Vec3(1, 0, 0));
            geometry.push_back( decorativeFrameM );

            std::ostringstream streamObjM;
            streamObjM << std::string("M") + std::to_string(int(mbx)); 
            std::string textM = streamObjM.str();
            DecorativeText decorativeTextM(textM);
            SimTK::Transform textOffsetM(SimTK::Rotation(), SimTK::Vec3(0.0, -0.01, 0.0));
            decorativeTextM.setTransform(G_X_M * textOffsetM);
            decorativeTextM.setScaleFactors(SimTK::Vec3(0.008, 0.008, 0.008));
            decorativeTextM.setColor(SimTK::Vec3(1, 0, 0));
            geometry.push_back(decorativeTextM);
	    */

            /*
            // X_PF, X_BM lines
            DecorativeLine decorativeLinePF(G_X_P.p(), G_X_F.p());
            decorativeLinePF.setLineThickness(5);
            if(mobod.getNumU(someState) == 3){
                decorativeLinePF.setColor(SimTK::Vec3(SimTK::Blue));
	        }else if(mobod.getNumU(someState) == 5){
                decorativeLinePF.setColor(SimTK::Vec3(SimTK::Blue));
	        }else if(mobod.getNumU(someState) == 2){
                decorativeLinePF.setColor(SimTK::Vec3(SimTK::Blue));
            }else if (mobod.getNumU(someState) == 1){
                decorativeLinePF.setColor(SimTK::Vec3(SimTK::Blue));
            }else{
                decorativeLinePF.setColor(SimTK::Vec3(SimTK::Blue));
            }
            geometry.push_back(decorativeLinePF);

            DecorativeLine decorativeLineFM(G_X_F.p(), G_X_M.p());
            decorativeLineFM.setLineThickness(3);
            if(mobod.getNumU(someState) == 3){
                decorativeLineFM.setColor(SimTK::Vec3(SimTK::Red));
	        }else if(mobod.getNumU(someState) == 5){
                decorativeLineFM.setColor(SimTK::Vec3(SimTK::Magenta));
            }else if (mobod.getNumU(someState) == 2){
                decorativeLineFM.setColor(SimTK::Vec3(SimTK::Cyan));
            }else if (mobod.getNumU(someState) == 1){
                decorativeLineFM.setColor(SimTK::Vec3(SimTK::Green));
            }else{
                decorativeLineFM.setColor(SimTK::Vec3(SimTK::Magenta));
            }
            geometry.push_back(decorativeLineFM);
            
            DecorativeLine decorativeLineBM(G_X_B.p(), G_X_M.p());
            decorativeLineBM.setLineThickness(4);
            if(mobod.getNumU(someState) == 3){
                decorativeLineBM.setColor(SimTK::Vec3(SimTK::Black));
	        }else if(mobod.getNumU(someState) == 5){
                decorativeLineBM.setColor(SimTK::Vec3(SimTK::Black));
            }else if (mobod.getNumU(someState) == 2){
                decorativeLineBM.setColor(SimTK::Vec3(SimTK::Black));
            }else if (mobod.getNumU(someState) == 1){
                decorativeLineBM.setColor(SimTK::Vec3(SimTK::Black));
            }else{
                decorativeLineBM.setColor(SimTK::Vec3(SimTK::Black));
            }
            geometry.push_back(decorativeLineBM);
            // */


        }

        }

}

ParaMolecularDecorator::~ParaMolecularDecorator()
{
    points.clear();
    lines.clear();
}













