#include "World.hpp"



// TODO write a pdb writer for all the Compounds
// TODO move this in Topology since they work only for one Compound
void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix)
{
  double mult = 10000*advanced.getTime(); // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix)
{
  double mult = 10000*advanced.getTime(); // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime)
{
  double mult = 10000*aTime; // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime)
{
  double mult = 10000*aTime; // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(SimTK::PdbStructure pdb, const char *FN)
{
  std::filebuf fb;
  fb.open(FN, std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); //automatically multiplies by ten (nm to A)
  fb.close();
}

/** Print a Compound Cartesian coordinates as given by
     * Compound::calcAtomLocationInGroundFrame **/
void World::printPoss(const SimTK::Compound& c, SimTK::State& advanced)
{
    SimTK::Vec3 vertex;
    std::cout<<"Positions:"<<std::endl;
    for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
        vertex   = c.calcAtomLocationInGroundFrame(advanced, aIx);
        std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
          <<"["<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<"]"<<std::endl;
    }
}

void World::printVels(const SimTK::Compound& c, SimTK::State& advanced)
{
    SimTK::Vec3 vel;
    std::cout<<"Velocities:"<<std::endl;
    for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
        vel      = c.calcAtomVelocityInGroundFrame(advanced, aIx);
        std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
          <<"["<<vel[0]<<" "<<vel[1]<<" "<<vel[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
}

void World::printPossVels(const SimTK::Compound& c, SimTK::State& advanced)
{
    SimTK::Vec3 vertex, vel;
    std::cout<<"Positions:"<<std::endl;
    for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
        vertex   = c.calcAtomLocationInGroundFrame(advanced, aIx);
        std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
          <<"["<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<"Velocities:"<<std::endl;
    for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
        vel      = c.calcAtomVelocityInGroundFrame(advanced, aIx);
        std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
          <<"["<<vel[0]<<" "<<vel[1]<<" "<<vel[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
}

/** Constructor. Initializes the following pretty much empty objects:
 *  - CompoundSystem,
 *      - SimbodyMatterSubsystem, GeneralForceSubsystem, DecorationSubsystem,
 *        Visualizer, Visualizer::Reporter, DuMMForceFieldSubsystem,
 *  - Integrator with a TimeStepper on top **/
World::World(int worldIndex, bool isVisual, SimTK::Real visualizerFrequency)
{
    // Get an index from a higher caller
    ownWorldIndex = worldIndex;

    // Initialize Fixman torque flag
    _useFixmanTorque = false;
  
    // Molmodel System derived from Simbody System
    compoundSystem = new SimTK::CompoundSystem;

    // Simbody subsystems (Minimum requirements)
    matter = new SimTK::SimbodyMatterSubsystem(*compoundSystem);
    forces = new SimTK::GeneralForceSubsystem(*compoundSystem);

    // Initialize Molmodel default ForceSubsystem (DuMM)
    forceField = new SimTK::DuMMForceFieldSubsystem(*compoundSystem);

    // Intialize an integrator and a TimeStepper to manage it
    integ = new SimTK::VerletIntegrator(*compoundSystem);
    ts = new SimTK::TimeStepper(*compoundSystem, *integ);

    // Set the visual flag and if true initialize a Decorations Subsystem,
    // a Visualizer and a Simbody EventReporter which interacts with the
    // Visualizer
    this->visual = isVisual;
    if(visual){
        decorations = new SimTK::DecorationSubsystem(*compoundSystem);
        visualizer = new SimTK::Visualizer(*compoundSystem);
        visualizerReporter = new SimTK::Visualizer::Reporter(*visualizer
                , visualizerFrequency);
        compoundSystem->addEventReporter( visualizerReporter );

        // Initialize a DecorationGenerator
        paraMolecularDecorator = new ParaMolecularDecorator(
                compoundSystem,
                matter,
                forceField,
                forces
        );

        visualizer->addDecorationGenerator(paraMolecularDecorator);

    }

    // Statistics
    moleculeCount = -1;
    nofSamples = 0;

    // Thermodynamics
    this->temperature = -1; // this leads to unusal behaviour hopefully

}

/** Creates Gmolmodel topologies objects and based on amberReader forcefield
 * adds parameters: defines Biotypes; - adds BAT parameters to DuMM. Also
 * creates decorations for visualizers **/
void World::AddMolecule(
        readAmberInput *amberReader,
        std::string rbFN,
        std::string flexFN,
        std::string regimenSpec,
        std::string argRoot)
{
    // Statistics
    moleculeCount++; // Used for unique names of molecules
 
    // Add a new molecule (Topology object which inherits Compound)
    // to the vector of molecules.
    // TODO: Why resName and moleculeName have to be the same?
    std::string moleculeName = regimenSpec + std::to_string(moleculeCount);
    Topology *top = new Topology(moleculeName);
    topologies.emplace_back(top);

    // Set atoms properties from a reader: number, name, element, initial
    // name, force field type, charge, coordinates, mass, LJ parameters
    (topologies.back())->SetGmolAtomPropertiesFromReader(amberReader);

    // Set bonds properties from reader: bond indeces, atom neighbours
    (topologies.back())->SetGmolBondingPropertiesFromReader(amberReader);

    // Set atoms Molmodel types (Compound::SingleAtom derived) based on
    // their valence
    //(topologies.back())->SetGmolAtomsMolmodelTypes();
    (topologies.back())->SetGmolAtomsMolmodelTypesTrial();

    // Add parameters from amberReader
    (topologies.back())->bAddAllParams(amberReader, *forceField);

    // Build the graph representing molecule's topology
    (topologies.back())->buildGraphAndMatchCoords(*forceField, std::stoi(argRoot));

    (topologies.back())->loadTriples();

    // Set flexibility according to the flexibility file
    (topologies.back())->setFlexibility(regimenSpec, flexFN);

    // Print Molmodel types
    (topologies.back())->PrintMolmodelAndDuMMTypes(*forceField);

    // Allocate the vector of coordinates (DCD)
    Xs.resize(Xs.size() + topologies.back()->getNAtoms());
    Ys.resize(Ys.size() + topologies.back()->getNAtoms());
    Zs.resize(Zs.size() + topologies.back()->getNAtoms());
  
    // Add Topology to CompoundSystem and realize topology
    compoundSystem->adoptCompound( *(topologies.back()) );
    std::cout<<"World::AddMolecule CompoundSystem adoptCompound "<< std::endl;
    (topologies.back())->setCompoundIndex(
            SimTK::CompoundSystem::CompoundIndex(
             compoundSystem->getNumCompounds() - 1));

    // Add the new molecule to Decorators's vector of molecules
    if(visual){
        // We need copy here.
        paraMolecularDecorator->AddMolecule(topologies.back());
    }

}

// Get the number of molecules
int World::getNofMolecules(void)
{
    return (this->moleculeCount + 1);
}

/** Calls CompoundSystem.modelOneCompound which links the Compounds to the
 * Simbody subsystems and realizes Topology. To be called after setting all
 * Compounds properties. **/
void World::modelTopologies(std::string GroundToCompoundMobilizerType)
{
    // Model the Compounds one by one in case we want to attach different types
    // of Mobilizers to the Ground in the feature.
    for ( unsigned int i = 0; i < this->topologies.size(); i++){
        //SimTK::String GroundToCompoundMobilizerType = "Free";
        //SimTK::String GroundToCompoundMobilizerType = "Weld";
        //SimTK::String GroundToCompoundMobilizerType = "Cartesian";


        this->rootMobility = GroundToCompoundMobilizerType;

        std::cout<<"World::ModelTopologies call to CompoundSystem::modelCompound " << i
                 << " grounded with mobilizer " << GroundToCompoundMobilizerType << std::endl;

        compoundSystem->modelOneCompound(
                SimTK::CompoundSystem::CompoundIndex(i),
                GroundToCompoundMobilizerType);


        // Realize Topology
        //compoundSystem->realizeTopology(); // restore MULMOL
        //((this->topologies)[i])->loadMobodsRelatedMaps(); // restore MULMOL

    }

    // Realize Topology
    //compoundSystem->realizeTopology();
}


const SimTK::State& World::realizeTopology(void)
{
        const SimTK::State& returnState = compoundSystem->realizeTopology();

	    //for ( unsigned int i = 0; i < this->topologies.size(); i++){
	    //    ((this->topologies)[i])->loadMobodsRelatedMaps();
	    //}

	return returnState;
}

void World::loadMobodsRelatedMaps(void)
{
	for ( unsigned int i = 0; i < this->topologies.size(); i++){
	    ((this->topologies)[i])->loadMobodsRelatedMaps();
        ((this->topologies)[i])->printMaps();
	}
}


/** Set up Fixman torque **/
void World::addFixmanTorque()
{
    // Set flag
    _useFixmanTorque = true;

    // Alloc memory for FixmanTorque implementation and add to forces
    FixmanTorqueImpl = new FixmanTorque(compoundSystem, *matter);
    ExtForce = new Force::Custom(*forces, FixmanTorqueImpl);
}

/** Check if the Fixman torque flag is set **/
bool World::isUsingFixmanTorque(void)
{
    return _useFixmanTorque;
}

/** Get writble pointer to FixmanTorque implementation **/
FixmanTorque * World::updFixmanTorque(void)
{
    return FixmanTorqueImpl;
}

/** Get pointer to FixmanTorque implementation **/
FixmanTorque * World::getFixmanTorque(void) const
{
    return FixmanTorqueImpl;
}

// ----------------------
// --- Thermodynamics ---
// ----------------------

/** Get the World temperature **/
SimTK::Real World::getTemperature(void)
{
    return this->temperature;
}

/** Set this World temperature but also ths samplers and 
Fixman torque temperature. **/
void World::setTemperature(SimTK::Real argTemperature)
{
    // Set the temperature for the samplers
    for(unsigned int samplerIx = 0; samplerIx < samplers.size(); samplerIx++){
        samplers[samplerIx]->setTemperature(argTemperature);
    }

    // Set the temperature for this World
    this->temperature = argTemperature;

    // Set the temperature for the Fixman torque also
    if(_useFixmanTorque){
        FixmanTorqueImpl->setTemperature(this->temperature);
    }
}
//...............

//...................
// --- Simulation ---
//...................

/** Get/Set seed for reproducibility. **/
void World::setSeed(int whichSampler, unsigned long long int argSeed)
{
    samplers[whichSampler]->setSeed(argSeed);
}

unsigned long long int World::getSeed(int whichSampler)
{
    return samplers[whichSampler]->getSeed();
}


/** Amber like scale factors. **/
void World::setAmberForceFieldScaleFactors()
{
    forceField->setVdw12ScaleFactor(0.0);
    forceField->setVdw13ScaleFactor(0.0);
    forceField->setVdw14ScaleFactor(0.5); // RESTORE from OpenMM
    //forceField->setVdw14ScaleFactor(0.0); // for OpenMM
    forceField->setVdw15ScaleFactor(1.0);

///* RESTORE SAFETY  forceField->setCoulomb12ScaleFactor(0.0);
    forceField->setCoulomb13ScaleFactor(0.0);
    forceField->setCoulomb14ScaleFactor(0.8333333333); // RESTORE from OpenMM
    //forceField->setCoulomb14ScaleFactor(0.0); // for OpenMM
    forceField->setCoulomb15ScaleFactor(1.0);
    //forceField->setVdwMixingRule(
    //       SimTK::DuMMForceFieldSubsystem::LorentzBerthelot); */
 
   /* DANGER ! no electrostatics
    forceField->setCoulomb12ScaleFactor(0.0);
    forceField->setCoulomb13ScaleFactor(0.0);
    forceField->setCoulomb14ScaleFactor(0.0);
    forceField->setCoulomb15ScaleFactor(0.0);
    // */
}

/** Set a global scaling factor for all the terms in the forcefield **/
void World::setGlobalForceFieldScaleFactor(SimTK::Real scaleFactor)
{
  
    forceField->setBondStretchGlobalScaleFactor(scaleFactor);
    forceField->setBondBendGlobalScaleFactor(scaleFactor);
    forceField->setBondTorsionGlobalScaleFactor(scaleFactor);
    forceField->setAmberImproperTorsionGlobalScaleFactor(scaleFactor);
  
    forceField->setVdw12ScaleFactor(scaleFactor);
    forceField->setVdw13ScaleFactor(scaleFactor);
    forceField->setVdw14ScaleFactor(scaleFactor);
    forceField->setVdw15ScaleFactor(scaleFactor);
    forceField->setVdwGlobalScaleFactor(scaleFactor);
  
    forceField->setCoulomb12ScaleFactor(scaleFactor);
    forceField->setCoulomb13ScaleFactor(scaleFactor);
    forceField->setCoulomb14ScaleFactor(scaleFactor);
    forceField->setCoulomb15ScaleFactor(scaleFactor);
    forceField->setCoulombGlobalScaleFactor(scaleFactor);

}

/** Set GBSA implicit solvent scale factor. **/
void World::setGbsaGlobalScaleFactor(SimTK::Real scaleFactor)
{
    forceField->setGbsaGlobalScaleFactor(scaleFactor);
}

/** Get a writeble pointer to the DuMM force field **/
SimTK::DuMMForceFieldSubsystem * World::updForceField()
{
    return forceField;
}

//...................
// --- Statistics ---
//...................

/** How many samples do we have so far **/
int World::getNofSamples(void)
{
    // Zero it every time the user asks
    nofSamples = 0;
 
    // Gather samples from all the samplers
    for(size_t i = 0; i < samplers.size(); i++){
        nofSamples += (samplers[i])->getNofSamples();
    }

    return this->nofSamples;
}

/** How many Samplers does this World have. **/
int World::getNofSamplers(void)
{
    return samplers.size();
}

/** Add a sampler to this World using the specialized struct
for samplers names. **/
int World::addSampler(SamplerName samplerName)
{
    if(samplerName == HMC){

        BaseSampler *p = new HamiltonianMonteCarloSampler(
                compoundSystem, matter, topologies[0],
                forceField, forces, ts
                );
        samplers.emplace_back(p);

    }

    return samplers.size();
}

// Get a sampler based on its position in the samplers vector
// TODO Use ampler polymorphism
const BaseSampler * World::getSampler(int which)
{
    return samplers[which];
}

// Get a writable sampler based on its position in the samplers vector
// TODO Use Sampler polymorphism
BaseSampler * World::updSampler(int which)
{
    return samplers[which];
}


/** Get a const reference to a molecule **/
const Topology& World::getTopology(int moleculeNumber) const{
    return *(topologies[moleculeNumber]);
}

/** Get a writble reference to the last molecule. **/
Topology& World::updTopology(int moleculeNumber){
    //return *(topologies.back());
    return *(topologies[moleculeNumber]);
}

/** Return a 2D vector representing all the coordinates of this World.
 * The first dimension represents the molecules (topologies) and the second
 * dimension (inner) represents the coordinates. The second inner dimension
 * type is pair of bSpecificAtom* and a Vec3. Thus, besides coordinates, it
 * contains all the information in bSpecificAtom as well. The bottleneck here
 * is the calcAtomLocationInGroundFrame from Compound.
 **/
std::vector < std::vector < std::pair <bSpecificAtom *, SimTK::Vec3> > >
        World::getAtomsLocationsInGround(const SimTK::State & state)
{
    std::vector< std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > >
            returnVector;

    // Iterate through topologies
    for ( unsigned int i = 0; i < this->topologies.size(); i++){
        std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>> currentTopologyInfo;
        // Iterate through atoms
        for(int j = 0; j < (topologies[i])->getNumAtoms(); j++){
            SimTK::Compound::AtomIndex compoundAtomIndex
                = topologies[i]->bAtomList[j].getCompoundAtomIndex();

            SimTK::Vec3 location = (topologies[i])
                    ->calcAtomLocationInGroundFrame(state, compoundAtomIndex);

            currentTopologyInfo.emplace_back(
                    &((topologies[i])->bAtomList[j]), location);
        }
        returnVector.emplace_back(currentTopologyInfo);
    }

    return returnVector;
}

/** Put coordinates into bAtomLists of Topologies.
 * When provided with a State, calcAtomLocationInGroundFrame
 * realizes Position and uses matter to calculate locations **/
void World::updateAtomListsFromCompound(const SimTK::State &state)
{
    // Iterate through topologies
    for ( unsigned int i = 0; i < this->topologies.size(); i++){
        // Iterate through atoms
        for(int j = 0; j < topologies[i]->getNumAtoms(); j++){
            SimTK::Compound::AtomIndex compoundAtomIndex = topologies[i]->bAtomList[j].getCompoundAtomIndex();
            SimTK::Vec3 location = topologies[i]->calcAtomLocationInGroundFrame(state, compoundAtomIndex);
            topologies[i]->bAtomList[j].setX(location[0]);
            topologies[i]->bAtomList[j].setY(location[1]);
            topologies[i]->bAtomList[j].setZ(location[2]);
        }
    }
}


// TODO Move in Topology
/** Set Compound, MultibodySystem and DuMM configurations according to
some other World's atoms **/ 
SimTK::State& World::setAtomsLocationsInGround(
        SimTK::State& someState,
        std::vector<
                std::vector<
                        std::pair<bSpecificAtom *, SimTK::Vec3> > > otherWorldsAtomsLocations)
{


    // Get the total no of bodies in this world (each World has its own
    // SimbodyMatterSubsystem)
    int totalNofBodies = matter->getNumBodies();

    // Declare arrays of Transforms that we need
    SimTK::Transform G_X_T;
    SimTK::Transform T_X_root[totalNofBodies]; // related to CompoundAtom.frameInMobilizedBodyFrame s

    // Iterate through molecules/topologies
    for(unsigned int i = 0; i < otherWorldsAtomsLocations.size(); i++) {

        // When in Debug mode
        if (visual == true) {
            paraMolecularDecorator->setAtomTargets(otherWorldsAtomsLocations[i]);
        }

        // Get the Ground to Top Transform
        SimTK::Transform G_X_T = topologies[i]->getTopLevelTransform();

        // Fully flexible regimen. realizeTopology is not needed
        if(topologies[i]->getRegimen().at(0) == 'I'){
            std::cout << "setAtomsLoc World IC" << std::endl << std::flush;

            // Create atomTargets
            std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;
            for(unsigned int j = 0; j < otherWorldsAtomsLocations[i].size(); j++){
                SimTK::Compound::AtomIndex atomIndex = ((otherWorldsAtomsLocations[i][j]).first)->atomIndex;
                SimTK::Vec3 location = ((otherWorldsAtomsLocations[i][j]).second);
                atomTargets.insert(std::pair<SimTK::Compound::AtomIndex, SimTK::Vec3>(atomIndex, location));
            }

            // As an alternative we can pool BAT from otherWorlds Qs
            // Match - only matchDefaultConfiguration should be necessary
            topologies[i]->matchDefaultAtomChirality(atomTargets, 0.01, false);
            topologies[i]->matchDefaultBondLengths(atomTargets);
            topologies[i]->matchDefaultBondAngles(atomTargets);
            topologies[i]->matchDefaultDirections(atomTargets);
            topologies[i]->matchDefaultDihedralAngles(atomTargets, SimTK::Compound::DistortPlanarBonds);
            topologies[i]->matchDefaultTopLevelTransform(atomTargets);
            // Seems to be redundant:
            //topologies[i]->matchDefaultConfiguration(atomTargets, SimTK::Compound::Match_Exact, true, 150.0);

            // Get transforms and locations: P_X_M, root_X_atom.p()
            G_X_T = topologies[i]->getTopLevelTransform();
            SimTK::Transform P_X_M[totalNofBodies]; // related to X_PFs

            // Iterate through atoms - get P_X_M for all the bodies
            for (SimTK::Compound::AtomIndex aIx(0); aIx < topologies[i]->getNumAtoms(); ++aIx){
                // Get body
                const auto mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);

                // Get P_X_M
                P_X_M[int(mbx)] = Transform(Rotation(), atomTargets[aIx]); // NEW
            }

            // Set X_FM and Q : backs up to stage Time
            for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
                SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);

                SimTK::Transform currentP_X_F = mobod.getInboardFrame(someState);
                SimTK::Transform invCurrentP_X_F = ~currentP_X_F;
                SimTK::Transform neededF_X_M = invCurrentP_X_F * P_X_M[int(mbx)];

                (mobod).setQToFitTransform(someState, neededF_X_M);

            }

        // END IC regimen

        //}else if((topologies[i]->getRegimen() == "TD") || (topologies[i]->getRegimen() == "RB")){
        }else{
            std::cout << "setAtomsLoc World TD" << std::endl << std::flush;
            // Create atomTargets
            std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;
            std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > currentTopology = otherWorldsAtomsLocations[i];
            for(unsigned int j = 0; j < currentTopology.size(); j++){
                SimTK::Compound::AtomIndex atomIndex = ((currentTopology[j]).first)->atomIndex;
                SimTK::Vec3 location = ((currentTopology[j]).second);
                atomTargets.insert(std::pair<SimTK::Compound::AtomIndex, SimTK::Vec3>(atomIndex, location));
            }

            // Match Default Configuration ( => default state)
            topologies[i]->matchDefaultAtomChirality(atomTargets, 0.01, false);
            topologies[i]->matchDefaultBondLengths(atomTargets);
            topologies[i]->matchDefaultBondAngles(atomTargets);
            topologies[i]->matchDefaultDirections(atomTargets);
            topologies[i]->matchDefaultDihedralAngles(atomTargets, SimTK::Compound::DistortPlanarBonds);
            topologies[i]->matchDefaultTopLevelTransform(atomTargets);
            topologies[i]->matchDefaultConfiguration(atomTargets, SimTK::Compound::Match_Exact, true, 150.0);

            // Get transforms and locations: P_X_F and root_X_atom.p()
            // M0 and Mr are actually used for F
            G_X_T = topologies[i]->getTopLevelTransform();
            SimTK::Transform invG_X_T;
            invG_X_T = ~G_X_T;
            SimTK::Transform M_X_pin = SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis); // Moves rotation from X to Z
            SimTK::Transform P_X_F[matter->getNumBodies()]; // related to X_PFs
            SimTK::Transform B_X_M[matter->getNumBodies()]; // related to X_BMs
            // SimTK::Transform T_X_root[matter->getNumBodies()]; // related to CompoundAtom.frameInMobilizedBodyFrame s
            SimTK::Transform T_X_Proot; // NEW
            SimTK::Transform root_X_M0[matter->getNumBodies()];
            SimTK::Angle inboardBondDihedralAngles[matter->getNumBodies()]; // related to X_FMs
            //SimTK::Real inboardBondLengths[matter->getNumBodies()]; // related to X_FMs
            SimTK::Vec3 locs[topologies[i]->getNumAtoms()];

            // Iterate through atoms - get T_X_roots for all the bodies
            for (SimTK::Compound::AtomIndex aIx(0); aIx < topologies[i]->getNumAtoms(); ++aIx){
                if(topologies[i]->getAtomLocationInMobilizedBodyFrame(aIx) == 0){ // atom is at body's origin
                    // Get body
                    SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                    T_X_root[int(mbx)] = topologies[i]->calcDefaultAtomFrameInCompoundFrame(aIx);
                }
            }

            P_X_F[1] = G_X_T * T_X_root[1]; // TODO: this doesn't work on multiple molecules
            // Iterate through atoms - get P_X_F for all the bodies
            for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
                if(topologies[i]->getAtomLocationInMobilizedBodyFrame(aIx) == 0){ // atom is at body's origin

                    // Get body, parentBody, parentAtom
                    SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
                    const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
                    SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

                    SimTK::Compound::AtomIndex parentAIx = topologies[i]->getMbx2aIx()[parentMbx];

                    if(parentMobod.getMobilizedBodyIndex() != 0){ // parent not Ground
                        // Find the true CHEMICAL parent
                        bSpecificAtom *originSpecAtom = nullptr;
                        originSpecAtom = (*(topologies[i])).updAtomByAtomIx(aIx);
                        //std::cout << "CHEM origin " << originSpecAtom->inName << std::endl;

                        SimTK::Compound::AtomIndex chemParentAIx;
                        for(unsigned int k = 0; k < (originSpecAtom->neighbors).size(); k++) {
                            //std::cout << "CHEM neighbours " << originSpecAtom->neighbors[k]->inName << std::endl;

                            for(std::vector<bBond *>::iterator bondsInvolvedIter = (originSpecAtom->bondsInvolved).begin();
                                bondsInvolvedIter != (originSpecAtom->bondsInvolved).end(); ++bondsInvolvedIter){
                                if((*bondsInvolvedIter)->isThisMe(originSpecAtom->getNumber(),
                                                                  originSpecAtom->neighbors[k]->getNumber())){
                                    Compound::AtomIndex candidateChemParentAIx = originSpecAtom->neighbors[k]->getCompoundAtomIndex();
                                    if(topologies[i]->getAtomMobilizedBodyIndex(candidateChemParentAIx) == parentMbx){
                                        if(!(*bondsInvolvedIter)->isRingClosing()){
                                            chemParentAIx = candidateChemParentAIx;
                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        Transform X_parentBC_childBC =
                                (*topologies[i]).getDefaultBondCenterFrameInOtherBondCenterFrame(aIx, chemParentAIx);
                        /////////////////////////////

                        T_X_Proot = T_X_root[parentMbx];

                        // Get inboard dihedral angle and put in root_X_M0 !!!!!!!
                        inboardBondDihedralAngles[int(mbx)] = topologies[i]->bgetDefaultInboardDihedralAngle(aIx);
                        //inboardBondLengths[int(mbx)] = topologies[i]->bgetDefaultInboardBondLength(aIx);
                        root_X_M0[int(mbx)] = SimTK::Transform(
                            SimTK::Rotation(inboardBondDihedralAngles[int(mbx)], SimTK::XAxis)
                            // , SimTK::Vec3(0, 0, 0) // Constructor takes care of this
                        );

                        // Get Proot_X_M0
                        SimTK::Transform T_X_M0 = T_X_root[int(mbx)] * root_X_M0[int(mbx)];
                        SimTK::Transform Proot_X_T = ~T_X_Proot;
                        SimTK::Transform Proot_X_M0 = Proot_X_T * T_X_M0;
                        //P_X_F[int(mbx)] = Proot_X_M0 * M_X_pin; // Get X_PF no chem

                        // Get X_PF, X_BM and X_FM with CHEMICAL joint placement
                        Transform oldX_PF = Proot_X_M0 * M_X_pin;
                        Transform oldX_BM = M_X_pin;
                        Transform oldX_MB = ~oldX_BM;
                        Transform oldX_FM = Rotation(inboardBondDihedralAngles[int(mbx)], ZAxis);

                        B_X_M[int(mbx)] = X_parentBC_childBC * M_X_pin;
                        P_X_F[int(mbx)] = oldX_PF * oldX_FM * oldX_MB * B_X_M[int(mbx)];
                        /////////////////////////////////

                    } //END if parent not Ground
                }
            }

            // Set transforms inside the bodies = root_X_atom.p; Set locations for everyone
            for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
                SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                if(topologies[i]->getAtomLocationInMobilizedBodyFrame(aIx) != 0){ // atom is not at body's origin
                    SimTK::Transform G_X_root = G_X_T * T_X_root[int(mbx)];
                    SimTK::Vec3 G_vchild = atomTargets[aIx];

                    SimTK::Transform root_X_child = alignFlipAndTranslateFrameAlongXAxis(G_X_root, G_vchild);
                    topologies[i]->bsetFrameInMobilizedBodyFrame(aIx, root_X_child);

                    locs[int(aIx)] = root_X_child.p();
                }
                else{
                    locs[int(aIx)] = SimTK::Vec3(0);
                }
            }

            // Set stations and AtomPLacements for atoms in DuMM
            for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
                SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                    SimTK::DuMM::AtomIndex dAIx = topologies[i]->getDuMMAtomIndex(aIx);

                    // Set station_B
                    forceField->bsetAtomStationOnBody( dAIx, locs[int(aIx)] );
                    forceField->bsetAllAtomStationOnBody( dAIx, locs[int(aIx)] ); // full

                    // Set included atom
                    //std::cout << "World setAtomLocations: updIncludedAtomStation(" << dAIx << ")" << std::endl;
                    forceField->updIncludedAtomStation(dAIx) = (locs[int(aIx)]);
                    forceField->updAllAtomStation(dAIx) = (locs[int(aIx)]); // full

                    // Atom placements in clusters
                    forceField->bsetAtomPlacementStation(dAIx, mbx, locs[int(aIx)] );
            }

            // Set X_PF and Q
            (topologies[i])->printMaps();
            for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx) {
                SimTK::MobilizedBody &mobod = matter->updMobilizedBody(mbx);

                // Get body, parentBody, parentAtom
                const SimTK::MobilizedBody &parentMobod = mobod.getParentMobilizedBody();
                SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
                SimTK::Compound::AtomIndex parentAIx = topologies[i]->getMbx2aIx()[parentMbx];
                SimTK::Compound::AtomIndex aIx = topologies[i]->getMbx2aIx()[mbx];
                std::cout << "mbx aIx parentAIx " << mbx << " " << aIx << " " << parentAIx << std::endl;

                if((!aIx.isValid()) && (!parentAIx.isValid())){
                    ;
                }else if((aIx.isValid()) && (!parentAIx.isValid())){ // Root atom
                    ((SimTK::MobilizedBody::Translation &) mobod).setDefaultInboardFrame(P_X_F[1]);
                    ((SimTK::MobilizedBody::Translation &) mobod).setDefaultOutboardFrame(Transform()); // NEWMOB
                }else{

                    int aNumber, parentNumber;
                    aNumber = (topologies[i]->updAtomByAtomIx(aIx))->getNumber();
                    parentNumber = (topologies[i]->updAtomByAtomIx(parentAIx))->getNumber();

                    SimTK::BondMobility::Mobility mobility;
                    mobility = (topologies[i]->getBond(aNumber, parentNumber)).getBondMobility();

                    //if(mobod.getNumU(someState) == 6){ // Free mobilizer
                    if (mobility == SimTK::BondMobility::Free) {
                        ((SimTK::MobilizedBody::Free &) mobod).setDefaultInboardFrame(P_X_F[int(mbx)]);
                        ((SimTK::MobilizedBody::Free &) mobod).setDefaultOutboardFrame(B_X_M[int(mbx)]); // NEWMOB

                        //}else if(mobod.getNumU(someState) == 0){ // Weld mobilizer
                    } else if (mobility == SimTK::BondMobility::Rigid) {
                        ((SimTK::MobilizedBody::Weld &) mobod).setDefaultInboardFrame(P_X_F[int(mbx)]);
                        ((SimTK::MobilizedBody::Weld &) mobod).setDefaultOutboardFrame(B_X_M[int(mbx)]); // CHEM

                        //}else if(mobod.getNumU(someState) == 1){ // Pin mobilizer
                    } else if (mobility == SimTK::BondMobility::Torsion) {
                        ((SimTK::MobilizedBody::Pin &) mobod).setDefaultInboardFrame(P_X_F[int(mbx)]);
                        ((SimTK::MobilizedBody::Pin &) mobod).setDefaultOutboardFrame(B_X_M[int(mbx)]); // CHEM
                        ((SimTK::MobilizedBody::Pin &) mobod).setDefaultQ(0); // CHEM
                        //((SimTK::MobilizedBody::Pin &) mobod).setDefaultQ(inboardBondDihedralAngles[int(mbx)]); // no chem

                        //}else if(mobod.getNumU(someState) == 2) { // Pin mobilizer
                    } else if (mobility == SimTK::BondMobility::Translation) {
                        ((SimTK::MobilizedBody::Translation &) mobod).setDefaultInboardFrame(P_X_F[int(mbx)]);
                        ((SimTK::MobilizedBody::Translation &) mobod).setDefaultOutboardFrame(B_X_M[int(mbx)]); // CHEM
                        ((SimTK::MobilizedBody::Translation &) mobod).setDefaultQ(SimTK::Vec3(0)); // CHEM

                    } else if (mobility == SimTK::BondMobility::Cylinder) {
                        ((SimTK::MobilizedBody::Cylinder &) mobod).setDefaultInboardFrame(P_X_F[int(mbx)]);
                        ((SimTK::MobilizedBody::Cylinder &) mobod).setDefaultOutboardFrame(B_X_M[int(mbx)]); // CHEM
                        ((SimTK::MobilizedBody::Cylinder &) mobod).setDefaultQ(Vec2(0, 0)); // CHEM
                        //((SimTK::MobilizedBody::Cylinder &) mobod).setDefaultQ( // no chem
                        // inboardBondDihedralAngles[int(mbx)], inboardBondLengths[int(mbx)]); // no chem

                        //} else if(mobod.getNumU(someState) == 3) { // Ball mobilizer
                    } else if (mobility == SimTK::BondMobility::UniversalM) {
                        ((SimTK::MobilizedBody::Universal &) mobod).setDefaultInboardFrame(P_X_F[int(mbx)]);
                        ((SimTK::MobilizedBody::Universal &) mobod).setDefaultOutboardFrame(B_X_M[int(mbx)]); // CHEM

                    } else if (mobility == SimTK::BondMobility::BallM) {
                        //std::cout << "mass props" << mobod.getBodyMassProperties(someState) << std::endl;
                        //if(mobod.getBodyMassProperties(someState).getUnitInertia().getMoments() == 0){ // Cartesian
                        //	std::cout << "atom in mobod " << (topologies[i]->getMbx2aIx())[mbx] << std::endl;;
                        //    //mobod.setQ();
                        //}else{
                        ((SimTK::MobilizedBody::Ball &) mobod).setDefaultInboardFrame(P_X_F[int(mbx)]);
                        ((SimTK::MobilizedBody::Ball &) mobod).setDefaultOutboardFrame(B_X_M[int(mbx)]); // CHEM

                        SimTK::Rotation R_FM;
                        R_FM.setRotationFromAngleAboutX(0.0);
                        R_FM.setRotationFromAngleAboutY(0.0);
                        //R_FM.setRotationFromAngleAboutZ(inboardBondDihedralAngles[int(mbx)]); // no chem
                        R_FM.setRotationFromAngleAboutZ(0); // CHEM
                        ((SimTK::MobilizedBody::Ball &) mobod).setDefaultRotation(R_FM);
                        //}
                    } else if (mobility == SimTK::BondMobility::Spherical) {
                        ((SimTK::MobilizedBody::SphericalCoords &) mobod).setDefaultInboardFrame(P_X_F[int(mbx)]);
                        ((SimTK::MobilizedBody::SphericalCoords &) mobod).setDefaultOutboardFrame(
                                B_X_M[int(mbx)]); // CHEM
                        ((SimTK::MobilizedBody::SphericalCoords &) mobod).setDefaultQ(SimTK::Vec3(0, 0, 0));
                    }
                }
            }

            // Set mass properties for mobilized bodies
            for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
                SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
                DuMM::ClusterIndex clusterIx = forceField->bgetMobodClusterIndex(mbx);
                SimTK::MassProperties massProperties = forceField->calcClusterMassProperties(clusterIx);
                mobod.setDefaultMassProperties(massProperties);
            }

            this->compoundSystem->realizeTopology();
            someState = compoundSystem->updDefaultState();

            // CHECK
/*            for(unsigned int j = 0; j < otherWorldsAtomsLocations[i].size(); j++){
                SimTK::Compound::AtomIndex aIx = ((otherWorldsAtomsLocations[i][j]).first)->atomIndex;
                SimTK::Vec3 location = ((otherWorldsAtomsLocations[i][j]).second);
                std::cout << "setAtomsLoc atomTargets from previous World i aIx "
                    << j << " " << aIx << " " << location << std::endl;
            }

            for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
                SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                SimTK::DuMM::AtomIndex dAIx = topologies[i]->getDuMMAtomIndex(aIx);

                // Check station_B
                std::cout << "setAtomsLoc aIx dumm.station_B " << aIx
                    << " " << forceField->getAtomStationOnBody(dAIx) << std::endl;

            }*/
            // CHECK END



        } // END TD regimen and all regimens

    } // END iterating through molecules/topologies

    this->compoundSystem->realize(someState, SimTK::Stage::Position);

    updateAtomListsFromCompound(someState);

    return someState;

}


/** Print information about Simbody systems. For debugging purpose. **/
void World::PrintSimbodyStateCache(SimTK::State& someState){
    std::cout << " System Stage: " << someState.getSystemStage() << std::endl;
    for(int i = 0; i < someState.getNumSubsystems(); i++){
        std::cout << " Subsystem " << i << " Name: " << someState.getSubsystemName(SimTK::SubsystemIndex(i))
            << " Stage: " << someState.getSubsystemStage(SimTK::SubsystemIndex(i))
            << " Version: " << someState.getSubsystemVersion(SimTK::SubsystemIndex(i)) << std::endl;
    }
}

/** Fill Worlds Cartesian coordinates buffers. 
To be called before use of getXs, getYs or getZs **/
void World::updateCoordBuffers(void)
{
    int allAtIx = -1;
    for(unsigned int tIx = 0; tIx < topologies.size(); tIx++){
        for(int aIx = 0; aIx < topologies[tIx]->getNAtoms(); aIx++){
            allAtIx++;
            Xs[allAtIx] = (topologies[tIx]->bAtomList[aIx]).getX();
            Ys[allAtIx] = (topologies[tIx]->bAtomList[aIx]).getY();
            Zs[allAtIx] = (topologies[tIx]->bAtomList[aIx]).getZ();
        }
    }

}

/** Get the coordinates from buffers **/
std::vector<SimTK::Real> World::getXs(void)
{
    return Xs;
}

std::vector<SimTK::Real> World::getYs(void)
{
    return Ys;
}

std::vector<SimTK::Real> World::getZs(void)
{
    return Zs;
}

/** Destructor **/
World::~World(){
    if(this->visual == true){
        //delete paraMolecularDecorator;

        delete decorations;
        delete visualizer;
        //delete vizReporter;
    }
    delete ts;
    delete integ;
    delete forceField;
    if(_useFixmanTorque){
        delete ExtForce;
    }
    delete matter;
    delete forces;
    delete compoundSystem;

    for(unsigned int i = 0; i < topologies.size(); i++){
        delete topologies[i];
    }
    for(unsigned int i = 0; i < samplers.size(); i++){
        delete samplers[i];
    }

    Xs.clear();
    Ys.clear();
    Zs.clear();
}

/** Return own CompoundSystem **/
CompoundSystem *World::getCompoundSystem() const {
    return compoundSystem;
}

/** Set own Compound system **/
// TODO find a solution for the old one
void World::setCompoundSystem(CompoundSystem *compoundSystem) {
    World::compoundSystem = compoundSystem;
}



