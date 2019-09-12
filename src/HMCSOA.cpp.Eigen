/**@file
Implementation of TestHMCSOA class. **/
#include "HMCSOA.hpp"

// Constructor
TestHMCSOA::TestHMCSOA(SimTK::CompoundSystem *argCompoundSystem,
                                     SimTK::SimbodyMatterSubsystem *argMatter,
                                     SimTK::Compound *argResidue,
                                     SimTK::DuMMForceFieldSubsystem *argDumm,
                                     SimTK::GeneralForceSubsystem *argForces,
                                     SimTK::TimeStepper *argTimeStepper)
    : MonteCarloSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
    , Sampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
{
    this->useFixman = true;  
    this->fix_n = this->fix_o = 0.0;
    this->residualEmbeddedPotential = 0.0;
    nofSamples = 0;

    //TO BE DELETED
    prevM = SimTK::Matrix(matter->getNumMobilities(), matter->getNumMobilities());
    //TO BE DELETED

}

// Destructor
TestHMCSOA::~TestHMCSOA()
{
}

// Compute O(n2) the square root of the mass matrix using using Eige - doesn't work !!n
void TestHMCSOA::calcNumSqrtMUpper(SimTK::State& someState, SimTK::Matrix& SqrtMUpper)
{
    int nu = someState.getNU();
    assert((SqrtMUpper.nrow() == nu) && (SqrtMUpper.ncol() == nu) && "calcSqrtMUpper: passed matrix doesn't have nu x nu size.");

    // Calc sqrt(MInv) and put it in SqrtM
    SimTK::Matrix SqrtMInvU(nu, nu);
    this->calcSqrtMInvU(someState, SqrtMInvU);

    // Put sqrt(MInv) in Eigen
    Eigen::MatrixXd EiSqrtMInvU(nu, nu);
    for(int i=0; i<nu; i++){
        for(int j=0; j<nu; j++){
            EiSqrtMInvU(i, j) = SqrtMInvU(i, j);
        }
    }
    std::cout << "SqrtMInvU: " << std::endl << SqrtMInvU << std::endl;
    std::cout << "EiSqrtMInvU: " << std::endl << EiSqrtMInvU << std::endl;

    // Compute the inverse of sqrt(M) = inv(sqrt(MInv)) with Eigen
    Eigen::MatrixXd EiSqrtMUpper(nu, nu);
    EiSqrtMUpper = EiSqrtMInvU.inverse();

    // Put sqrt(M) back in Simbody style matrix SqrtM
    for(int i=0; i<nu; i++){
        for(int j=0; j<nu; j++){
            SqrtMUpper(i, j) = EiSqrtMUpper(i, j);
        }
    }
    std::cout << "EiSqrtMUpper: " << std::endl << EiSqrtMUpper << std::endl;
    std::cout << "SqrtMUpper: " << std::endl << SqrtMUpper << std::endl;

    // ---- Verify with Eigen ----------
    // Get M
    SimTK::Matrix M(nu, nu);
    matter->calcM(someState, M);
    //std::cout << "M: " << M << std::endl;
    SimTK::Matrix MInv(nu, nu);
    matter->calcMInv(someState, MInv);
    std::cout << "MInv: " << MInv << std::endl;

    // Transfer into Eigen
    Eigen::MatrixXd EiVer(nu, nu);
    for(int i=0; i<nu; i++){
        for(int j=0; j<nu; j++){
            EiVer(i, j) = SqrtMUpper(i, j);
        }
    }
    std::cout << "EiMInv: " << EiSqrtMInvU.transpose() * EiSqrtMInvU << std::endl;
    //std::cout << "EiM: " << EiVer.transpose() * EiVer << std::endl;
    
}

// Calculate O(n2) the square root of the mass matrix inverse
void TestHMCSOA::calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv)
{
    int nu = someState.getNU();
    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvU: passed matrix doesn't have nu x nu size.");

    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);

    // Zero the matrix and the temporary vector
    for (int i=0; i < nu; ++i){
        V[i] = 0;
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = 0;
        }
    }

    // Calculate the values inside the matrix
    for (int i=0; i < nu; ++i){
        V[i] = 1;
        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = SqrtMInvV[j];
        }
        V[i] = 0;
    }

    // ---- Verify with Eigen ----------
    // Get M
    //SimTK::Matrix MInv(nu, nu);
    //matter->calcMInv(someState, MInv);
    //std::cout << "MInv: " << MInv << std::endl;

    // Transfer into Eigen
    //Eigen::MatrixXd EiMInv(nu, nu);
    //Eigen::MatrixXd EiSqrtMInv(nu, nu);
    //for(int i=0; i<nu; i++){
    //    for(int j=0; j<nu; j++){
    //        EiSqrtMInv(i, j) = SqrtMInv(i, j);
    //    }
    //}
    //std::cout << "Eigen MInv calc: " << EiSqrtMInv.transpose() * EiSqrtMInv << std::endl;

    // Diagonalization
    ////Eigen::EigenSolver<Eigen::MatrixXd> EiSoM(EiM);
    ////Eigen::MatrixXcd EiD = EiSoM.eigenvalues().asDiagonal();
    ////Eigen::MatrixXcd EiV = EiSoM.eigenvectors();

}

void TestHMCSOA::calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv)
{
    int nu = someState.getNU();
    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvL: passed matrix doesn't have nu x nu size.");

    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);

    // Zero the matrix and the temporary vector
    for (int i=0; i < nu; ++i){
        V[i] = 0;
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = 0;
        }
    }

    // Calculate the values inside the matrix
    for (int i=0; i < nu; ++i){
        V[i] = 1;
        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
        for (int j=0; j < nu; ++j){
            SqrtMInv(j, i) = SqrtMInvV[j];
        }
        V[i] = 0;
    }

    // ---- Verify with Eigen ----------
    // Get M
    //SimTK::Matrix MInv(nu, nu);
    //matter->calcMInv(someState, MInv);
    //std::cout << "MInv: " << MInv << std::endl;

    // Transfer into Eigen
    //Eigen::MatrixXd EiMInv(nu, nu);
    //Eigen::MatrixXd EiSqrtMInv(nu, nu);
    //for(int i=0; i<nu; i++){
    //    for(int j=0; j<nu; j++){
    //        EiSqrtMInv(i, j) = SqrtMInv(i, j);
    //    }
    //}
    //std::cout << "Eigen MInv calc: " << EiSqrtMInv * EiSqrtMInv.transpose() << std::endl;

}

// Set set kinetic energy
void TestHMCSOA::setSetKE(SimTK::Real inpKE)
{
    this->ke_set = inpKE;
}

// Set old kinetic energy
void TestHMCSOA::setOldKE(SimTK::Real inpKE)
{
    this->ke_o = inpKE;
}

// Initialize variables (identical to setTVector)
void TestHMCSOA::initialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature, bool argUseFixman)
{
    randomEngine.seed( std::time(0) );
    //compoundSystem->realizeTopology(); // ELIMINATE REALIZE TOPOLOGY
    //SimTK::State state = compoundSystem->updDefaultState();
    timeStepper->initialize(compoundSystem->getDefaultState());
    setTemperature(argTemperature); // Needed for Fixman

    // Initialize x
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }
    setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    this->useFixman = argUseFixman;  
    if(useFixman){
        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());
    }

    // Initialize velocities to temperature
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int j=0; j < nu; ++j){
        V[j] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;
    system->realize(someState, SimTK::Stage::Velocity);
    setOldKE(matter->calcKineticEnergy(someState));
    setSetKE(getOldKE());
    this->etot_o = getOldPE() + getOldKE() + getOldFixman();
    this->etot_set = this->etot_o;
 
    // TO BE DELETED 
    prevTheta = SimTK::Vector(someState.getNQ());
    // TO BE DELETED 

    //timeStepper->initialize(someState);
    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Velocity, false);

    #ifdef HARMONICOSCILLATOR
    /////////////////////////
    // Harmonic oscillator
    /////////////////////////
    HO_randomEngine.seed( std::time(0) );
    // Initialize x
    for(int i = 0; i < HO_D; i++){HO_x[i] = HO_xini[i] = 1.0;}
    HO_PE_set = HO_PE_xprop = HO_PE_x = HarmonicOscillatorPE(HO_x);
    // Initialize v
    HO_InitializeVelocity(HO_v, getTemperature());
    HO_KE_set = HO_KE_xprop = HO_KE_x = HarmonicOscillatorKE(HO_v);
    HO_etot_x = HO_PE_x + HO_KE_x;
    HO_etot_set = HO_etot_xprop = HO_etot_x;
    #endif
}

// Initialize variables (identical to setTVector)
void TestHMCSOA::reinitialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature)
{
    setTemperature(argTemperature); // Needed for Fixman

    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Topology, false);

    // Initialize x
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }
    setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    if(useFixman){
        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());
    }

    // Initialize velocities to temperature
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int j=0; j < nu; ++j){
        V[j] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;
    system->realize(someState, SimTK::Stage::Velocity);
    setOldKE(matter->calcKineticEnergy(someState));
    setSetKE(getOldKE());
    this->etot_o = getOldPE() + getOldKE() + getOldFixman();
    this->etot_set = this->etot_o;
  
    //std::cout <<  "Sampler after reinitialize State Cache Info: " << std::endl;
    //PrintSimbodyStateCache(someState);

    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Velocity, false);

    #ifdef HARMONICOSCILLATOR
    /////////////////////////
    // Harmonic oscillator
    /////////////////////////
    // Initialize x
    for(int i = 0; i < HO_D; i++){HO_x[i] = HO_xini[i] = 1.0;}
    HO_PE_set = HO_PE_xprop = HO_PE_x = HarmonicOscillatorPE(HO_x);
    // Initialize v
    HO_InitializeVelocity(HO_v, getTemperature());
    HO_KE_set = HO_KE_xprop = HO_KE_x = HarmonicOscillatorKE(HO_v);
    HO_etot_x = HO_PE_x + HO_KE_x;
    HO_etot_set = HO_etot_xprop = HO_etot_x;
    #endif
}


// Assign random conformation
// In torsional dynamics the first body has 7 Q variables for 6 dofs - one
// quaternion (q) and 3 Cartesian coordinates (x). updQ will return: 
// [qw, qx, qy, qz, x1, x2, x3]
void TestHMCSOA::propose(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
    randomEngine.seed(4294653137UL); // for reproductibility

    system->realize(someState, SimTK::Stage::Position);
    // Initialize x - not necessary
    int t = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        TVector[t] = SetTVector[t];

        // A bit of study. TO BE DELETED
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        //std::cout << "Mobod " << t+1 << " :" << std::endl;
        for(int k = 0; k < mobod.getNumU(someState); k++){
            SimTK::SpatialVec HCol = mobod.getHCol(someState, SimTK::MobilizerUIndex(k));
        }
        // A bit of study end. TO BE DELETED

        t++;
    }
    setOldPE(getSetPE());
    setOldFixman(getSetFixman());

    // Assign velocities according to Maxwell-Boltzmann distribution
    // and set Old kinetic energy
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int i=0; i < nu; ++i){
        V[i] = gaurand(randomEngine);
    }
    //std::cout << "Q: " << someState.getQ() << std::endl;
    //std::cout << "Before stepTo PE: " << forces->getMultibodySystem().calcPotentialEnergy(someState) << std::endl;
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    //std::cout << "TestHMCSOA::propose SqrtMInvV: " << SqrtMInvV << std::endl;

    //SimTK::Real temperatureBoost = 2.58; // sqrt(2000/300) : brings temperature from 300 to 1000
    //SqrtMInvV *= sqrtRT * temperatureBoost; // Set stddev according to temperature
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature

    someState.updU() = SqrtMInvV;

    // TO BE DELETED
    //someState.updU() = 0.0;
    // TO BE DELETED

    // Set old kinetic energy
    system->realize(someState, SimTK::Stage::Velocity);
    //std::cout << "FixmanTorque 1: " << "Qs = " << someState.getQ() << std::endl;
    setOldKE(matter->calcKineticEnergy(someState));
    //std::cout << "FixmanTorque 2: " << "Qs = " << someState.getQ() << std::endl;
    this->etot_o = getOldPE() + getOldKE() + getOldFixman();


    // Propagate through phase space (integrate)
    //std::cout << "Before stepTo time: " << someState.getTime() << std::endl;
    //std::cout << "Q = " << someState.getQ() << std::endl;
    this->timeStepper->stepTo(someState.getTime() + (timestep*nosteps));
    //std::cout << "Q = " << someState.getQ() << std::endl;
    //std::cout << "After stepTo time: " << someState.getTime() << std::endl;

    // Compute Fixman torque
    //SimTK::Vector V3(nu);
    //SimTK::Vector V4(nu);
    //SimTK::Real* D0 = new SimTK::Real(1.0);
    //matter->calcFixmanTorque(someState, V3, V4, D0);
    //std::cout << "Fixman torque " ;
    //for(int i = 0; i < nu; i++){
    //    std::cout << std::setprecision(10) << V4[i] << " ";
    //}
    //std::cout << std::endl;
    //delete D0;
    // end - Compute Fixman torque
 
    //PrintSimbodyStateCache(someState);
    //std::cout << "After  stepTo time: " << someState.getTime() << std::endl;
    //writePdb(*residue, someState, "pdbs", "sb_", 8, "HMCprop", nofSamples);

    //std::cout << "After  stepTo Q: " << someState.getQ() << std::endl;
    //std::cout << "After  stepTo U: " << someState.getU() << std::endl;
    //std::cout << "After  stepTo PE: " << forces->getMultibodySystem().calcPotentialEnergy(someState) << std::endl;

    // Get M
    //matter->calcM(someState, M);
    //std::cout << "After stepTo:" << std::setprecision(3) << std::endl;
    //for(int i=0; i<nu; i++){
    //    std::cout << "M: ";
    //    for(int j=0; j<nu; j++){
    //        std::cout << M(i, j) << " ";
    //    }
    //    std::cout << std::endl;
    //}

    #ifdef HARMONICOSCILLATOR
    /////////////////////////
    // Harmonic oscillator
    /////////////////////////

    for(int i = 0; i < HO_D; i++){HO_x[i] = HO_xini[i];}
    HO_PE_x = HO_PE_set;
    HO_InitializeVelocity(HO_v, getTemperature());
    HO_KE_x = HarmonicOscillatorKE(HO_v);
    HO_etot_x = HO_PE_x + HO_KE_x;

    HO_VelocityVerlet(timestep, 100);
    HO_PE_xprop = HarmonicOscillatorPE(HO_xprop);
    HO_KE_xprop = HarmonicOscillatorKE(HO_v);
    HO_etot_xprop = HO_PE_xprop + HO_KE_xprop;
    //////////////////////////
    #endif

}

// The update step in Monte Carlo methods consists in:
// Acception - rejection step
void TestHMCSOA::update(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);

    propose(someState, timestep, nosteps);

    SimTK::Real pe_o  = getOldPE();
    if(useFixman){
        SimTK::Real fix_o = getOldFixman();
    }
    SimTK::Real ke_o  = getOldKE();

    if(useFixman){
        fix_n = calcFixman(someState);
    }else{
        fix_n = 0.0;
    }
    SimTK::Real pe_n = getPEFromEvaluator(someState); // OPENMM
    system->realize(someState, SimTK::Stage::Velocity);
    SimTK::Real ke_n = matter->calcKineticEnergy(someState);

    // Apply Metropolis criterion
    SimTK::Real etot_o, etot_n;
    assert(!std::isnan(pe_n));
    if(useFixman){
        etot_n = pe_n + ke_n + fix_n;
        etot_o = pe_o + ke_o + fix_o;
    }else{
        etot_n = pe_n + ke_n;
        etot_o = pe_o + ke_o;
    }

    etot_o;
    etot_n;

    //std::cout <<  "Sampler after energies calculations State Cache Info: " << std::endl;
    //PrintSimbodyStateCache(someState);

    std::cout << "UDot = " << someState.getUDot() << std::endl;
    std::cout<<std::setprecision(5)<<std::fixed;
    std::cout << "pe_o " << pe_o + getREP() << " ke_o " << ke_o << " fix_o " << fix_o << " rep " << getREP()
        << " pe_n " << pe_n  + getREP() << " ke_n " << ke_n << " fix_n " << fix_n
        //<< " rand_no " << rand_no << " RT " << RT << " exp(-(etot_n - etot_o) " << exp(-(etot_n - etot_o) / RT)
        //<< " etot_n " << etot_n  + getREP() << " etot_o " << etot_o + getREP()
        ;

    //if (( (pe_n + fix_n) < (pe_o + fix_o) ) || (rand_no < exp(-( (pe_n + fix_n) - (pe_o + fix_o) ) / RT))){ // Accept
    //if (( (pe_n) < (pe_o) ) || (rand_no < exp(-( (pe_n) - (pe_o) ) / RT))){ // Accept
    //if ((etot_n < etot_o) || (rand_no < exp(-(etot_n - etot_o)/RT))){ // Accept

    if(1){ // Always accept
    //if ((etot_n < etot_o) || (rand_no < exp(-(etot_n - etot_o)/RT))){ // Accept
    //if ( (!std::isnan(pe_n)) && (etot_n < etot_o) || (rand_no < exp(-(etot_n - etot_o)/RT)) ){ // Accept
        std::cout << " acc 1 " ;
        setSetTVector(someState);
        //sendConfToEvaluator(); // OPENMM
        setSetPE(pe_n);
        setSetFixman(fix_n);
        setSetKE(ke_n);
        this->etot_set = getSetPE() + getSetFixman() + getSetKE();
        //someState.updU() = 0.0;
        //setOldKE(0.0);

    }else{ // Reject
        std::cout << " acc 0 " ;
        assignConfFromSetTVector(someState);
        //someState.updU() = 0.0;
        //setOldKE(0.0);

    }

    std::cout << " pe_os " << getSetPE() + getREP() << " ke_os " << getSetKE() << " fix_os " << getSetFixman()
        << " pe_n " << pe_n << " ke_n " << ke_n << " fix_n " << fix_n
        << std:: endl;
    //std::cout << "Number of times the force field was evaluated: " << dumm->getForceEvaluationCount() << std::endl;

    ++nofSamples;


    // SOA
    /////////////////////////////////////////////
    // Get Fixman torque numerically. [Jain 1997] Eq. 2.11 and 2.14.
    // TODO: detM quatities are actually log(det(M))
    /////////////////////////////////////////////


    int NBODIESplusG = matter->getNumBodies();
    int SPATIAL_DOFSplusG = NBODIESplusG * 6;
    int SPATIAL_DOFS = SPATIAL_DOFSplusG - 6;
    int NU = someState.getNU();
    SimTK::Matrix M(NU, NU);
    SimTK::Matrix MInv(NU, NU);
    SimTK::Matrix dM(NU, NU);
    SimTK::Matrix dMdThetaK(NU, NU);
    SimTK::Matrix MInvdMdThetaK(NU, NU);
    SimTK::Vector dTheta(NU);
    SimTK::Real detM;
    SimTK::Real dDetM;
    SimTK::Real dDetMdThetaK;
    SimTK::Real numDetM;
    SimTK::Real dNumDetM;
    SimTK::Real dNumDetMdThetaK;
    SimTK::Real ThetaK;
    SimTK::Real dThetaK;
    SimTK::Real MInvdMdThetaKTrace = 0.0;



    // Compute quatities at time t
    SimTK::Vector V1(NU);
    SimTK::Vector V2(NU);
    SimTK::Real* D0 = new SimTK::Real(1.0);
    matter->calcDetM(someState, V1, V2, D0);
    
    std::cout << "FixmanTorque 8: " << std::setprecision(10) << "Qs    = " << someState.getQ() << std::endl;
    std::cout << "FixmanTorque  : " << std::setprecision(10) << "QDots = " << someState.getQDot() << std::endl;
    std::cout << "FixmanTorque  : " << std::setprecision(10) << "U     = " << someState.getU() << std::endl;

    detM = std::log(*D0);
    delete D0;
    numDetM = std::log(calcNumDetM(someState));

    matter->calcM(someState, M);
    matter->calcMInv(someState, MInv);

    // Compute some differentials
    dM = M - prevM;
    dDetM = detM - prevDetM;
    dNumDetM = numDetM - prevNumDetM;
    //dTheta = someState.getQ() - prevTheta;
    dTheta = someState.getQDot();
    //PrintBigMat(MInv, NU, NU, 2, "MInv"); 

    //for (int k = 0; k < someState.getNQ(); k++){
    for (int k = 0; k < someState.getNU(); k++){
        std::cout << "Numerical Fixman Torque for U " << k << std::endl;
    
        // Compute derivatives
        //dDetMdThetaK = dDetM / dTheta[k];
        //dNumDetMdThetaK = dNumDetM / dTheta[k];
        //dMdThetaK = dM / dTheta[k];
        dDetMdThetaK = dDetM / someState.getU()[k];
        dNumDetMdThetaK = dNumDetM / someState.getU()[k];
        dMdThetaK = dM / someState.getU()[k];

        //PrintBigMat(M, NU, NU, 2, "M"); 
        //PrintBigMat(dM, NU, NU, 4, "dM");
        //std::cout << std::setprecision(10) << "dThetaK " << dTheta[int(mbx)] << std::endl << std::setprecision(2);
    
        MInvdMdThetaK = MInv * dMdThetaK;
        for(int i = 0; i < NU; i++){
            MInvdMdThetaKTrace += MInvdMdThetaK(i, i);
        }
        //PrintBigMat(MInvdMdThetaK, NU, NU, 4, "MInvdMdThetaK");
        std::cout << "HMC FixmanTorque: " << std::setprecision(10) << "1/2 dNumDetMdThetaK = " << 0.5 * dNumDetMdThetaK << std::endl;
        std::cout << "HMC FixmanTorque: " << "1/2 dDetMdThetaK = " << 0.5 * dDetMdThetaK << std::endl;
        std::cout << "HMC FixmanTorque: 1/2 Trace(MInvdMdThetaK):" << 0.5 * MInvdMdThetaKTrace << std::endl << std::setprecision(2);
    }

    // Save t-1 quantities
    for(int i = 0; i < NU; i++){
        for(int j = 0; j < NU; j++){prevM(i, j) = M(i, j);}
    }
    prevDetM = detM;
    prevNumDetM = numDetM;
    prevTheta = someState.getQ();

    /////////////////////////////////////////////
    // Get Newton Euler Operators H, Phi, and MkTot. [Rodriguez 1991] Eq. 3.2
    /////////////////////////////////////////////
    //Get Jacobian J = Phi H
    SimTK::Matrix Jstar;
    SimTK::Matrix J;
    SimTK::Matrix MkTot(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    matter->calcSystemJacobian(someState, J);
    Jstar = ~J;

    // Get MkTot
    SimTK::Matrix M3_4a(NU, NU, 0.0);
    for (SimTK::MobilizedBodyIndex mbx(0); mbx < NBODIESplusG; ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::SpatialInertia Mk = mobod.getBodySpatialInertiaInGround(someState);
        //std::cout << "FixmanTorque matrix: " << "Mk = " << Mk.toSpatialMat() << std::endl;
        int ti = -1; // 6x6
        int tj = -1; // 6x6
        for(int i = 0; i < 2; i++){ // i for SpatialMatrix
            for(int k = 0; k < 3; k++){ // k for 3x3 
                for(int j = 0; j < 2; j++){ // j for SpatialMatrix
                    for(int l = 0; l < 3; l++){ // l for 3x3
                        tj++;
                        tj = tj % 6;
                        if(!tj){ti++;}
                        if(std::isinf(Mk.toSpatialMat()[i][j][k][l]) || std::isnan(Mk.toSpatialMat()[i][j][k][l])){
                            MkTot[int(mbx) * 6 + ti][int(mbx) * 6 + tj] = 0.0;
                        }else{
                            MkTot[int(mbx) * 6 + ti][int(mbx) * 6 + tj] = Mk.toSpatialMat()[i][j][k][l];
                        }
                    }
                }
            }
        }
    }
    //PrintBigMat(MkTot, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "MkTot");

    // Get H based on H_PB_G (what Jain calls Hstar) 
    SimTK::Matrix H(NU, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix Hstar(SPATIAL_DOFSplusG, NU, 0.0);
    SimTK::Matrix H_FM(NU, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix H_FMstar(SPATIAL_DOFSplusG, NU, 0.0);
    SimTK::Matrix Phi(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix Phistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);

    int tj = -1;
    for (SimTK::MobilizedBodyIndex mbx(0); mbx < NBODIESplusG; ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

        for(int k = 0; k < mobod.getNumU(someState); k++){
            tj++;
            SimTK::SpatialVec HCol = mobod.getHCol(someState, SimTK::MobilizerUIndex(k));
            SimTK::SpatialVec H_FMCol = mobod.getH_FMCol(someState, SimTK::MobilizerUIndex(k));
            //std::cout << "mobod " << int(mbx) << " HCol " << k << " = " << HCol << std::endl;
            int ti = -1;
            for(int i = 0; i < 2; i++){
                for(int j = 0; j < 3; j++){
                    ti++;
                    Hstar[int(mbx)*6 + ti][tj] = HCol[i][j];
                    H_FMstar[int(mbx)*6 + ti][tj] = H_FMCol[i][j];
                }
            }
        }
        
    }

    H = Hstar.transpose();
    H_FM = H_FMstar.transpose();
    PrintBigMat(Hstar, SPATIAL_DOFSplusG, NU, 2, "Hstar"); 
    //PrintBigMat(H_FMstar, SPATIAL_DOFSplusG, NU, 2, "H_FMstar"); 

    // Get Phi
    for(int phi_I = 0; phi_I < NBODIESplusG; phi_I++){
        // Fill diagonal element
        for(int i = 0; i < 6; i++){
                Phistar[phi_I*6 + i][phi_I*6 + i] = 1.0;
        }
        // Fill lower triangle
        for(int phi_J = 0; (phi_J < NBODIESplusG) && (phi_J < phi_I); phi_J++){

            const SimTK::MobilizedBody& Jmobod = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(phi_J));
            const SimTK::MobilizedBody& Imobod = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(phi_I));
            //const SimTK::MobilizedBody& Imobod = Jmobod.getParentMobilizedBody();

            const SimTK::Transform X_GI = Imobod.getBodyTransform(someState);
            const SimTK::Transform X_GJ = Jmobod.getBodyTransform(someState);

            const SimTK::Transform X_IBM = Imobod.getOutboardFrame(someState);
            const SimTK::Transform X_JBM = Jmobod.getOutboardFrame(someState);

            const SimTK::Transform X_IFM = Imobod.getMobilizerTransform(someState);
            const SimTK::Transform X_JFM = Jmobod.getMobilizerTransform(someState);

            SimTK::Mat33 PhiLCross = SimTK::crossMat( X_GJ.p() - X_GI.p() );

            SimTK::Matrix PhiElem(6, 6, 0.0);
            for(int i = 0; i < 6; i++){
                    PhiElem[i][i] = 1.0;
            }
            for(int i = 3; i < 6; i++){
                for(int j = 0; j < 3; j++){
                    PhiElem[i][j] = PhiLCross[i - 3][j];
                }
            }

            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    Phistar[phi_I*6 + i][phi_J*6 + j] = PhiElem[i][j];
                    //Phistar[Imobod.getMobilizedBodyIndex()*6 + i][Jmobod.getLevelInMultibodyTree()*6 + j] = PhiElem[i][j];
                }
            }
    
            //std::cout << "phi tree levels " << Imobod.getLevelInMultibodyTree() << " " << Jmobod.getLevelInMultibodyTree() << std::endl;
            //PrintBigMat(PhiElem, 6, 6, 3, "PhiElem");
        }
    }
    Phi = Phistar.transpose();
    //PrintBigMat(Phistar, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "FixmanTorque matrix: Phistar"); 

    // Verify [Jain 1997] Eq. 3.4a
    //M3_4a = (~J) * MkTot * J;
    //M3_4a = (~Hstar * ~Phistar) * MkTot * (Phistar * Hstar);
    //std::cout << "M - M3_4a = " << M - M3_4a << std::endl;

    /////////////////////////////////////////////
    // Get operators Groundless truncated versions (Jain's versions ?)
    /////////////////////////////////////////////
    SimTK::Matrix J_Jain(SPATIAL_DOFS, NU, 0.0);
    SimTK::Matrix Phistar_Jain(SPATIAL_DOFS, SPATIAL_DOFS, 1.0);
    SimTK::Matrix Hstar_Jain(SPATIAL_DOFS, NU, 0.0);
    SimTK::Matrix Jstar_Jain(NU, SPATIAL_DOFS, 0.0);
    SimTK::Matrix Phi_Jain(SPATIAL_DOFS, SPATIAL_DOFS, 1.0);
    SimTK::Matrix H_Jain(NU, SPATIAL_DOFS, 0.0);
    SimTK::Matrix MkTot_Jain(SPATIAL_DOFS, SPATIAL_DOFS, 0.0);

    for(int i = 0; i < SPATIAL_DOFS; i++){
        for(int j = 0; j < NU; j++){
            Hstar_Jain[i][j] = Hstar[i + 6][j];
            J_Jain[i][j] = J[i + 6][j];
        }
    }
    for(int i = 0; i < SPATIAL_DOFS; i++){
        for(int j = 0; j < SPATIAL_DOFS; j++){
            MkTot_Jain[i][j] = MkTot[i + 6][j + 6];
            Phistar_Jain[i][j] = Phistar[i + 6][j + 6];
        }
    }

    Phi_Jain = Phistar_Jain.transpose();
    H_Jain = Hstar_Jain.transpose();
    Jstar_Jain = H_Jain.transpose();

    /////////////////////////////////////////////
    // Get SOA Operators [Rodriguez 1991, 1992]
    /////////////////////////////////////////////
    SimTK::Matrix P(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix P_Jain(SPATIAL_DOFS, SPATIAL_DOFS, 0.0);
    SimTK::Matrix D(NU, NU, 0.0);
    SimTK::Matrix DInv(NU, NU, 0.0);
    SimTK::Matrix D_Jain(NU, NU, 0.0);
    SimTK::Matrix G(SPATIAL_DOFSplusG, NU, 0.0);
    SimTK::Matrix G_Jain(SPATIAL_DOFS, NU, 0.0);
   
    // Get P, D, and G 
    for (SimTK::MobilizedBodyIndex mbx(0); mbx < NBODIESplusG; ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::ArticulatedInertia Pk = matter->getArticulatedBodyInertia(someState, mbx);
        SimTK::Mat<6, 6> PkMat66;
        SpatialMat2Mat66(Pk.toSpatialMat(), PkMat66);
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                P[int(mbx) * 6 + i][int(mbx) * 6 + j] = PkMat66[i][j];
            }
        }
    }
    for(int i = 0; i < SPATIAL_DOFS; i++){
        for(int j = 0; j < SPATIAL_DOFS; j++){
            P_Jain[i][j] = P[i + 6][j + 6];
        }
    }

    D = H * P * Hstar;
    D_Jain = H_Jain * P_Jain * Hstar_Jain;
    NumericalInverse(D, DInv, NU, NU);
    G = P * Hstar * DInv;
    G_Jain = P_Jain * Hstar_Jain * DInv;
    
    // Get Eta operators and K
    SimTK::Matrix EtaPhistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix EtaPhistar_Jain(SPATIAL_DOFS, SPATIAL_DOFS, 0.0);
    SimTK::Matrix EtaPhi(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix EtaPhi_Jain(SPATIAL_DOFS, SPATIAL_DOFS, 0.0);

    SimTK::Matrix K(SPATIAL_DOFSplusG, NU, 0.0);
    SimTK::Matrix K_Jain(SPATIAL_DOFS, NU, 0.0);

    for(int i = 1; i < NBODIESplusG; i++){
        for(int k = 0; k < 6; k++){
            for(int l = 0; l < 6; l++){
                EtaPhistar[i*6 + k][(i*6 - 6) + l] = Phistar[i*6 + k][(i*6 - 6) + l];
            }
        }
    }
    for(int i = 0; i < SPATIAL_DOFS; i++){
        for(int j = 0; j < SPATIAL_DOFS; j++){
            EtaPhistar_Jain[i][j] = EtaPhistar[i + 6][j + 6];
        }
    }

    EtaPhi = EtaPhistar.transpose();
    EtaPhi_Jain = EtaPhistar_Jain.transpose();
    K = EtaPhi * G;
    K_Jain = EtaPhi_Jain * G_Jain;
    //PrintBigMat(EtaPhistar, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "EtaPhistar"); 
    //PrintBigMat(K, SPATIAL_DOFSplusG, NU, 2, "K"); 

    // Verify SOA using the mass matrix
    //SimTK::Matrix I(NU, NU, 0.0);
    //for(int i = 0; i < NU; i++){
    //    I[i][i] = 1.0;
    //}
    //SimTK::Matrix MSqParan(NU, NU, 0.0);
    //MSqParan = I + (H * Phi * K);
    //PrintBigMat( (H * Phi) * MkTot * (Phistar * Hstar), NU, NU, 2, "Mcheck 3.4a" ); // Eq 3.4a
    //PrintBigMat( MSqParan * D * MSqParan.transpose(), NU, NU, 2, "Mcheck 3.5a 1");

    // Check EtaPhistar 1992 Rodriguez eq. 0.6
    //SimTK::Matrix I(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    //SimTK::Matrix I_EtaPhistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    //SimTK::Matrix I_EtaPhistarInv(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);

    //for(int i = 0; i < SPATIAL_DOFSplusG; i++){
    //    I[i][i] = 1.0;
    //}
    //I_EtaPhistar = I - EtaPhistar;
    //NumericalInverse(I_EtaPhistar, I_EtaPhistarInv, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG);
    //PrintBigMat(I_EtaPhistar, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "I_EtaPhistar"); 
    //PrintBigMat(Phistar - I_EtaPhistarInv, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Check EtaPhistar"); 

    // Get Psi and EtaPsi
    SimTK::Matrix Psi(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix Psistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix EtaPsi(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix EtaPsistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);

    EtaPsi = EtaPhi - (K * H);
    EtaPsistar = EtaPsi.transpose();

    // Fill Psi
    Psi = EtaPsi;
    for(int i = 0; i < NBODIESplusG; i++){
        for(int k = 0; k < 6; k++){
            Psi[i*6 + k][i*6 + k] = 1.0;
        }
        for(int j = i + 2; j < NBODIESplusG; j++){
            SimTK::Mat66 Psi_IJ, Psi_IJm, Psi_JmJ;
            for(int k = 0; k < 6; k++){
                for(int l = 0; l < 6; l++){
                    Psi_IJ[ k][l] = Psi[  (i*6)      + k ][  (j*6)      + l];
                    Psi_IJm[k][l] = Psi[  (i*6)      + k ][ ((j*6) - 6) + l];
                    Psi_JmJ[k][l] = Psi[ ((j*6) - 6) + k ][  (j*6)      + l];
                }
            }

            Psi_IJ = Psi_IJm * Psi_JmJ;

            for(int k = 0; k < 6; k++){
                for(int l = 0; l < 6; l++){
                    Psi[ (i*6) + k ][ (j*6) + l] = Psi_IJ[k][l];
                }
            }
            
        }
    }
    Psistar = Psi.transpose();

    /////////////////////////////////////////////
    // Get Omega using [Jain 1997] Eq. 4.13
    /////////////////////////////////////////////
    SimTK::Matrix Omega(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    Omega = Psistar * Hstar * DInv * H * Psi;

    // Check eq. 4.14b
    //PrintBigMat(Phi * MkTot * Omega, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Phi M Omega");
    //PrintBigMat((Phi - Psi) + (P * Omega), SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "(Phi - Psi) + P Omega");

    /////////////////////////////////////////////
    // Get Upsilon and Psi bar using [Jain 1997] recursion
    /////////////////////////////////////////////
    SimTK::Matrix Upsilon(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Mat66 UpsK(0.0);
    SimTK::Mat66 PsistarKKm(0.0);
    SimTK::Mat66 PsiKKm(0.0);
    SimTK::Mat66 DInv0(0.0);
    SimTK::Mat11 DInvK(0.0);
    SimTK::Mat66 Hstar0(0.0);
    SimTK::Mat66 H0(0.0);
    SimTK::Mat61 HstarK(0.0);
    SimTK::Mat16 HK(0.0);

    UpsK = 0;
    for(int k = 1; k < (NBODIESplusG); k++){
        //std::cout << "BODY K = " << k << std::endl;
        // Get PsistarKKm PsiKKm DInvK
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                //std::cout << "PsistarKKm i j " << i << " " << j << " PsistarKKm[][] " << ((k + 1)*6) + i << " " << (k*6) + j << std::flush << std::endl; 
                PsistarKKm[i][j] = Psistar[ (k*6) + i ][ ((k-1)*6) + j ];
                //PsiKKm[i][j] = Psi[ ((k + 1)*6) + i ][ (k*6) + j ];
            }
        }
        PsiKKm = PsistarKKm.transpose();
        //PrintBigMat(PsistarKKm, 6, 6, 2, "PsistarKKm");
        //PrintBigMat(PsiKKm, 6, 6, 2, "PsiKKm");

        // Get HstarK and HK
        if(k == 0){
            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    //std::cout << "Hstar0 i j " << i << " " << j << " = 0 " << std::flush << std::endl; 
                    //std::cout << "DInv0 i j " << i << " " << j << " = 0 " << std::flush << std::endl; 
                    Hstar0[i][j] =  0.0;;
                    DInv0[i][j] =  0.0;;
                }
            }
            H0 = Hstar0.transpose();
        }else if(k == 1){
            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    //std::cout << "Hstar0 i j " << i << " " << j << " Hstar[][] " << (k*6) + i << " " << j << std::flush << std::endl; 
                    //std::cout << "DInv0 i j " << i << " " << j << " DInv[][] " << ((k-1)*6) + i << " " << j << std::flush << std::endl; 
                    Hstar0[i][j] =  Hstar[ (k*6) + i ][ j ];
                    DInv0[i][j] =  DInv[ ((k-1)*6) + i ][ j ];
                }
            }
            H0 = Hstar0.transpose();
        }else if(k > 1){
            //std::cout << "DInvK i j " << 0 << " " << 0 << " DInv " << 5 + (k - 1) << " " << 5 + (k - 1) << std::flush << std::endl;
            DInvK[0][0] =  DInv[ 5 + (k - 1) ][5 + (k - 1)];
            for(int i = 0; i < 6; i++){
                //std::cout << "HstarK i j " << i << " " << 0 << " Hstar " << (6*k) + i << " " << 5 + (k - 1) << std::flush << std::endl;
                HstarK[i][0] =  Hstar[ (6*k) + i ][5 + (k - 1)];
            }
            HK = HstarK.transpose();
        }

        if(k <= 1){
            UpsK = (PsistarKKm * UpsK * PsiKKm) + (Hstar0 * DInv0 * H0);
        }else{
            UpsK = (PsistarKKm * UpsK * PsiKKm) + (HstarK * DInvK * HK);
        }

        //std::cout << "Upsilon(" << k << ") = " << std::endl;
        //std::cout << UpsK << std::endl;
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                Upsilon[(k*6) + i][(k*6) + j] = UpsK[i][j];
            }
        }

    }

    
    // Get PsiBar and Check 4.15
    SimTK::Matrix IdPsi(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    for(int i = 0; i < SPATIAL_DOFSplusG; i++){IdPsi[i][i] = 1.0;}
    SimTK::Matrix PsiBar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix PsiBarstar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);

    PsiBar = Psi - IdPsi;
    PsiBarstar = PsiBar.transpose();

    //PrintBigMat(Omega, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Omega");
    //PrintBigMat(Upsilon + (PsiBarstar * Upsilon) + (Upsilon * PsiBar), SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Omega(Upsilon)");

    // Check 4.17
    //PrintBigMat(P, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "P");
    //PrintBigMat(Upsilon, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Upsilon");
    //PrintBigMat(Hstar, SPATIAL_DOFSplusG, NU, 2, "Hstar"); 

    /////////////////////////////////////////////
    // Get Fixman torque using [Jain 1997] Eq. 4.17 and setting 
    // h(i) to Hstar(k)
    /////////////////////////////////////////////
    SimTK::Mat66 PKK(0.0);
    SimTK::Mat66 UKK(0.0);
    SimTK::Mat66 PYH(0.0);
    /*
    SimTK::Mat66 Hi(0.0);
    SimTK::Vec3 RotVec(0.0);
    SimTK::Vec3 RotVec_G(0.0);
    SimTK::Mat33 crossH(0.0);

    for(int k = 1; k < (NBODIESplusG); k++){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(k));
        std::cout << "BODY K = " << k << std::endl;
        // Get PsistarKKm PKK
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                PKK[i][j] = P[ (k*6) + i ][ (k*6) + j ];
                UKK[i][j] = Upsilon[ (k*6) + i ][ (k*6) + j ];
            }
        }

        // Get HstarK
        if(k == 0){
            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    Hstar0[i][j] =  0.0;
                }
            }
        }else if(k == 1){
            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    Hstar0[i][j] =  Hstar[ (k*6) + i ][ j ];
                }
            }
        }else if(k > 1){
            for(int i = 0; i < 6; i++){
                HstarK[i][0] =  Hstar[ (6*k) + i ][5 + (k - 1)];
            }
        }

        if(k <= 1){
            for(int ExtDof = 0; ExtDof < 6; ExtDof++){
                for (int i = 0; i < 3; i++){
                    RotVec[i] = Hstar0[i][ExtDof];
                }
                const SimTK::Rotation R_GB = mobod.getBodyRotation(someState);
                RotVec_G = R_GB * RotVec;
    
                crossH = SimTK::crossMat(RotVec_G);
                Hi = SimTK::Mat66(0.0);
    
                for (int i = 0; i < 3; i++){
                    for (int j = 0; j < 3; j++){
                        Hi[i][j] = Hi[i +3][j + 3] = crossH[i][j];
                    }
                }
    
                PYH = PKK * UKK * Hi;
                std::cout << std::setprecision(10) << "Tr{PYH} ExtDof " << ExtDof << " = " << PYH.trace() << std::endl << std::setprecision(2);
            }
        }else{
            for (int i = 0; i < 3; i++){
                RotVec[i] = HstarK[i][0];
            }
            const SimTK::Rotation R_GB = mobod.getBodyRotation(someState);
            RotVec_G = R_GB * RotVec;

            //crossH = SimTK::crossMat(RotVec_G);
            crossH = SimTK::crossMat(RotVec);
            Hi = SimTK::Mat66(0.0);

            for (int i = 0; i < 3; i++){
                for (int j = 0; j < 3; j++){
                    Hi[i][j] = Hi[i + 3][j + 3] = crossH[i][j];
                }
            }
            //std::cout << "Hi: " << Hi << std::endl;

            PYH = PKK * UKK * Hi;
            std::cout << std::setprecision(10) << "Tr{PYH} = " << PYH.trace() << std::endl << std::setprecision(2);
        }

    }
    */




    /////////////////////////////////////////////
    // Equations of motion (dynamics) [Jain 1995] Eq. 5
    /////////////////////////////////////////////
    const SimTK::MultibodySystem& multibodySystem = forces->getMultibodySystem();

    // Declare Applied forces = mobility and body
    SimTK::Vector mobilityForces(matter->getNumBodies());
    SimTK::Vector_< SimTK::SpatialVec > rigidBodyForces(matter->getNumBodies());

    // Declare bias forces
    SimTK::Vector_< SimTK::SpatialVec > totalCentrifugalForces(matter->getNumBodies());
    //SimTK::Vector_< SimTK::SpatialVec > totalCentrifugalForces(NU);
    
    // Declare Constraint forces
    SimTK::Vector constraintGeneralizedForces(matter->getNumBodies());
    SimTK::Vector_< SimTK::SpatialVec > constraintSpatialForces(matter->getNumBodies());

    // Declare Unknown forces
    SimTK::Vector_< SimTK::SpatialVec > mobReactF_AtMInG(matter->getNumBodies());
    SimTK::Vector mobReactF_AtMInGAsVector(matter->getNumBodies() * 6);
    SimTK::Matrix mobReactF_AtMInGAsMatrix(SPATIAL_DOFSplusG, NU, 0.0);
    SimTK::Vector motionForces(matter->getNumBodies()); 
    SimTK::Vector_< SimTK::SpatialVec > constraintForces(matter->getNumBodies()); // TODO

    // Get Applied forces = mobility and body
    mobilityForces = multibodySystem.getMobilityForces(someState, SimTK::Stage::Dynamics);
    rigidBodyForces = multibodySystem.getRigidBodyForces(someState, SimTK::Stage::Dynamics);

    // Declare bias forces


    // Get Constraint forces
    matter->findConstraintForces(someState, constraintSpatialForces, constraintGeneralizedForces);	

    // Get Unknown forces
    matter->findMotionForces(someState, motionForces);
    matter->calcMobilizerReactionForces(someState, mobReactF_AtMInG);

    // Print loop
    for(int k = 0; k < (NBODIESplusG); k++){
        std::cout << "body " << k << std::endl;
        SimTK::MobilizedBodyIndex mbx = SimTK::MobilizedBodyIndex(k);
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

        // Print Applied forces = mobility and body
        std::cout <<  std::setprecision(10) << "MobilityForces: " << mobilityForces[k] << std::endl;

        SimTK::Vector rigidBodyForcesAsVector_K(6);
        SOA_SpatialVec2Vector(rigidBodyForces[k], rigidBodyForcesAsVector_K);
        PrintSpatialVec(rigidBodyForces[k], 10, "RigidBodyForce:");
       
        // Print Bias forces
        totalCentrifugalForces[k] = matter->getTotalCentrifugalForces(someState, mbx);
        PrintSpatialVec(totalCentrifugalForces[k], 10, "TotalCentrifugalForces:");

        PrintSpatialVec(matter->getTotalCoriolisAcceleration(someState, mbx), 10, "TotalCoriolisAcceleration:");
        PrintSpatialVec(matter->getGyroscopicForce(someState, mbx), 10, "GyroscopicForce:");

        // Print Constraint forces
        SimTK::Vector constraintSpatialForcesAsVector_K(6);
        SOA_SpatialVec2Vector(constraintSpatialForces[k], constraintSpatialForcesAsVector_K);
        PrintBigMat(constraintSpatialForcesAsVector_K, 6, 10, "constraintSpatialForcesAsVector_K:");
        std::cout << std::setprecision(10) << "constraintGeneralizedForces: " << constraintGeneralizedForces[k] << std::endl;

        // Print Unknown forces
        SimTK::Vector mobReactF_AtMInGAsVector_K(6);
        SOA_SpatialVec2Vector(mobReactF_AtMInG[k], mobReactF_AtMInGAsVector_K);
        PrintBigMat(mobReactF_AtMInGAsVector_K, 6, 10, "MobilizerReactionForce_AtMInGAsVector_K:");
        std::cout << std::setprecision(10) << "MotionForces: " << motionForces[k] << std::endl;

        // Get a block from a convertion matrix (H, H_FM,) 
        //int HBlockNRows = 6, HBlockNCols;
        //(k <= 1) ? (HBlockNCols = 6) : (HBlockNCols = 1);
        //SimTK::Matrix HBlock(HBlockNRows, HBlockNCols);
        //SOA_GetHstarLikeElement(Hstar, k, HBlock);
        std::cout << "===================" << std::endl;
    }
   
    SimTK::Vector generalizedForces(NU);
    //matter->multiplyBySystemJacobianTranspose(someState, mobReactF_AtMInG, generalizedForces);
    matter->multiplyBySystemJacobianTranspose(someState, totalCentrifugalForces, generalizedForces);
    PrintBigMat(generalizedForces, NU, 10, "Resulting generalized forces = SystemJacobianT * TotalCentrifugalForces:");
 
    // Put UDot Vector in a matrix form
    //SimTK::Matrix UDotAsMatrix(NU, 1, 0.0);
    //for(int u = 0; u < (NU); u++){
    //    UDotAsMatrix[u][0] = someState.getUDot()[u];
    //}
    //PrintBigMat( (M * UDotAsMatrix), NU, 1, 5, "M * UDot");

    SimTK::Vector Ma(NU);
    matter->multiplyByM(someState, someState.getUDot(), Ma);
    PrintBigMat( Ma, NU, 1, 10, "M * UDot");





    /////////////////////////////////////////////
    // Eq. (49) 1995 Jain: ThetaDot* Mtheta ThetaDot = 2 H VCross* Phi M V
    /////////////////////////////////////////////
    std::cout << "Equation (49) Jain 1995" << std::endl;

    // First get V and VCross: V = col{V(k)}, VCross = diag{VCross(k)}, VCross(k) = 
    //SimTK::Vector V_FM(6 * NBODIESplusG);
    SimTK::Vector V_GB(6 * NBODIESplusG);
    SimTK::Matrix crossV_GB((6 * NBODIESplusG), (6 * NBODIESplusG), 0.0);
    for(int k = 0; k < (NBODIESplusG); k++){
        SimTK::MobilizedBodyIndex mbx = SimTK::MobilizedBodyIndex(k);
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

        //SimTK::SpatialVec V_FMk = mobod.getMobilizerVelocity(someState); // = getH_FMCol(i in getNumU()) * getOneU()
        SimTK::SpatialVec V_GBk = mobod.getBodyVelocity(someState); // 
        //SimTK::Transform  X_GBk = mobod.getBodyTransform(someState);

        //std::cout << "mobod " << int(mbx) << " V_GB( "<< int(mbx) << " ) = " << V_GBk << std::endl;
        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 3; j++){
                V_GB[(6 * k) + (i * 3) + j] = V_GBk[i][j];
            }
        }

        // Set crossVk according to Eq (34) Jain 1995
        SimTK::SpatialMat crossV_GBk;
        crossV_GBk[0][0] = SimTK::crossMat(V_GBk[0]);
        crossV_GBk[1][0] = SimTK::crossMat(V_GBk[1]);
        crossV_GBk[0][1] = SimTK::Mat33(0.0);
        crossV_GBk[1][1] = SimTK::crossMat(V_GBk[0]);
        //std::cout << "crossV_GB( "<< int(mbx) << " ) = " << crossV_GBk << std::endl;

        // Set Diagonal V_GBCross
        for(int m = 0; m < 2; m++){
            for(int n = 0; n < 2; n++){
                for(int o = 0; o < 3; o++){
                    for(int p = 0; p < 3; p++){
                        crossV_GB[(6 * k) + (m * 3) + o][(6 * k) + (n * 3) + p] = crossV_GBk[m][n][o][p];
                    }
                }
            }
        }

    }
    //PrintBigMat(V_GB - (Phistar * Hstar * someState.getU()), NU, 10, "V_GB - (Phistar * Hstar * U)");
    //PrintBigMat( V_GB, SPATIAL_DOFSplusG, 2, "V_GB");
    //PrintBigMat( crossV_GB, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "crossV_GB");

    SimTK::Vector Eq49RHS(NU); // Eq. 49 right hand side
    Eq49RHS = 2.0 * (( ( (H * (~crossV_GB)) * Phi ) * MkTot ) * V_GB);
    std::cout << std::setprecision(10) << std::endl; 
    std::cout << "Right hand side of (49) = " << Eq49RHS << std::endl;
 
    //



    // Now get MThetaIs for Eq. 4.11 1997 Jain
    // Build HboldI
    //SimTK::Vector TMT(NBODIESplusG); // ThetaDot* x MThetaI x ThetaDot
    SimTK::Vector TMT(NU, 0.0); // ThetaDot* x MThetaI x ThetaDot
    int UIx = -1;
    for (SimTK::MobilizedBodyIndex mbx(0); mbx < NBODIESplusG; ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        int k = int(mbx);
        std::cout << "mobod " << k << std::endl;

                // Get PsistarKKm PKK
                for(int i = 0; i < 6; i++){
                    for(int j = 0; j < 6; j++){
                        PKK[i][j] = P[ (k*6) + i ][ (k*6) + j ];
                        UKK[i][j] = Upsilon[ (k*6) + i ][ (k*6) + j ];
                    }
                } // 
        
        for(int mobodUIx = 0; mobodUIx < mobod.getNumU(someState); mobodUIx++){
            UIx++;
            std::cout << "UIx " << UIx << " mobility " << mobodUIx << std::endl;

            SimTK::Matrix HBoldIDelta(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
            SimTK::Matrix MThetaI(NU, NU, 0.0);
            SimTK::SpatialMat Hi(0);
            SimTK::Mat66 HiAsMat(0);
            SimTK::Vec4 h4;
            SimTK::Vec3 h, h_G;
        
            //if(k == 1){
            if(k >= 0){
                SimTK::SpatialVec HCol = mobod.getHCol(someState, SimTK::MobilizerUIndex(mobodUIx));
                //SimTK::SpatialVec H_FMCol = mobod.getH_FMCol(someState, SimTK::MobilizerUIndex(k));
                //std::cout << "H_FMCol: " << H_FMCol << std::endl;
                h_G[0] = HCol[0][0];
                h_G[1] = HCol[0][1];
                h_G[2] = HCol[0][2];
            
            }else{ // reciprocal to if(k > 1)
                // Get body's inboard mobilizer frame M measured and expressed in the parent body's corresponding outboard frame F
                SimTK::Transform X_FMk = mobod.getMobilizerTransform(someState);
        
                // Get other various transforms
                SimTK::Transform X_GPk;
                SimTK::Transform X_PFk;
                SimTK::Transform X_GFk;
       
                if(k > 0){ 
                    X_GPk = mobod.getParentMobilizedBody().getBodyTransform(someState);
                }
                X_PFk = mobod.getDefaultInboardFrame();
                X_GFk= X_GPk * X_PFk;
                //PrintBigMat(X_FMk.toMat44(), 4, 4, 3, "X_FM");
        
                // Joint rotation axis expressed in parents outboard frame F
                //SimTK::Vec3 rotAngles = X_FMk.R().convertThreeAxesRotationToThreeAngles(SimTK::BodyOrSpaceType::BodyRotationSequence, SimTK::XAxis, SimTK::YAxis, SimTK::ZAxis);
                //std::cout << "rotation angles on XYZ: " << rotAngles << std::endl;
                h4 = X_FMk.R().convertRotationToAngleAxis();
        
                // Get rotation axis. First value of h4 is the angle 
                h[0] = h4[1];
                h[1] = h4[2];
                h[2] = h4[3];
        
                // Reexpress h in Ground frame : it's the same as get HCol
                h_G = X_GFk.R() * h;
            } // end of else reciprocal of if > 1

            if(h_G.norm() > SimTK::TinyReal){ 
                h_G = h_G.normalize();
            }
            std::cout << "h_G for mobod " << k << ": " << h_G << std::endl;
    
            SimTK::Mat33 crossH = SimTK::crossMat(h_G);
            //SimTK::Mat33 crossH = SimTK::crossMat(h);
    
            Hi[0][0] = crossH;
            Hi[1][1] = crossH;
            //PrintSpatialMat(Hi, 2, "Hi");
            //std::cout << "Hi " << Hi << std::endl;
    
            SOA_SpatialMat2Mat66(Hi, HiAsMat);
            //PrintBigMat(HiAsMat, 6, 6, 3, "HiAsMat");
    
            // Set Diagonal HboldI
            for(int m = 0; m < 2; m++){
                for(int n = 0; n < 2; n++){
                    for(int o = 0; o < 3; o++){
                        for(int p = 0; p < 3; p++){
                            HBoldIDelta[(6 * k) + (m * 3) + o][(6 * k) + (n * 3) + p] = Hi[m][n][o][p];
                        }
                    }
                }
            }
            //std::cout << "HBoldIDelta for mobod " << int(mbx) << std::endl;
            //PrintBigMat(HBoldIDelta, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "HBoldIDelta:");
    
            MThetaI = H * Phi * ( (HBoldIDelta * Phi * MkTot) - (MkTot * Phistar * HBoldIDelta) ) * Phistar * Hstar;

            //std::cout << "MThetaI for mobod " << int(mbx) << std::endl;
            //PrintBigMat(MThetaI, NU, NU, 2, "MThetaI:");
    
            // U and QDot are the same in our situation because we only set last u
            TMT[UIx] = (~someState.getU()) * (MThetaI * someState.getU());
            std::cout << std::setprecision(5) << std::endl; 
    
            PYH = PKK * UKK * HiAsMat;
            std::cout << std::setprecision(10) << "Tr{PYH} = " << PYH.trace() << std::endl << std::setprecision(2);

            // Check 
            SimTK::Matrix MInvMThetaI(NU, NU, 0.0);
            SimTK::Real Tr_MInvMThetaI = 0.0;
            MInvMThetaI = MInv * MThetaI;
            for(int i = 0; i < NU; i++){
                Tr_MInvMThetaI += MInvMThetaI[i][i];
            }
            std::cout << std::setprecision(10) << std::endl; 
            std::cout << "0.5 * MInv * MThetaI " <<  0.5 * Tr_MInvMThetaI << std::endl;

         } // mobod mobility
            
    } // mobod
    
    std::cout << std::setprecision(10) << std::endl; 
    std::cout << "TMT: " << TMT << std::endl;

    //std::cout << "Eq. (5) should be equal to Coriolis term" << dM * someState.getU() - 0.5 * TMT << std::endl;

    //PrintBigMat(HboldI, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "HboldI:");
    //PrintBigMat(MThetaI, NU, NU, 2, "MThetaI:");

    /*
    // Check Eq. 411
    SimTK::Matrix MPhistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    MPhistar = MkTot * Phistar;

    SimTK::Matrix HboldPhiM(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix MPhistarHbold(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix SqParan(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix dMdThetaK4_11(NU, NU);

    HboldPhiM = Hbold * Phi * MkTot;
    MPhistarHbold = MkTot * Phistar * Hbold;

    SqParan = HboldPhiM - MPhistarHbold;

    dMdThetaK4_11 = H * Phi * SqParan * Phistar * Hstar;

    //PrintBigMat(HboldPhiM, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "HboldPhiM");
    //PrintBigMat(MPhistarHbold, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "MPhistarHbold");

    //PrintBigMat(SqParan, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "SqParan");
    //PrintBigMat(dMdThetaK, NU, NU, 4, "dMdThetaK");
    //PrintBigMat(dMdThetaK4_11, NU, NU, 4, "dMdThetaK4_11");

    //SimTK::Matrix I(NU, NU, 0.0);
    //for(int i = 0; i < NU; i++){
    //    I[i][i] = 1.0;
    //}
    //SimTK::Matrix MSqParan(NU, NU, 0.0);
    //MSqParan = I - (H * Psi * K);
    //PrintBigMat( MSqParan.transpose() * DInv * MSqParan, NU, NU, 2, "MInvcheck 3.5c");
    */

    /*
    SimTK::Mat66 IK;
    for(int i = 0; i < 6; i++){IK[i][i] = 1.0;}
    SimTK::Mat66 Hstar0;
    SimTK::Mat66 H0;
    SimTK::Mat61 HstarK;
    SimTK::Mat16 HK;

    for(int k = 1; k < NBODIESplusG; k++){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(k));
        SimTK::SpatialVec HCol = mobod.getHCol(someState, SimTK::MobilizerUIndex(k));

        SimTK::Mat66 PsiK, GK;

        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                PsiK[i][j] = Phi[       (k * 6) + i ][ ((k - 1) * 6) + j ];
                GK[i][j]   =   G[ ((k - 1) * 6) + i ][ ((k - 1) * 6) + j ];
            }
        }
 
        if( (k - 1) <= 1){
            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    //std::cout << "Hstar0 i j " << i << " " << j << " Hstar[][] " << ((k - 1) * 6) + i << " " << j << std::endl; 
                    Hstar0[i][j] =  Hstar[ ((k - 1) * 6) + i ][ j ];
                }
            }
        }else{
            for(int i = 0; i < 6; i++){
                //std::cout << "HstarK i j " << i << " " << 0 << " HstarK " << (6 * (k - 1)) + i << " " << 5 + (k - 2) << std::endl;
                HstarK[i][0] =  Hstar[ (6 * (k - 1)) + i ][5 + (k - 2)];
            }
        }

        if( (k - 1) <= 1){
            //std::cout << "Hstar0 " << Hstar0 << std::endl;
            H0 = Hstar0.transpose();       
            PsiK = PsiK * (IK - (GK * H0));
        }else{
            HK = HstarK.transpose();       
            std::cout << "HstarK " << HstarK << std::endl;
            PsiK = PsiK * (IK - (GK * HK));
        }

    }
    */



    // END SOA

    #ifdef HARMONICOSCILLATOR
    /////////////////////////
    // Harmonic oscillator
    /////////////////////////
    //if (( (pe_n) < (pe_o) ) || (rand_no < exp(-( (pe_n) - (pe_o) ) / RT))){ // Accept
    if ((HO_etot_xprop < HO_etot_x) || (rand_no < exp(-(HO_etot_xprop - HO_etot_x)/RT))){ // Accept
        for(int i = 0; i < HO_D; i++){HO_xini[i] = HO_xprop[i];}
        HO_PE_set = HO_PE_xprop;
        HO_KE_set = HO_KE_xprop;
        HO_etot_set = HO_PE_xprop + HO_KE_xprop;
        HO_etot_set = HO_etot_x;
    }else{ // Reject
        for(int i = 0; i < HO_D; i++){HO_x[i] = HO_xini[i];}
    }
    std::cout << "HO 1 pe_os " << HO_PE_set << " ke_os " << HO_KE_set << " etot_os " << HO_etot_set << std::endl;
    #endif

}







