/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Chris Dembia                                                      *
 * Contributors: Michael Sherman                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "Simbody.h"
#include "TaskSpace.hpp"

using namespace SimTK;

//==============================================================================
// Jacobian
//==============================================================================
void TaskSpace::Jacobian::updateCache(Matrix& cache) const
{
    m_tspace->getMatterSubsystem().calcStationJacobian(
            getState(),
            m_tspace->getMobilizedBodyIndices(),
            m_tspace->getStations(),
            cache);
}

const TaskSpace::JacobianTranspose& TaskSpace::Jacobian::transpose() const
{
    return m_tspace->getJacobianTranspose(getState());
}

Vector TaskSpace::Jacobian::operator*(const Vector& u) const
{
    Vector_<Vec3> Ju;
    m_tspace->getMatterSubsystem().multiplyByStationJacobian(getState(),
            m_tspace->getMobilizedBodyIndices(), m_tspace->getStations(),
            u, Ju);

    // Convert to a Vector.
    Vector out(3 * Ju.size());
    for (int i = 0; i < Ju.size(); ++i) {
        out[3 * i] = Ju[i][0];
        out[3 * i + 1] = Ju[i][1];
        out[3 * i + 2] = Ju[i][2];
    }

    return out;
}


//==============================================================================
// JacobianTranspose
//==============================================================================
void TaskSpace::JacobianTranspose::updateCache(Matrix& cache) const
{
    cache = transpose().value().transpose();
}

const TaskSpace::Jacobian& TaskSpace::JacobianTranspose::transpose() const
{
    return m_tspace->getJacobian(getState());
}

Vector TaskSpace::JacobianTranspose::operator*(const Vector_<Vec3>& f_GP) const
{
    Vector f;
    m_tspace->getMatterSubsystem().multiplyByStationJacobianTranspose(
            getState(),
            m_tspace->getMobilizedBodyIndices(),
            m_tspace->getStations(),
            f_GP,
            f);
    return f;
}

Vector TaskSpace::JacobianTranspose::operator*(const Vector& f_GP) const
{
    unsigned int nIn = f_GP.size();
    SimTK_APIARGCHECK1_ALWAYS(nIn % 3 == 0,
            "TaskSpace::JacobianTranspose", "operator*",
            "Length of f_GP, which is %i, is not divisible by 3.", nIn);

    unsigned int nOut = nIn / 3;

    // Create the Vector_<Vec3>.
    // TODO debug, or look for methods that already do this.
    Vector_<Vec3> my_f_GP(nOut);
    for (unsigned int i = 0; i < nOut; ++i)
    {
        // getAs is just a recast; doesn't copy.
        my_f_GP[i] = Vec3::getAs(&f_GP[3 * i]);
    }

    // Perform the multiplication.
    return operator*(my_f_GP);
}

Vector TaskSpace::JacobianTranspose::operator*(const Vec3& f_GP) const
{
   return operator*(Vector_<Vec3>(1, f_GP));
}

Matrix TaskSpace::JacobianTranspose::operator*(const Matrix& f_GP) const
{
    unsigned int nrow = getState().getNU();
    unsigned int ncol = f_GP.ncol();

    Matrix out(nrow, ncol);
    for (unsigned int j = 0; j < ncol; ++j)
    {
        // TODO is this cast inefficient? Is it copying?
        out(j) = operator*(Vector(f_GP(j)));
    }

    return out;
}

Matrix TaskSpace::JacobianTranspose::operator*(
        const TaskSpace::Inertia& Lambda) const
{
    // TOOD could be more efficient.
    return operator*(Lambda.value());
}

Matrix TaskSpace::JacobianTranspose::operator*(
        const TaskSpace::DynamicallyConsistentJacobianInverseTranspose& JBarT) const
{
    return operator*(JBarT.value());
}

//==============================================================================
// Inertia
//==============================================================================
void TaskSpace::Inertia::updateCache(Matrix& cache) const
{
    FactorLU inertiaInverse(m_tspace->getInertiaInverse(getState()).value());
    inertiaInverse.inverse(cache);
}

const TaskSpace::InertiaInverse& TaskSpace::Inertia::inverse() const
{
    return m_tspace->getInertiaInverse(getState());
}

Vector TaskSpace::Inertia::operator*(const Vector& a) const
{
    return value() * a;
}

Vector TaskSpace::Inertia::operator*(const Vec3& a) const
{
    return operator*(Vector(a));
}

//==============================================================================
// InertiaInverse
//==============================================================================
void TaskSpace::InertiaInverse::updateCache(Matrix& cache) const
{
    const SimbodyMatterSubsystem& matter = m_tspace->getMatterSubsystem();

    const JacobianTranspose& JT = m_tspace->JT(getState());
    // TODO const Matrix& JT = m_tspace.JT().value();
    const Matrix& J = m_tspace->J(getState()).value();
    /* TODO 
    // TODO cache the result.

    unsigned int nst = m_tspace.getNumScalarTasks();
    unsigned int nu = m_tspace.getState().getNU();

    Matrix J = m_tspace.getJacobian().value();

    Matrix MInvJt(nu, nst);

    for (unsigned int j = 0; j < nst; ++j)
    {
        matter.multiplyByMInv(m_tspace.getState(), J.transpose()(j), MInvJt(j));
    }

    updCacheValue() = J * MInvJt;
    */

    unsigned int nt = m_tspace->getNumTasks();
    unsigned int nst = m_tspace->getNumScalarTasks();
    unsigned int nu = getState().getNU();

    Matrix& inertiaInverse = cache;
    inertiaInverse.resize(nst, nst);

    // Create temporary variables.
    Vector Jtcol(nu);
    Vector MInvJtcol(nu);
    Vector_<Vec3> JMInvJt_j(nt);

    // f_GP is used to pluck out one column at a time of Jt. Exactly one
    // element at a time of f_GP will be 1, the rest are 0.
    Vector f_GP(nst, Real(0));

    for (unsigned int j = 0; j < nst; ++j)
    {
        f_GP[j] = 1;
        Jtcol = JT * f_GP;
        f_GP[j] = 0;

        matter.multiplyByMInv(getState(), Jtcol, MInvJtcol);

        // TODO replace with operator.
        inertiaInverse(j) = J * MInvJtcol;
        /* TODO
        matter.multiplyByStationJacobian(m_tspace.getState(),
                m_tspace.getMobilizedBodyIndices(), m_tspace.getStations(),
                MInvJtcol, JMInvJt_j);

        inertiaInverse(j) = JMInvJt_j;
        */
    }
}

const TaskSpace::Inertia& TaskSpace::InertiaInverse::inverse() const
{
    return m_tspace->getInertia(getState());
}


//==============================================================================
// DynamicallyConsistentJacobianInverse
//==============================================================================
void TaskSpace::DynamicallyConsistentJacobianInverse::updateCache(Matrix& cache)
    const
{
    const JacobianTranspose& JT = m_tspace->getJacobianTranspose(getState());
    const Inertia& Lambda = m_tspace->getInertia(getState());

    // TODO inefficient?
    Matrix JtLambda = JT * Lambda;

    unsigned int nst = m_tspace->getNumScalarTasks();
    unsigned int nu = getState().getNU();

    Matrix& Jbar = cache;
    Jbar.resize(nu, nst);

    for (unsigned int j = 0; j < nst; ++j)
    {
        m_tspace->getMatterSubsystem().multiplyByMInv(getState(),
                JtLambda(j), Jbar(j));
    }
}

const TaskSpace::DynamicallyConsistentJacobianInverseTranspose&
TaskSpace::DynamicallyConsistentJacobianInverse::transpose() const
{
    return m_tspace->getDynamicallyConsistentJacobianInverseTranspose(getState());
}

Vector TaskSpace::DynamicallyConsistentJacobianInverse::operator*(
        const Vector& vec) const
{
    const JacobianTranspose& JT = m_tspace->getJacobianTranspose(getState());
    const Inertia& Lambda = m_tspace->getInertia(getState());

    // TODO where is this even used? TODO test this.

    Vector JBarvec;
    m_tspace->getMatterSubsystem().multiplyByMInv(getState(),
            JT * (Lambda * vec), JBarvec);
    return  JBarvec;
}

Matrix TaskSpace::DynamicallyConsistentJacobianInverse::operator*(
        const Matrix& mat) const
{
    unsigned int nrow = getState().getNU();
    unsigned int ncol = mat.ncol();

    Matrix out(nrow, ncol);
    for (unsigned int j = 0; j < ncol; ++j)
    {
        out(j) = operator*(mat(j).getAsVector());
    }

    return out;
}


//==============================================================================
// DynamicallyConsistentJacobianInverseTranspose
//==============================================================================
void TaskSpace::DynamicallyConsistentJacobianInverseTranspose::updateCache(
            Matrix& cache) const
{

    cache = transpose().value().transpose();
}

const TaskSpace::DynamicallyConsistentJacobianInverse&
TaskSpace::DynamicallyConsistentJacobianInverseTranspose::transpose() const
{
    return m_tspace->getDynamicallyConsistentJacobianInverse(getState());
}

Vector TaskSpace::DynamicallyConsistentJacobianInverseTranspose::operator*(
        const Vector& g) const
{
    // TODO inefficient. can we have an MInvT operator??
    return value() * g;
}


//==============================================================================
// InertialForces
//==============================================================================
void TaskSpace::InertialForces::updateCache(Vector& cache) const
{
    Vector jointSpaceInertialForces;
    m_tspace->getMatterSubsystem().calcResidualForceIgnoringConstraints(
            getState(), Vector(0), Vector_<SpatialVec>(0), Vector(0),
            jointSpaceInertialForces);

    Vector JDotu;
    m_tspace->getMatterSubsystem().calcBiasForStationJacobian(
            getState(),
            m_tspace->getMobilizedBodyIndices(), m_tspace->getStations(),
            JDotu);

    const DynamicallyConsistentJacobianInverseTranspose& JBarT =
        m_tspace->JBarT(getState());
    const Vector& b = jointSpaceInertialForces;
    const Inertia& Lambda = m_tspace->Lambda(getState());

    cache = JBarT * b - Lambda * JDotu;
}

Vector TaskSpace::InertialForces::operator+(const Vector& f) const
{
    return value() + f;
}


//==============================================================================
// NullspaceProjection
//==============================================================================
void TaskSpace::NullspaceProjection::updateCache(Matrix& cache) const
{
    cache = transpose().value().transpose();
}

const TaskSpace::NullspaceProjectionTranspose&
TaskSpace::NullspaceProjection::transpose() const
{
    return m_tspace->getNullspaceProjectionTranspose(getState());
}

Vector TaskSpace::NullspaceProjection::operator*(const Vector& vec)
    const
{
    return vec - (m_tspace->JBar(getState()) * (m_tspace->J(getState()) * vec));
}


//==============================================================================
// NullspaceProjectionTranspose
//==============================================================================
void TaskSpace::NullspaceProjectionTranspose::updateCache(Matrix& cache) const
{
    cache = 1 - (m_tspace->JT(getState()) * m_tspace->JBarT(getState()));
}

const TaskSpace::NullspaceProjection&
TaskSpace::NullspaceProjectionTranspose::transpose() const
{
    return m_tspace->getNullspaceProjection(getState());
}

Vector TaskSpace::NullspaceProjectionTranspose::operator*(const Vector& vec)
    const
{
    return vec - (m_tspace->JT(getState()) * (m_tspace->JBarT(getState()) * vec));
}


//==============================================================================
// TaskSpace
//==============================================================================
// TODO account for applied forces? velocities?
Vector TaskSpace::calcInverseDynamics(const State& s, const Vector& taskAccelerations) const
{
    return  Lambda(s) * taskAccelerations + mu(s);
}


















//==============================================================================
//                             TASKS MEASURE
//==============================================================================
// This Measure is added to the modelRobot and is used to manage the tasks
// it is supposed to achieve and to return as its value the control torques
// that should be applied to the realRobot (that is, the simulated one).
// This should only be instantiated with T=Vector.
template <class T>
class TasksMeasure : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(TasksMeasure, Measure_<T>);

    TasksMeasure(SimTK::GeneralForceSubsystem& force, SimTK::SimbodyMatterSubsystem& matter, MobilizedBodyIndex body, Vec3 station) 
    :   Measure_<T>(force, 
                    new Implementation(matter, body, station),
                    AbstractMeasure::SetHandle()) {}

    const Vec3& getTarget() const { return getImpl().m_desiredTaskPosInGround; }
    Vec3& updTarget() { return updImpl().m_desiredTaskPosInGround; }
    void setTarget(Vec3 pos) { updImpl().m_desiredTaskPosInGround = pos; }

    void toggleGravityComp() {
        updImpl().m_compensateForGravity = !isGravityCompensationOn();}
    void togglePoseControl() {
        updImpl().m_controlPose = !isPoseControlOn();}
    void toggleTask() {updImpl().m_controlTask = !getImpl().m_controlTask;}
    void toggleEndEffectorSensing() 
    {   updImpl().m_endEffectorSensing = !getImpl().m_endEffectorSensing;}

    bool isGravityCompensationOn() const 
    {   return getImpl().m_compensateForGravity; }
    bool isPoseControlOn() const 
    {   return getImpl().m_controlPose; }
    bool isEndEffectorSensingOn() const 
    {   return getImpl().m_endEffectorSensing; }
    bool isTaskPointFollowingOn() const
    {   return getImpl().m_controlTask; }
    const Vec3& getTaskPointInEndEffector() const 
    {   return getImpl().m_taskPointInEndEffector; }

    SimTK_MEASURE_HANDLE_POSTSCRIPT(TasksMeasure, Measure_<T>);
};


template <class T>
class TasksMeasure<T>::Implementation : public Measure_<T>::Implementation {
public:
    Implementation(SimTK::SimbodyMatterSubsystem& matter, MobilizedBodyIndex body, Vec3 station)
    :   Measure_<T>::Implementation(T(), 1),
        m_tspace1(matter),
        m_matter(matter),
        m_body(body),
        m_endEffectorStation(station),
        m_proportionalGain(225),
        m_derivativeGain(30),
        m_jointPositionGain(100),
        m_jointDampingGain(20),
        m_controlPose(true),
        m_controlTask(false),
        m_endEffectorSensing(false),
        m_desiredTaskPosInGround(Vec3(0.4, -0.1, 1)) // Z is up
    {       
        //TODO: should have end effector body
        m_tspace1.addStationTask(body, station);
    }

    // Implementations of virtual methods.
    Implementation* cloneVirtual() const override
    {   return new Implementation(*this); }
    int getNumTimeDerivativesVirtual() const override {return 0;}
    Stage getDependsOnStageVirtual(int order) const override
    {   return Stage::Velocity; }

    // This is the task space controller. It returns joint torques tau as the
    // value of the enclosing Measure.
    void calcCachedValueVirtual(const State& s, int derivOrder, T& tau) const override
    {
        SimTK_ASSERT1_ALWAYS(derivOrder==0,
            "TasksMeasure::Implementation::calcCachedValueVirtual():"
            " derivOrder %d seen but only 0 allowed.", derivOrder);

        // Shorthands.
        // -----------
        const TaskSpace& p1 = m_tspace1;

        const int mnq = s.getNQ();
        const int mnu = s.getNU();
        tau.resize(mnu);

        const Real& kd = m_derivativeGain;
        const Real& kp = m_proportionalGain;

        // The desired task position is in Ground. We need instead to measure it
        // from the real robot's pelvis origin so that we can translate it into the 
        // model's pelvis-centric viewpoint.
        const Transform& X_GP = Transform(); // m_modelRobot.getSampledPelvisPose(s); // this is origin
        const Vec3 x1_des = ~X_GP*m_desiredTaskPosInGround; // measure in P

        // Compute control law in task space (F*).
        // ---------------------------------------
        Vec3 xd_des(0);
        Vec3 xdd_des(0);

        const auto EE = m_matter.getMobilizedBody(m_body).findStationLocationInGround(s, m_endEffectorStation);

        // Get the model's estimate of the end effector location in Ground, which
        // is also the pelvis origin.
        Vec3 x1, x1d;
        p1.findStationLocationAndVelocityInGround(s,
                TaskSpace::StationTaskIndex(0),
                m_modelRobot.getEndEffectorStation(), x1, x1d);

        // if (m_endEffectorSensing) {
        //     // Since the controller model has the pelvis origin fixed at (0,0,0),
        //     // we need to know the real robot's pelvis location so we can measure
        //     // the real robot's end effector from its pelvis location. We don't
        //     // have to modify x1d because we want the end effector stationary
        //     // in Ground, not in the pelvis.
        //     const Vec3& x1_G = m_modelRobot.getSampledEndEffectorPos(s);
        //     x1 = ~X_GP*x1_G; // measure end effector in pelvis frame
        // }

        // Units of acceleration.
        Vec3 Fstar1 = xdd_des + kd * (xd_des - x1d) + kp * (x1_des - x1);

        // Compute task-space force that achieves the task-space control.
        // F = Lambda Fstar + p
        Vector F1 = p1.Lambda(s) * Fstar1 + p1.mu(s);

        // Combine the reaching task with the gravity compensation and pose 
        // control to a neutral q=0 pose with u=0 also.
        const Vector& q = s.getQ();
        const Vector& u = s.getU();
        // const Real k = m_jointPositionGain;
        // const Real c = m_jointDampingGain;
        Vector Mu(mnu), Mq(mnu);
        m_matter.multiplyByM(s, u, Mu);
        m_matter.multiplyByM(s, q, Mq);

        // tau.setToZero();
        // const Real pFac = m_controlPose?1.:0.;
        // if (m_controlTask) {
        //     tau += p1.JT(s) * F1;
        //     tau += p1.NT(s) * (gFac*p1.g(s) - pFac*k*Mq - c*Mu); // damping always
        // } else {
        //     tau += gFac*p1.g(s) - (pFac*k*Mq + c*Mu);
        // }

        // // Cut tau back to within effort limits.
        // // TODO: can't use these limits with one-foot support!
        // const Vector& effortLimits = m_modelRobot.getEffortLimits();
        // for (int i=0; i < mnu; ++i) {
        //     const Real oldtau = tau[i], effort = 10*effortLimits[i]; // cheating
        //     if (std::abs(oldtau) <= effort) continue;
        //     const Real newtau = clamp(-effort, oldtau, effort);
        //     tau[i] = newtau;
        // }

    }

    // TaskSpace objects require some State resources; this call is the time
    // for doing that so forward on to the TaskSpace.
    void realizeMeasureTopologyVirtual(State& modelState) const override {
        m_tspace1.realizeTopology(modelState);
    }
private:
friend class TasksMeasure<T>;

    TaskSpace m_tspace1;
    const SimbodyMatterSubsystem& m_matter;
    MobilizedBodyIndex m_body;
    Vec3 m_endEffectorStation;

    const Real m_proportionalGain;     // task space
    const Real m_derivativeGain;
    const Real m_jointPositionGain;    // joint space
    const Real m_jointDampingGain;

    bool m_controlPose;
    bool m_controlTask;
    bool m_endEffectorSensing;
    Vec3 m_desiredTaskPosInGround;
};

//==============================================================================
//                    REACHING AND GRAVITY COMPENSATION
//==============================================================================
// This is a task-space controller that tries to move the end effector to 
// a particular target location, and applies gravity compensation and some
// joint damping as lower-priority tasks.
//
// The controller has its own internal Atlas model which in general does not
// match the "real" Atlas perfectly. Each time it is asked to
// generate control torques it reads the sensors on the real Atlas, updates
// its internal model to match. It then generates torques that would be right
// for the internal model, but returns them to be applied to the real Atlas.
class ReachingAndGravityCompensation : public Force::Custom::Implementation {
public:
    ReachingAndGravityCompensation(SimTK::GeneralForceSubsystem& force, SimTK::SimbodyMatterSubsystem& matter, MobilizedBodyIndex body, Vec3 station) : m_modelTasks(force, matter, body, station)
    {
    }

    // // Call this after the real robot has been initialized to set up
    // // a mapping between the joint coordinates in the model robot and
    // // the corresponding ones in the real robot.
    // void mapModelToRealRobot(const State& realState) {
    //     const URDFJoints& modelJoints = m_modelRobot.getURDFRobot().joints;
    //     const URDFJoints& realJoints  = m_realRobot.getURDFRobot().joints;

    //     m_model2realU.resize(m_modelState.getNU()); 
    //     m_model2realQ.resize(m_modelState.getNQ());

    //     for (int mj=0; mj < (int)modelJoints.size(); ++mj) {
    //         const URDFJointInfo& modelInfo = modelJoints.getJoint(mj);
    //         const URDFJointInfo& realInfo = realJoints.getJoint(modelInfo.name);
    //         const MobilizedBody& modelMobod = modelInfo.mobod;
    //         const MobilizedBody& realMobod = realInfo.mobod;
    //         const int mnu = modelMobod.getNumU(m_modelState), 
    //                 mnq = modelMobod.getNumQ(m_modelState),
    //                 mu0 = modelMobod.getFirstUIndex(m_modelState),
    //                 mq0 = modelMobod.getFirstQIndex(m_modelState);
    //         if (mnu==0)
    //             continue; // this is fixed in the model; might not be in real robot

    //         const int rnu = realMobod.getNumU(realState), 
    //                 rnq = realMobod.getNumQ(realState),
    //                 ru0 = realMobod.getFirstUIndex(realState),
    //                 rq0 = realMobod.getFirstQIndex(realState);
    //         SimTK_ASSERT1_ALWAYS(mnu==rnu && mnq==rnq,
    //             "ReachingAndGravityCompensation::mapModelToRealRobot(): "
    //             "joint '%s' dof mismatch.", modelInfo.name.c_str());
    //         for (int mu=0; mu < mnu; ++mu)
    //             m_model2realU[mu0+mu] = ru0+mu;
    //         for (int mq=0; mq < mnq; ++mq) 
    //             m_model2realQ[mq0+mq] = rq0+mq;
    //     }

    //     std::cout<<"m2rU="<<m_model2realU<<std::endl;
    //     std::cout<<"m2rQ="<<m_model2realQ<<std::endl;
    // }

    // const Vec3& getTarget() const {return m_modelTasks.getTarget();}
    // Vec3& updTarget() {return m_modelTasks.updTarget();}
    // void setTarget(Vec3 pos) {m_modelTasks.setTarget(pos);}

    // void toggleGravityComp() {m_modelTasks.toggleGravityComp();}
    // void togglePoseControl() {m_modelTasks.togglePoseControl();}
    // void toggleTask() {m_modelTasks.toggleTask();}
    // void toggleEndEffectorSensing() {m_modelTasks.toggleEndEffectorSensing();}

    // bool isGravityCompensationOn() const 
    // {   return m_modelTasks.isGravityCompensationOn(); }

    // This method calculates the needed control torques and adds them into
    // the given mobilityForces Vector which will be applied to the real Atlas.
    // The supplied State is from the real Atlas and will be used to read its
    // sensors.
    void calcForce(const SimTK::State&                realState,
                   SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                   SimTK::Vector_<SimTK::Vec3>&       particleForces,
                   SimTK::Vector&                     mobilityForces) const override
    {
        // // Sense the real robot and use readings to update model robot.
        // // ------------------------------------------------------------
        const int mnq = m_modelState.getNQ(), mnu = m_modelState.getNU();
        // const Vector& sensedQ = m_realRobot.getSampledAngles(realState); // real + noise
        // const Vector& sensedU = m_realRobot.getSampledRates(realState); // real + noise
        // for (int i=0; i < mnq; ++i)
        //     realState.updQ()[QIndex(i)] = sensedQ[m_model2realQ[i]];
        // for (int i=0; i < mnu; ++i)
        //     realState.updU()[UIndex(i)] = sensedU[m_model2realU[i]];

        // // We have to know the pose of the real robot's pelvis so we can figure
        // // out the pelvis-relative location of the end effector, and the effective
        // // gravity direction since the model robot has its pelvis frame welded to
        // // its Ground frame.

        // const Transform& X_GP = m_realRobot.getSampledPelvisPose(realState);
        // m_modelRobot.setSampledPelvisPose(m_modelState, X_GP);

        // // Optional: if real robot end effector location can be sensed, it can
        // // be used to improve accuracy. Otherwise, estimate the end effector
        // // location using the model robot.
        // const Vec3& sensedEEPos = m_realRobot.getSampledEndEffectorPos(realState);
        // m_modelRobot.setSampledEndEffectorPos(m_modelState, sensedEEPos);

        // // Calculate model kinematics.
        // m_modelRobot.realize(m_modelState, Stage::Velocity);
        // // Obtain joint torques from task controller.
        const Vector& tau = m_modelTasks.getValue(realState);

        // Apply model joint torques to their corresponding real robot dofs.
        for (int i=0; i < mnu; ++i)
            mobilityForces[m_model2realU[i]] += tau[i];
    }

    // This controller does not contribute potential energy to the system.
    Real calcPotentialEnergy(const SimTK::State& state) const override
    { return 0; }

private:
    TasksMeasure<Vector> m_modelTasks;
    mutable State        m_modelState;   // Temporary: State for the model robot.

    // Map from model robot coordinates to real robot coordinates.
    Array_<int>          m_model2realQ;
    Array_<int>          m_model2realU;
};