#include <gtest/gtest.h>
#include <vector>

#include "SimTKsimbody.h"

using namespace SimTK;
using namespace std;

static const int NUM_BODIES = 10;
static const Real BOND_LENGTH = 0.5;
static const int ITERATIONS = 4;
static const Real TOL = 1e-4;

auto testFitting(const MultibodySystem& mbs,
                 State& state,
                 const vector<MobilizedBodyIndex>& bodyIxs,
                 const vector<vector<Vec3>>& stations,
                 const vector<vector<Vec3>>& targetLocations,
                 Real minError,
                 Real maxError,
                 Real endDistance) -> bool {
    // Find the best fit.
    const Real reportedError =
        ObservedPointFitter::findBestFit(mbs, state, bodyIxs, stations, targetLocations, TOL);
    const bool result = (reportedError <= maxError && reportedError >= minError);
    EXPECT_TRUE(result) << "reported error=" << reportedError << " not in [" << minError << ", " << maxError
                        << "]";

    // Verify that the error was calculated correctly.
    Real error = 0.0;
    int numStations = 0;
    mbs.realize(state, Stage::Position);
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
    for (int i = 0; i < (int)bodyIxs.size(); ++i) {
        MobilizedBodyIndex id = bodyIxs[i];
        numStations += (int)stations[i].size();
        for (int j = 0; j < (int)stations[i].size(); ++j) {
            error +=
                (targetLocations[i][j] - matter.getMobilizedBody(id).getBodyTransform(state) * stations[i][j])
                    .normSqr();
        }
    }
    error = std::sqrt(error / numStations);
    EXPECT_LT(std::abs(1.0 - (error / reportedError)), 0.0001); // should match to machine precision

    // Verify that the ends are the correct distance apart.
    if (endDistance >= 0) {
        Real distance = (matter.getMobilizedBody(bodyIxs[0]).getBodyOriginLocation(state)
                         - matter.getMobilizedBody(bodyIxs[bodyIxs.size() - 1]).getBodyOriginLocation(state))
                            .norm();
        EXPECT_LT(std::abs(1.0 - (endDistance / distance)), TOL);
    }

    return result;
}


static void testObservedPointFitter(bool useConstraint) {
    int failures = 0;
    for (int iter = 0; iter < ITERATIONS; ++iter) {
        // Build a system consisting of a chain of bodies with occasional side chains, and
        // a variety of mobilizers.
        MultibodySystem mbs;
        SimbodyMatterSubsystem matter(mbs);
        Body::Rigid body = Body::Rigid(MassProperties(1, Vec3(0), Inertia(1)));
        body.addDecoration(Transform(), DecorativeSphere(.1));

        MobilizedBody* lastBody = &matter.Ground();
        MobilizedBody* lastMainChainBody = &matter.Ground();

        vector<MobilizedBody*> bodies;
        Random::Uniform random(0.0, 1.0);
        random.setSeed(iter);

        for (int i = 0; i < NUM_BODIES; ++i) {
            bool mainChain = random.getValue() < 0.5;
            MobilizedBody* parent = (mainChain ? lastMainChainBody : lastBody);
            int type = (int)(random.getValue() * 4);
            MobilizedBody* nextBody;
            if (type == 0) {
                MobilizedBody::Cylinder cylinder(*parent,
                                                 Transform(Vec3(0, 0, 0)),
                                                 body,
                                                 Transform(Vec3(0, BOND_LENGTH, 0)));
                nextBody = &matter.updMobilizedBody(cylinder.getMobilizedBodyIndex());
            } else if (type == 1) {
                MobilizedBody::Slider slider(*parent,
                                             Transform(Vec3(0, 0, 0)),
                                             body,
                                             Transform(Vec3(0, BOND_LENGTH, 0)));
                nextBody = &matter.updMobilizedBody(slider.getMobilizedBodyIndex());
            } else if (type == 2) {
                MobilizedBody::Ball ball(*parent,
                                         Transform(Vec3(0, 0, 0)),
                                         body,
                                         Transform(Vec3(0, BOND_LENGTH, 0)));
                nextBody = &matter.updMobilizedBody(ball.getMobilizedBodyIndex());
            } else {
                MobilizedBody::Pin pin(*parent,
                                       Transform(Vec3(0, 0, 0)),
                                       body,
                                       Transform(Vec3(0, BOND_LENGTH, 0)));
                nextBody = &matter.updMobilizedBody(pin.getMobilizedBodyIndex());
            }
            bodies.push_back(nextBody);
            if (mainChain) {
                lastMainChainBody = nextBody;
            }
            lastBody = nextBody;
        }

        mbs.realizeTopology();
        State state = mbs.getDefaultState();

        matter.setUseEulerAngles(state, true);
        mbs.realizeModel(state);

        // Choose a random initial conformation.
        vector<Real> targetQ(state.getNQ(), Real(0));
        for (MobilizedBodyIndex mbx(1); mbx < matter.getNumBodies(); ++mbx) {
            const MobilizedBody& mobod = matter.getMobilizedBody(mbx);
            for (int i = 0; i < mobod.getNumQ(state); ++i) {
                const QIndex qx0 = mobod.getFirstQIndex(state);
                state.updQ()[qx0 + i] = targetQ[qx0 + i] = 2.0 * random.getValue();
            }
        }

        mbs.realize(state, Stage::Position);

        // Select some random stations on each body.
        vector<vector<Vec3>> stations(NUM_BODIES);
        vector<vector<Vec3>> targetLocations(NUM_BODIES);
        vector<MobilizedBodyIndex> bodyIxs;
        for (int i = 0; i < NUM_BODIES; ++i) {
            MobilizedBodyIndex id = bodies[i]->getMobilizedBodyIndex();
            bodyIxs.push_back(id);
            int numStations = 1 + (int)(random.getValue() * 4);
            for (int j = 0; j < numStations; ++j) {
                Vec3 pos((2.0 * random.getValue()) - 1.0,
                         (2.0 * random.getValue()) - 1.0,
                         (2.0 * random.getValue()) - 1.0);
                stations[i].push_back(pos);
                targetLocations[i].push_back(bodies[i]->getBodyTransform(state) * pos);
            }
        }

        // Add a constraint fixing the distance between the first and last bodies
        Real distance = -1;
        if (useConstraint) {
            Real distance = (bodies[0]->getBodyOriginLocation(state)
                             - bodies[NUM_BODIES - 1]->getBodyOriginLocation(state))
                                .norm();
            Constraint::Rod(*bodies[0], Vec3(0), *bodies[NUM_BODIES - 1], Vec3(0), 1.001 * distance);
        }
        state = mbs.realizeTopology();
        matter.setUseEulerAngles(state, true);
        mbs.realizeModel(state);

        // Try fitting it
        State initState = state;
        const bool initialFit =
            testFitting(mbs, state, bodyIxs, stations, targetLocations, 0.0, 0.03, distance);
        EXPECT_TRUE(initialFit) << "noisy fitting failed at iter " << iter << " q=" << initState.getQ()
                                << "\n";

        // Now add random noise to the target locations, and see if it can still fit decently.
        Random::Gaussian gaussian(0.0, 0.15);
        for (int i = 0; i < (int)targetLocations.size(); ++i) {
            for (int j = 0; j < (int)targetLocations[i].size(); ++j) {
                targetLocations[i][j] += Vec3(gaussian.getValue(), gaussian.getValue(), gaussian.getValue());
            }
        }

        // Start from same config as before
        state = initState;
        const bool noisyFit = testFitting(mbs, state, bodyIxs, stations, targetLocations, 0.1, 0.5, distance);
        EXPECT_TRUE(noisyFit) << "noisy fitting failed at iter " << iter << " q=" << initState.getQ() << "\n";
    }

    EXPECT_EQ(failures, 0);
}

TEST(Simbody_ObservedPointFitter_Unconstrained, FindsReasonableFit) {
    testObservedPointFitter(false);
}

TEST(Simbody_ObservedPointFitter_Constrained, RespectsEndDistanceConstraintAndFits) {
    testObservedPointFitter(true);
}
