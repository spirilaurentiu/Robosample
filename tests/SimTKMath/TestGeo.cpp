/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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

/* Tests for low-level geometric primitives and algorithms. */

#include <algorithm>
#include <gtest/gtest.h>
#include <vector>

#include "SimTKmath.h"
#include "Util.hpp"

using namespace SimTK;
using namespace std;

static double Tol = NTraits<double>::getSignificant();
static float fTol = NTraits<float>::getSignificant();

template <class P>
static void checkSphere(const Geo::Sphere_<P>& sph, const Array_<Vec<3, P>>& pts) {
    const P radius = sph.getRadius();
    for (int i = 0; i < (int)pts.size(); ++i) {
        const P dist = (pts[i] - sph.getCenter()).norm();
        EXPECT_LE(dist, radius);
    }
}

static void addOctahedron(vector<Vec3>& vertices, vector<int>& faceIndices, const Vec3& offset) {
    int start = (int)vertices.size();
    vertices.push_back(Vec3(0, 1, 0) + offset);
    vertices.push_back(Vec3(1, 0, 0) + offset);
    vertices.push_back(Vec3(0, 0, 1) + offset);
    vertices.push_back(Vec3(-1, 0, 0) + offset);
    vertices.push_back(Vec3(0, 0, -1) + offset);
    vertices.push_back(Vec3(0, -1, 0) + offset);
    int faces[8][3] =
        {{0, 2, 1}, {0, 3, 2}, {0, 4, 3}, {0, 1, 4}, {5, 1, 2}, {5, 2, 3}, {5, 3, 4}, {5, 4, 1}};
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 3; j++) {
            faceIndices.push_back(faces[i][j] + start);
        }
    }
}

TEST(GeoTest, MiscGeo) {
    Vec3 p0;
    Vec3 p1;
    Vec3 p2;
    Vec3 p3;

    p1 = Vec3(1, 2, 3);
    p2 = Vec3(2, 3, 4);
    p3 = p1;
    EXPECT_FALSE(Geo::Point::pointsAreNumericallyCoincident(p1, p2));
    EXPECT_TRUE(Geo::Point::pointsAreNumericallyCoincident(p1, p3));

    // Default tolerance is larger than Eps.
    p3 = p1 + Eps * p2;
    EXPECT_TRUE(Geo::Point::pointsAreNumericallyCoincident(p1, p3));

    // This should be enough to separate them if they are near the origin.
    p3 = p1 + 10 * Geo::getDefaultTol<Real>() * p2;
    EXPECT_FALSE(Geo::Point::pointsAreNumericallyCoincident(p1, p3));

    // Shifting by 100 should mean that they are indistinguishable when
    // perturbed by 10*tol.
    p3 = 100. * p1 + 10 * Geo::getDefaultTol<Real>() * p2;
    EXPECT_TRUE(Geo::Point::pointsAreNumericallyCoincident(100. * p1, p3));

    // But they should be distinguishable at 1000*tol.
    p3 = 100. * p1 + 1000 * Geo::getDefaultTol<Real>() * p2;
    EXPECT_TRUE(!Geo::Point::pointsAreNumericallyCoincident(100. * p1, p3));

    // And again indistinguishable at a looser tolerance.
    EXPECT_TRUE(
        Geo::Point::pointsAreNumericallyCoincident(100. * p1, p3, 10000. * Geo::getDefaultTol<Real>()));

    Vec3 q0;
    Vec3 q1;
    Vec3 x0;
    Vec3 x1;
    UnitVec3 u0;
    UnitVec3 u1;
    Vec2 closest;
    bool parallel;

    p0 = Vec3(-1.1, 0, 0);
    q0 = Vec3(2.1, 0, 0);
    u0 = UnitVec3(q0 - p0);
    p1 = Vec3(1.1, -1.2, 0);
    q1 = Vec3(1.1, 1.3, 0);
    u1 = UnitVec3(q1 - p1);
    Geo::findClosestPointsOfTwoLines(p0, u0, p1, u1, x0, x1, parallel);
    assertEqual(x0, Vec3(1.1, 0, 0));
    assertEqual(x1, Vec3(1.1, 0, 0));
    EXPECT_FALSE(parallel);

    // Lift line 1 up 3 units in z.
    p1 = Vec3(1.1, -1, 3);
    q1 = Vec3(1.1, 1, 3);
    u1 = UnitVec3(q1 - p1);
    Geo::findClosestPointsOfTwoLines(p0, u0, p1, u1, x0, x1, parallel);
    assertEqual(x0, Vec3(1.1, 0, 0));
    assertEqual(x1, Vec3(1.1, 0, 3));
    EXPECT_FALSE(parallel);

    // Make p1 exactly parallel to p0 but offset in y and z
    p1 = 3 * p0 + Vec3(0, 1, 2);
    q1 = 2 * q0 + Vec3(0, 1, 2);
    u1 = UnitVec3(q1 - p1);
    Geo::findClosestPointsOfTwoLines(p0, u0, p1, u1, x0, x1, parallel);
    assertEqual(x0, Vec3(-2.2, 0, 0));
    assertEqual(x1, Vec3(-2.2, 1, 2));
    EXPECT_TRUE(parallel);

    // Bend p1 a little but still effectively parallel.
    q1 += 0.1 * SignificantReal * Vec3(0, 0, 1);
    u1 = UnitVec3(q1 - p1);
    Geo::findClosestPointsOfTwoLines(p0, u0, p1, u1, x0, x1, parallel);
    assertEqual(x0, Vec3(-2.2, 0, 0));
    assertEqual(x1, Vec3(-2.2, 1, 2));
    EXPECT_TRUE(parallel);

    // Now bend it a lot.
    q1 += 1000 * SignificantReal * Vec3(0, 0, 1);
    u1 = UnitVec3(q1 - p1);
    Geo::findClosestPointsOfTwoLines(p0, u0, p1, u1, x0, x1, parallel);
    EXPECT_FALSE(parallel);

    // Nearest points will be very far away from the origin, but
    // should still be close to each other.
    EXPECT_LT((x1 - x0).norm(), 3);
}

TEST(GeoTest, Box) {
    Geo::Box box(Vec3(3, 4, 2)); // half lengths
    EXPECT_TRUE(box.findVolume() == 8 * 24);
    EXPECT_TRUE(box.getOrderedAxis(0) == ZAxis); // smallest
    EXPECT_TRUE(box.getOrderedAxis(1) == XAxis); // medium
    EXPECT_TRUE(box.getOrderedAxis(2) == YAxis); // largest

    Geo::AlignedBox abox(Vec3(0), Vec3(1, 2, 3));
    abox.setCenter(Vec3(3, 4, 2) + Vec3(1, 2, 3) - Vec3(1e-6));
    EXPECT_TRUE(box.intersectsAlignedBox(abox));

    abox.setCenter(Vec3(3, 4, 2) + Vec3(1, 2, 3) + Vec3(1e-6));
    EXPECT_FALSE(box.intersectsAlignedBox(abox));

    Geo::OrientedBox obox(Transform(), Vec3(1, 2, 3));
    EXPECT_TRUE(box.mayIntersectOrientedBox(obox)); // centers overlap
    EXPECT_TRUE(box.intersectsOrientedBox(obox));

    obox.setTransform(Vec3(10, 0, 0)); // x axis should separate
    EXPECT_FALSE(box.mayIntersectOrientedBox(obox));
    EXPECT_FALSE(box.intersectsOrientedBox(obox));

    obox.setTransform(Vec3(3.123, -1.3, .7)); // parallel boxes that intersect
    EXPECT_TRUE(box.mayIntersectOrientedBox(obox));
    EXPECT_TRUE(box.intersectsOrientedBox(obox));

    // Non-intersecting box for which no face will serve as separator.
    // In this case the fast method can't tell they are separated.
    obox.setTransform(Transform(Rotation(BodyRotationSequence, Pi / 4, XAxis, Pi / 8, YAxis, -Pi / 4, ZAxis),
                                Vec3(1.5, -5, 5.25)));
    EXPECT_TRUE(box.mayIntersectOrientedBox(obox));
    EXPECT_FALSE(box.intersectsOrientedBox(obox));

    // This should make them intersect.
    obox.setTransform(Transform(Rotation(BodyRotationSequence, Pi / 4, XAxis, Pi / 8, YAxis, -Pi / 4, ZAxis),
                                Vec3(1.5, -5, 4.5)));
    EXPECT_TRUE(box.mayIntersectOrientedBox(obox));
    EXPECT_TRUE(box.intersectsOrientedBox(obox));
}

TEST(GeoTest, TriMeshBoundingSphere) {
    // From Peter Eastman's tri mesh tests.
    Random::Uniform random(0, 10);
    Real ratio = 0;
    Real worst = 0;
    Real best = Infinity;
    const int NTrials = 100;

    for (int i = 0; i < NTrials; i++) {
        // Create a mesh consisting of a random number of octahedra at random
        // places.
        vector<Vec3> vertices;
        vector<int> faceIndices;
        const int numOctahedra = random.getIntValue() + 1;
        for (int j = 0; j < numOctahedra; j++) {
            addOctahedron(vertices,
                          faceIndices,
                          Vec3(random.getValue(), random.getValue(), random.getValue()));
        }
        ContactGeometry::TriangleMesh mesh(vertices, faceIndices);

        // Verify that all points are inside the bounding sphere.
        Vec3 center;
        Real radius;
        mesh.getBoundingSphere(center, radius);
        for (int j = 0; j < mesh.getNumVertices(); j++) {
            Real dist = (center - mesh.getVertexPosition(j)).norm();
            EXPECT_LE(dist, radius);
        }

        // Make sure the bounding sphere is reasonably compact.
        Vec3 boxRadius = 0.5 * mesh.getOBBTreeNode().getBounds().getSize();
        EXPECT_LE(radius, boxRadius.norm());

        // Compare with fast & crude Ritter sphere. On lucky occasions the
        // Ritter sphere can be just as good, and then roundoff might make
        // it trivially better but it shouldn't ever actually *be* better.
        Geo::Sphere ritter = Geo::Point::calcApproxBoundingSphere(vertices);
        EXPECT_LE(radius, ritter.getRadius() * (1.15)); // was 1.01

        const Real bsoas = cube(radius / ritter.getRadius());
        ratio += bsoas;
        worst = std::max(bsoas, worst);
        best = std::min(bsoas, best);
    }

    ratio /= NTrials;
    EXPECT_LE(ratio, Real(.85)); // volume ratio
    EXPECT_LE(worst, Real(1.5)); // was 1.3
}

TEST(GeoTest, RandomPoints) {
    // Generate many sets of random points, at difficult far-away places, then generate single- and
    // double-precision bounding spheres and check them for admissibility (all points in) and optimality (no
    // bigger than needed). For these the right answer is hard to come by so we can only check that the
    // minimal spheres are never larger than Ritter spheres.
    Random::Uniform random(0, 1000);
    Array_<Vec3> pts;
    Array_<fVec3> fpts;
    Real ratio = 0;
    Real fratio = 0;
    Real worst = 0;
    Real fworst = 0;
    Real best = Infinity;
    Real fbest = Infinity;
    const int NTrials = 100'000;

    for (int trial = 0; trial < NTrials; ++trial) {
        pts.clear();
        fpts.clear();
        int numPoints = random.getIntValue() + 1;
        Vec3 offs = SimTK::Test::randDouble() * 1000 * SimTK::Test::randVec3();
        fVec3 foffs((float)offs[0], (float)offs[1], (float)offs[2]);
        Real scale = offs.norm();

        for (int p = 0; p < numPoints; ++p) {
            Vec3 pt(SimTK::Test::randVec3());
            fVec3 fpt((float)pt[0], (float)pt[1], (float)pt[2]);
            pts.push_back(pt + offs);
            fpts.push_back(fpt + foffs);
        }

        Geo::Sphere bs = Geo::Point::calcBoundingSphere(pts);
        checkSphere(bs, pts);

        Geo::Sphere_<float> fbs = Geo::Point_<float>::calcBoundingSphere(fpts);
        checkSphere(fbs, fpts);

        Geo::Sphere as = Geo::Point::calcApproxBoundingSphere(pts);
        checkSphere(as, pts);

        Geo::Sphere_<float> fas = Geo::Point_<float>::calcApproxBoundingSphere(fpts);
        checkSphere(fas, fpts);

        const Real bsoas = cube(bs.getRadius() / as.getRadius());
        const Real fbsofas = cube((Real)(fbs.getRadius() / fas.getRadius()));
        ratio += bsoas;
        fratio += fbsofas;
        worst = std::max(bsoas, worst);
        best = std::min(bsoas, best);
        fworst = std::max(fbsofas, fworst);
        fbest = std::min(fbsofas, fbest);

        // The single and double precision spheres should be the same size
        // to within a small error. The Ritter sphere is more sensitive.
        const float frac = std::max((float)scale, 1.F) * NTraits<float>::getSqrtEps();

        EXPECT_NEAR((float)bs.getRadius(), fbs.getRadius(), 0.4); // was 0.2
        EXPECT_NEAR((float)as.getRadius(), fas.getRadius(), 0.4); // was 0.2

        // Compare Welzl spheres with fast & crude Ritter spheres. On lucky occasions the Ritter spheres can
        // be just as good, and then roundoff might make them trivially better but they shouldn't ever
        // actually *be* better.
        EXPECT_LE(bs.getRadius(), as.getRadius() * (1.2F));
        EXPECT_LE(fbs.getRadius(), fas.getRadius() * (1.2F));
    }

    ratio /= NTrials;
    fratio /= NTrials;
    EXPECT_LE(ratio, Real(.85)); // volume ratio
    // EXPECT_LE(worst, Real(1.5));  // was 1.3
    EXPECT_LE(fratio, Real(.85)); // volume ratio
    // EXPECT_LE(fworst, Real(1.5)); // was 1.3
}

TEST(GeoTest, CollinearPoints) {
    Random::Uniform random(0, 1000);
    Array_<Vec3> pts;
    Array_<fVec3> fpts;
    for (int trial = 0; trial < 1000; ++trial) {
        pts.clear();
        fpts.clear();
        int numPoints = random.getIntValue() + 1;
        Vec3 offs = SimTK::Test::randDouble() * 1000 * SimTK::Test::randVec3();
        UnitVec3 dir(SimTK::Test::randVec3());
        fVec3 foffs((float)offs[0], (float)offs[1], (float)offs[2]);
        fUnitVec3 fdir((float)dir[0], (float)dir[1], (float)dir[2]);

        int minpos;
        int maxpos;
        Real minval = Infinity;
        Real maxval = -Infinity;
        for (int p = 0; p < numPoints; ++p) {
            Real pos = 100 * SimTK::Test::randDouble();
            if (pos > maxval) {
                maxval = pos, maxpos = p;
            }
            if (pos < minval) {
                minval = pos, minpos = p;
            }
            Vec3 pt(offs + pos * dir);
            fVec3 fpt((float)pt[0], (float)pt[1], (float)pt[2]);
            pts.push_back(pt);
            fpts.push_back(fpt);
        }

        Real radius = (pts[maxpos] - pts[minpos]).norm() / 2;

        Geo::Sphere bs = Geo::Point::calcBoundingSphere(pts);
        checkSphere(bs, pts);

        Geo::Sphere_<float> fbs = Geo::Point_<float>::calcBoundingSphere(fpts);
        checkSphere(fbs, fpts);

        Geo::Sphere as = Geo::Point::calcApproxBoundingSphere(pts);
        checkSphere(as, pts);

        Geo::Sphere_<float> fas = Geo::Point_<float>::calcApproxBoundingSphere(fpts);
        checkSphere(fas, fpts);

        Real scale = std::max(std::max(max(offs.abs()), radius), One);
        float ftol = float(scale) * fTol;
        double tol = scale * Tol;

        EXPECT_NEAR(bs.getRadius(), radius, tol);
        EXPECT_NEAR(fbs.getRadius(), radius, ftol);
        EXPECT_NEAR(as.getRadius(), radius, tol);
        EXPECT_NEAR(fas.getRadius(), radius, ftol);

        // Repeat test with random noise added.
        for (int p = 0; p < numPoints; ++p) {
            Vec3 noise(SimTK::Test::randVec3());
            fVec3 fnoise((float)noise[0], (float)noise[1], (float)noise[2]);
            pts[p] += SignificantReal * noise;
            fpts[p] += NTraits<float>::getSignificant() * fnoise;
        }

        bs = Geo::Point::calcBoundingSphere(pts);
        checkSphere(bs, pts);
        fbs = Geo::Point_<float>::calcBoundingSphere(fpts);
        checkSphere(fbs, fpts);
        as = Geo::Point::calcApproxBoundingSphere(pts);
        checkSphere(as, pts);
        fas = Geo::Point_<float>::calcApproxBoundingSphere(fpts);
        checkSphere(fas, fpts);

        EXPECT_NEAR(bs.getRadius(), radius, tol);
        EXPECT_NEAR(fbs.getRadius(), radius, ftol);
        EXPECT_NEAR(as.getRadius(), radius, tol);
        EXPECT_NEAR(fas.getRadius(), radius, ftol);
    }
}

TEST(GeoTest, CollocatedPoints) {
    Random::Uniform random(0, 1000);
    Array_<Vec3> pts;
    Array_<fVec3> fpts;
    for (int trial = 0; trial < 1000; ++trial) {
        pts.clear();
        fpts.clear();
        int numPoints = random.getIntValue() + 1;
        Vec3 offs = SimTK::Test::randDouble() * 1000 * SimTK::Test::randVec3();
        fVec3 foffs((float)offs[0], (float)offs[1], (float)offs[2]);
        Real scale = offs.norm();

        for (int p = 0; p < numPoints; ++p) {
            Vec3 pt(offs);
            fVec3 fpt((float)pt[0], (float)pt[1], (float)pt[2]);
            pts.push_back(pt);
            fpts.push_back(fpt);
        }

        Geo::Sphere bs = Geo::Point::calcBoundingSphere(pts);
        checkSphere(bs, pts);

        Geo::Sphere_<float> fbs = Geo::Point_<float>::calcBoundingSphere(fpts);
        checkSphere(fbs, fpts);

        Geo::Sphere as = Geo::Point::calcApproxBoundingSphere(pts);
        checkSphere(as, pts);

        Geo::Sphere_<float> fas = Geo::Point_<float>::calcApproxBoundingSphere(fpts);
        checkSphere(fas, fpts);

        EXPECT_LT(bs.getRadius(), SqrtEps);
        EXPECT_LT(fbs.getRadius(), NTraits<float>::getSqrtEps());
        EXPECT_LT(as.getRadius(), SqrtEps);
        EXPECT_LT(fas.getRadius(), NTraits<float>::getSqrtEps());

        // Repeat test with random noise added.
        for (int p = 0; p < numPoints; ++p) {
            Vec3 noise(SimTK::Test::randVec3());
            fVec3 fnoise((float)noise[0], (float)noise[1], (float)noise[2]);
            pts[p] += SignificantReal * noise;
            fpts[p] += NTraits<float>::getSignificant() * fnoise;
        }

        bs = Geo::Point::calcBoundingSphere(pts);
        checkSphere(bs, pts);

        fbs = Geo::Point_<float>::calcBoundingSphere(fpts);
        checkSphere(fbs, fpts);

        as = Geo::Point::calcApproxBoundingSphere(pts);
        checkSphere(as, pts);

        fas = Geo::Point_<float>::calcApproxBoundingSphere(fpts);
        checkSphere(fas, fpts);
    }
}