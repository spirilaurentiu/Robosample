#include <gtest/gtest.h>
#include <sstream>
#include <string>

#include "SimTKcommon.h"

using namespace SimTK;

/**
 * @test SimTKCommon_PolygonalMesh_ManualConstruction, CorrectlyTracksTopologyAndGeometry
 * @details Validates the basic API for building a mesh by adding vertices and faces manually.
 */
TEST(SimTKCommon_PolygonalMesh_ManualConstruction, CorrectlyTracksTopologyAndGeometry) {
    PolygonalMesh mesh;

    EXPECT_EQ(mesh.getNumFaces(), 0);
    EXPECT_EQ(mesh.getNumVertices(), 0);

    // Verify vertex addition returns sequential indices
    EXPECT_EQ(mesh.addVertex(Vec3(0)), 0);
    EXPECT_EQ(mesh.addVertex(Vec3(0, 1, 0)), 1);
    EXPECT_EQ(mesh.addVertex(Vec3(0, 0, 1)), 2);

    // Build a face from the added vertices
    const Array_<int> vertices({1, 2, 0});
    EXPECT_EQ(mesh.addFace(vertices), 0);

    // Validate Connectivity
    EXPECT_EQ(mesh.getNumFaces(), 1);
    EXPECT_EQ(mesh.getNumVertices(), 3);
    EXPECT_EQ(mesh.getFaceVertex(0, 0), 1);
    EXPECT_EQ(mesh.getFaceVertex(0, 1), 2);
    EXPECT_EQ(mesh.getFaceVertex(0, 2), 0);

    // Validate Spatial Positions
    EXPECT_EQ(mesh.getVertexPosition(0), Vec3(0));
    EXPECT_EQ(mesh.getVertexPosition(1), Vec3(0, 1, 0));
    EXPECT_EQ(mesh.getVertexPosition(2), Vec3(0, 0, 1));
}

/**
 * @test SimTKCommon_PolygonalMesh_HandleSemantics, DistinguishesShallowAndDeepCopies
 * @details Ensures Simbody's internal reference counting and deep copy logic are preserved.
 */
TEST(SimTKCommon_PolygonalMesh_HandleSemantics, DistinguishesShallowAndDeepCopies) {
    PolygonalMesh mesh;
    static_cast<void>(mesh.addVertex(Vec3(0)));

    // Make sure copy and assignment are shallow.
    // Shallow copy construction: both handles point to the same implementation
    const PolygonalMesh mesh2(mesh);
    EXPECT_EQ(&mesh2.getImpl(), &mesh.getImpl());
    EXPECT_EQ(mesh.getImplHandleCount(), 2);

    PolygonalMesh mesh3;
    mesh3 = mesh; // Shallow assignment
    EXPECT_EQ(&mesh3.getImpl(), &mesh.getImpl());
    EXPECT_EQ(mesh.getImplHandleCount(), 3);

    // Deep copy: data is identical, but memory addresses of the implementation differ
    PolygonalMesh mesh4;
    mesh4.copyAssign(mesh);
    EXPECT_EQ(mesh4.getNumVertices(), mesh.getNumVertices());
    EXPECT_NE(&mesh4.getImpl(), &mesh.getImpl());
    EXPECT_EQ(mesh4.getImplHandleCount(), 1);
    EXPECT_EQ(mesh.getImplHandleCount(), 3);
}

/**
 * @test SimTKCommon_PolygonalMesh_ObjParser, CorrectlyHandlesComplexObjFormatting
 * @details Tests the OBJ loader against edge cases like line continuations and relative indexing.
 */
TEST(SimTKCommon_PolygonalMesh_ObjParser, CorrectlyHandlesComplexObjFormatting) {
    std::string fileContent;
    fileContent += "# This is a comment\n";
    fileContent += "v -1.0 1.0 2.0\n";
    fileContent += "v -2.0 2.0 3.0\n";
    fileContent += "v -3.0 3.0 \\\n"; // Line continuation test
    fileContent += "4.0\n";
    fileContent += "v -4.0 4.0 5.0\n";
    fileContent += "v -5.0 5.0 6.0\n";
    fileContent += "f 1 2 3\n";          // Standard indexing
    fileContent += "f 2// 3/4/5 4//2\n"; // Indices with skipped attributes
    fileContent += "f -3 -2/3/4 -1\n";   // Relative indexing
    fileContent += "f 1 3\\\n";          // Face line continuation
    fileContent += "5 2\n";

    PolygonalMesh mesh;
    std::stringstream stream(fileContent);
    mesh.loadObjFile(stream);

    ASSERT_EQ(mesh.getNumVertices(), 5);
    ASSERT_EQ(mesh.getNumFaces(), 4);

    // Verify coordinates: v_i = [-(i+1), i+1, i+2]
    for (int i = 0; i < mesh.getNumVertices(); ++i) {
        const Vec3& pos = mesh.getVertexPosition(i);
        EXPECT_DOUBLE_EQ(pos[0], static_cast<double>(-(i + 1)));
        EXPECT_DOUBLE_EQ(pos[1], static_cast<double>(i + 1));
        EXPECT_DOUBLE_EQ(pos[2], static_cast<double>(i + 2));
    }

    // Verify triangular faces (first 3 faces)
    for (int i = 0; i < 3; ++i) {
        ASSERT_EQ(mesh.getNumVerticesForFace(i), 3);
        EXPECT_EQ(mesh.getFaceVertex(i, 0), i);
        EXPECT_EQ(mesh.getFaceVertex(i, 1), i + 1);
        EXPECT_EQ(mesh.getFaceVertex(i, 2), i + 2);
    }

    // Verify the 4-vertex face created via continuation and custom indexing
    ASSERT_EQ(mesh.getNumVerticesForFace(3), 4);
    EXPECT_EQ(mesh.getFaceVertex(3, 0), 0);
    EXPECT_EQ(mesh.getFaceVertex(3, 1), 2);
    EXPECT_EQ(mesh.getFaceVertex(3, 2), 4);
    EXPECT_EQ(mesh.getFaceVertex(3, 3), 1);
}
