#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>

#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>
#include <OpenVolumeMesh/FileManager/FileManager.hh>

#include <Eigen/Dense>
#include <typeinfo>
#include "math.h"
#include "tests.h"
#include "other_fields.hpp"

#include <iostream>

using Quaternion = Eigen::Quaterniond;
using Vec3d = Eigen::Matrix<double, 3, 1>;
using Vec4d = Eigen::Matrix<double, 4, 1>;
using Vec9d = Eigen::Matrix<double, 9, 1>;
using Mat3d = Eigen::Matrix<double, 3, 3>;
using Mat4d = Eigen::Matrix<double, 4, 4>;
using Mat9d = Eigen::Matrix<double, 9, 9>;

using TM = OVM::TetrahedralGeometryKernel<OVM::Vec3d,
                                          OVM::TetrahedralMeshTopologyKernel>;

/*
This sets up different meshes to test edge cases
 */

/**
 * adds tet with the 4 vertices to the mesh
 */
void addTet(TM &tetmesh, OVM::VertexHandle &v0, OVM::VertexHandle &v1, OVM::VertexHandle &v2, OVM::VertexHandle &v3)
{
    tetmesh.add_cell(v0, v1, v2, v3, true);
}

// clears tetmesh and populates it
void populateMesh(TM &tetmesh)
{
    tetmesh.clear();
    OVM::VertexHandle v0 = tetmesh.add_vertex(OVM::Vec3d(0.0, 0.0, 0.0));
    OVM::VertexHandle v1 = tetmesh.add_vertex(OVM::Vec3d(1.0, 0.0, 0.0));
    OVM::VertexHandle v2 = tetmesh.add_vertex(OVM::Vec3d(0.0, 1.0, 0.0));
    OVM::VertexHandle v3 = tetmesh.add_vertex(OVM::Vec3d(0.0, 0.0, 1.0));
    OVM::VertexHandle v4 = tetmesh.add_vertex(OVM::Vec3d(-1.0, 0.0, 0.0));
    OVM::VertexHandle v5 = tetmesh.add_vertex(OVM::Vec3d(0.0, -1.0, 0.0));
    OVM::VertexHandle v6 = tetmesh.add_vertex(OVM::Vec3d(0.0, 0.0, -1.0));
    OVM::VertexHandle a1 = tetmesh.add_vertex(OVM::Vec3d(2.0, 2.0, 2.0));
    OVM::VertexHandle a2 = tetmesh.add_vertex(OVM::Vec3d(2.0, 2.0, -2.0));
    OVM::VertexHandle a3 = tetmesh.add_vertex(OVM::Vec3d(2.0, -2.0, 2.0));
    OVM::VertexHandle a4 = tetmesh.add_vertex(OVM::Vec3d(2.0, -2.0, -2.0));
    OVM::VertexHandle a5 = tetmesh.add_vertex(OVM::Vec3d(-2.0, 2.0, 2.0));
    OVM::VertexHandle a6 = tetmesh.add_vertex(OVM::Vec3d(-2.0, 2.0, -2.0));
    OVM::VertexHandle a7 = tetmesh.add_vertex(OVM::Vec3d(-2.0, -2.0, 2.0));
    OVM::VertexHandle a8 = tetmesh.add_vertex(OVM::Vec3d(-2.0, -2.0, -2.0));
    // int i = 80;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v0, v1, v2, v3) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v0, v1, v2, v3) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v6, v1, v2, v0) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v0, v5, v1, v3) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v0, v3, v2, v4) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v0, v4, v5, v3) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v6, v2, v4, v0) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v6, v4, v5, v0) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v0, v1, v5, v6) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v1, v2, v3, a1) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v2, v1, v6, a2) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v5, v1, v3, a3) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v1, v5, v6, a4) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v2, v4, v3, a5) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v4, v2, v6, a6) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v4, v5, v3, a7) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v5, v4, v6, a8) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v1, v2, a1, a2) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v3, v1, a1, a3) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v5, v1, a3, a4) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v1, v6, a2, a4) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v2, v4, a5, a6) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v4, v3, a5, a7) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v4, v5, a7, a8) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v6, v4, a6, a8) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v2, v3, a1, a5) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v6, v2, a2, a6) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v3, v5, a3, a7) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v5, v6, a4, a8) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v1, a2, a1, a3) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v1, a2, a3, a4) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v4, a5, a6, a7) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v4, a7, a6, a8) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v2, a1, a2, a5) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v2, a5, a2, a6) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v5, a4, a3, a7) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v5, a4, a7, a8) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v3, a3, a1, a5) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v3, a3, a5, a7) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v6, a2, a4, a6) << std::endl;
    // std::cout << i++ << "tet " << orient3dHelper(tetmesh, v6, a6, a4, a8) << std::endl;

    addTet(tetmesh, v0, v1, v2, v3);
    addTet(tetmesh, v6, v1, v2, v0);
    addTet(tetmesh, v0, v5, v1, v3);
    addTet(tetmesh, v0, v3, v2, v4);
    addTet(tetmesh, v0, v4, v5, v3);
    addTet(tetmesh, v6, v2, v4, v0);
    addTet(tetmesh, v6, v4, v5, v0);
    addTet(tetmesh, v0, v1, v5, v6);
    addTet(tetmesh, v1, v2, v3, a1);
    addTet(tetmesh, v2, v1, v6, a2);
    addTet(tetmesh, v5, v1, v3, a3);
    addTet(tetmesh, v1, v5, v6, a4);
    addTet(tetmesh, v2, v4, v3, a5);
    addTet(tetmesh, v4, v2, v6, a6);
    addTet(tetmesh, v4, v5, v3, a7);
    addTet(tetmesh, v5, v4, v6, a8);
    addTet(tetmesh, v1, v2, a1, a2);
    addTet(tetmesh, v3, v1, a1, a3);
    addTet(tetmesh, v5, v1, a3, a4);
    addTet(tetmesh, v1, v6, a2, a4);
    addTet(tetmesh, v2, v4, a5, a6);
    addTet(tetmesh, v4, v3, a5, a7);
    addTet(tetmesh, v4, v5, a7, a8);
    addTet(tetmesh, v6, v4, a6, a8);
    addTet(tetmesh, v2, v3, a1, a5);
    addTet(tetmesh, v6, v2, a2, a6);
    addTet(tetmesh, v3, v5, a3, a7);
    addTet(tetmesh, v5, v6, a4, a8);
    addTet(tetmesh, v1, a2, a1, a3);
    addTet(tetmesh, v1, a2, a3, a4);
    addTet(tetmesh, v4, a5, a6, a7);
    addTet(tetmesh, v4, a7, a6, a8);
    addTet(tetmesh, v2, a1, a2, a5);
    addTet(tetmesh, v2, a5, a2, a6);
    addTet(tetmesh, v5, a4, a3, a7);
    addTet(tetmesh, v5, a4, a7, a8);
    addTet(tetmesh, v3, a3, a1, a5);
    addTet(tetmesh, v3, a3, a5, a7);
    addTet(tetmesh, v6, a2, a4, a6);
    addTet(tetmesh, v6, a6, a4, a8);
}

void normalTetFinderCase(MetricField &field)
{
    OVM::Vec3d p = OVM::Vec3d(0.25, 0.87, 1.22);
    OVM::Vec3d q = OVM::Vec3d(-0.225, -0.425, 1.25);
    std::cout << "----------------- START normalTetFinderCase Case ------------------" << std::endl;
    // auto res = field.tetFinder(toVec3d(q), toVec3d(p));

    auto start_cell = field.startCellFinder(OVM::CellHandle(30), toVec3d(field.get_tetmesh().barycenter(OVM::CellHandle(30))),toVec3d(q));
    auto res = field.tetFinderTobias(start_cell, toVec3d(q), toVec3d(p));
    visualizePath(field, res);
    std::cout << "------------------ END normalTetFinderCase Case -------------------" << std::endl;
}

void throughVertex(MetricField &field)
{
    OVM::Vec3d q = OVM::Vec3d(0.25, 0.25, 0.25);
    OVM::Vec3d p = OVM::Vec3d(-0.25, -0.25, -0.25);
    std::cout << "-----------------START Through Vertex ------------------" << std::endl;
    // auto res = field.tetFinder(toVec3d(q), toVec3d(p));

    auto start_cell = field.startCellFinder(OVM::CellHandle(30), toVec3d(field.get_tetmesh().barycenter(OVM::CellHandle(30))),toVec3d(q));
    auto res = field.tetFinderTobias(start_cell, toVec3d(q), toVec3d(p));
    visualizePath(field, res);
    std::cout << "-----------------END Through Vertex ------------------" << std::endl;
}

void throughEdge(MetricField &field)
{
    OVM::Vec3d q = OVM::Vec3d(1.5, -0.5, 0.0);
    OVM::Vec3d p = OVM::Vec3d(-0.5, 1.5, 0.0);
    std::cout << "-----------------START Through Edge ------------------" << std::endl;
    // auto res = field.tetFinder(toVec3d(q), toVec3d(p));

    auto start_cell = field.startCellFinder(OVM::CellHandle(30), toVec3d(field.get_tetmesh().barycenter(OVM::CellHandle(30))),toVec3d(q));
    auto res = field.tetFinderTobias(start_cell, toVec3d(q), toVec3d(p));
    visualizePath(field, res);
    std::cout << "-----------------END Through Edge ------------------" << std::endl;
}

void throughFace(MetricField &field)
{
    OVM::Vec3d q = OVM::Vec3d(1.0, 1.0, 0.0);
    OVM::Vec3d p = OVM::Vec3d(0.75, -0.75, 0.0);
    std::cout << "-----------------START Through Face ------------------" << std::endl;
    // auto res = field.tetFinder(toVec3d(q), toVec3d(p));

    auto start_cell = field.startCellFinder(OVM::CellHandle(30), toVec3d(field.get_tetmesh().barycenter(OVM::CellHandle(30))),toVec3d(q));
    auto res = field.tetFinderTobias(start_cell, toVec3d(q), toVec3d(p));
    visualizePath(field, res);
    std::cout << "-----------------END Through Face ------------------" << std::endl;
}

void startPointOnVertex(MetricField &field)
{
    OVM::Vec3d q = OVM::Vec3d(1.0, 0.0, 0.0);
    OVM::Vec3d p = OVM::Vec3d(0.25, 0.87, 1.22);
    std::cout << "-----------------START Start on Vertex ------------------" << std::endl;
    // auto res = field.tetFinder(toVec3d(q), toVec3d(p));

    auto start_cell = field.startCellFinder(OVM::CellHandle(30), toVec3d(field.get_tetmesh().barycenter(OVM::CellHandle(30))),toVec3d(q));
    auto res = field.tetFinderTobias(start_cell, toVec3d(q), toVec3d(p));
    visualizePath(field, res);
    std::cout << "-----------------END Start on Vertex ------------------" << std::endl;
}

void startPointOnEdge(MetricField &field)
{
    OVM::Vec3d q = OVM::Vec3d(0.5, 0.5, 0.0);
    OVM::Vec3d p = OVM::Vec3d(0.25, 0.87, 1.22);
    std::cout << "-----------------START Start on Edge ------------------" << std::endl;
    // auto res = field.tetFinder(toVec3d(q), toVec3d(p));

    auto start_cell = field.startCellFinder(OVM::CellHandle(30), toVec3d(field.get_tetmesh().barycenter(OVM::CellHandle(30))),toVec3d(q));
    auto res = field.tetFinderTobias(start_cell, toVec3d(q), toVec3d(p));
    visualizePath(field, res);
    std::cout << "-----------------END Start on Edge ------------------" << std::endl;
}

void startPointOnFace(MetricField &field)
{
    OVM::Vec3d q = OVM::Vec3d(0.25, 0.25, 0.0);
    OVM::Vec3d p = OVM::Vec3d(0.25, 0.87, 1.22);
    std::cout << "-----------------START Start on Face ------------------" << std::endl;
    // auto res = field.tetFinder(toVec3d(q), toVec3d(p));

    auto start_cell = field.startCellFinder(OVM::CellHandle(30), toVec3d(field.get_tetmesh().barycenter(OVM::CellHandle(30))),toVec3d(q));
    auto res = field.tetFinderTobias(start_cell, toVec3d(q), toVec3d(p));
    visualizePath(field, res);
    std::cout << "-----------------END Start on Face ------------------" << std::endl;
}

void endPointOnVertex(MetricField &field)
{
    OVM::Vec3d p = OVM::Vec3d(1.0, 0.0, 0.0);
    OVM::Vec3d q = OVM::Vec3d(0.25, 0.87, 1.22);
    std::cout << "-----------------START End on Vertex ------------------" << std::endl;
    // auto res = field.tetFinder(toVec3d(q), toVec3d(p));

    auto start_cell = field.startCellFinder(OVM::CellHandle(30), toVec3d(field.get_tetmesh().barycenter(OVM::CellHandle(30))),toVec3d(q));
    auto res = field.tetFinderTobias(start_cell, toVec3d(q), toVec3d(p));
    visualizePath(field, res);
    std::cout << "-----------------END End on Vertex------------------" << std::endl;
}

void endPointOnEdge(MetricField &field)
{
    OVM::Vec3d p = OVM::Vec3d(0.5, 0.5, 0.0);
    OVM::Vec3d q = OVM::Vec3d(0.25, 0.87, 1.22);
    std::cout << "-----------------START End on Edge ------------------" << std::endl;
    // auto res = field.tetFinder(toVec3d(q), toVec3d(p));

    auto start_cell = field.startCellFinder(OVM::CellHandle(30), toVec3d(field.get_tetmesh().barycenter(OVM::CellHandle(30))),toVec3d(q));
    auto res = field.tetFinderTobias(start_cell, toVec3d(q), toVec3d(p));
    visualizePath(field, res);
    std::cout << "-----------------END End on Edge ------------------" << std::endl;
}

void endPointOnFace(MetricField &field)
{
    OVM::Vec3d p = OVM::Vec3d(0.25, 0.25, 0.0);
    OVM::Vec3d q = OVM::Vec3d(0.25, 0.87, 1.22);
    std::cout << "-----------------START End on Face ------------------" << std::endl;
    // auto res = field.tetFinder(toVec3d(q), toVec3d(p));

    auto start_cell = field.startCellFinder(OVM::CellHandle(30), toVec3d(field.get_tetmesh().barycenter(OVM::CellHandle(30))),toVec3d(q));
    auto res = field.tetFinderTobias(start_cell, toVec3d(q), toVec3d(p));
    visualizePath(field, res);
    std::cout << "-----------------END End on Face ------------------" << std::endl;
}

void tetFinderSameStartEndPoint(MetricField &field){

    auto p= OVM::Vec3d(0.75, 0.75, 0.75);
    auto q= OVM::Vec3d(0.75, 0.75, 0.75);

    
    std::cout << "-----------------START SameStartEndPoint ------------------" << std::endl;
    // auto res = field.tetFinder(toVec3d(q), toVec3d(p));

    auto start_cell = field.startCellFinder(OVM::CellHandle(30), toVec3d(field.get_tetmesh().barycenter(OVM::CellHandle(30))),toVec3d(q));
    auto res = field.tetFinderTobias(start_cell, toVec3d(q), toVec3d(p));
    visualizePath(field, res);
    std::cout << "-----------------END SameStartEndPoint ------------------" << std::endl;
}

void anotherCoolTest(MetricField &field){

    auto p= OVM::Vec3d(0, 1.75, -1.75);
    auto q= OVM::Vec3d(0, 1.75, 1.75);

    
    std::cout << "-----------------START anotherCoolTest ------------------" << std::endl;
    // auto res = field.tetFinder(toVec3d(q), toVec3d(p));

    auto start_cell = field.startCellFinder(OVM::CellHandle(30), toVec3d(field.get_tetmesh().barycenter(OVM::CellHandle(30))),toVec3d(q));
    auto res = field.tetFinderTobias(start_cell, toVec3d(q), toVec3d(p));
    visualizePath(field, res);
    std::cout << "-----------------END anotherCoolTest ------------------" << std::endl;
}

void visualizePath(MetricField &field, std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> lineSegments){

    TM visualized_path;

    TM ground_mesh = field.get_tetmesh();


    auto addCell = [&ground_mesh, &visualized_path](OVM::CellHandle cell, OVM::Vec3d p1, OVM::Vec3d p2){
        auto verts = ground_mesh.get_cell_vertices(cell);
        auto v1 = visualized_path.add_vertex(ground_mesh.vertex(verts[0]));
        auto v2 = visualized_path.add_vertex(ground_mesh.vertex(verts[1]));
        auto v3 = visualized_path.add_vertex(ground_mesh.vertex(verts[2]));
        auto v4 = visualized_path.add_vertex(ground_mesh.vertex(verts[3]));
        visualized_path.add_cell(v1,v2,v3,v4);

        visualized_path.add_edge(visualized_path.add_vertex(p1), visualized_path.add_vertex(p2));

    };

    for (const auto &segment : lineSegments)
    {
        auto cell = std::get<0>(segment);
        OVM::Vec3d point1 = toOVMVec3d(std::get<1>(segment));
        OVM::Vec3d point2 = toOVMVec3d(std::get<2>(segment));
        addCell(cell, point1, point2);
    }

    saveToFile(visualized_path, "visualized_path.ovm");
}


void showPoint()
{
    TM point;
    point.add_vertex(OVM::Vec3d(0.25, 0.87, 1.22));
    saveToFile(point, "point.ovm");
}
