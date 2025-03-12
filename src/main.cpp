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


namespace MF = MetricField;


int main(int argc, char* argv[])
{
	std::cout << "hello wdsforld" << std::endl;
    // ovmTest();
    
    // showPoint();
    // robustRayTriangleTest();
    TM tetmesh;
    MF::populateMesh(tetmesh);
    MF::MetricField field = MF::MetricField(tetmesh);
    MF::saveToFile(tetmesh, "testmesh.ovm");

    // MF::normalCase(field);
    MF::normalTetFinderCase(field);
    // MF::tetFinderCaseWeird(field);
    MF::throughVertex(field);
    MF::throughEdge(field);
    MF::throughFace(field);
    MF::startPointOnVertex(field);
    MF::startPointOnEdge(field);
    MF::startPointOnFace(field);
    MF::endPointOnVertex(field);
    MF::endPointOnEdge(field);
    MF::endPointOnFace(field);
    MF::tetFinderSameStartEndPoint(field);
    MF::anotherCoolTest(field);

    MF::failingTestWithBoundingBox();
    
    // saveToFile(tetmesh, "testmesh.ovm");

    // throughVertex(tetmesh);
    
    // saveToFile(tetmesh, "oneTet.ovm");
    // deformedTestCase();

    // Vec3d a = Vec3d(-1,-1, 0.5);
    // Vec3d b = Vec3d(1,1,0);
    // Vec3d c = Vec3d(0,0, 1);
    // Vec3d d = Vec3d(0,0,0.4);
    // Vec3d e = Vec3d(-1,0,0);
    // Vec3d f = Vec3d(0,-1,0);

    // std::cout << " cdab " << orient3dPaper(c.data(),d.data(),a.data(),b.data()) << std::endl;
    // std::cout << " edab " << orient3dPaper(e.data(),d.data(),a.data(),b.data()) << std::endl;
    // std::cout << " fdab " << orient3dPaper(f.data(),d.data(),a.data(),b.data()) << std::endl;

    // std::cout << " ab " << orient3dPaper(f.data(),d.data(),a.data(),b.data()) << std::endl;
}

