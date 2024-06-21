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





int main(int argc, char* argv[])
{
    
	std::cout << "hello wdsforld" << std::endl;
    // ovmTest();
    
    showPoint();
    // robustRayTriangleTest();
    TM tetmesh;
    populateMesh(tetmesh);
    MetricField field = MetricField(tetmesh);

    // normalCase(field);
    // throughVertex(field);
    // throughEdge(field);
    // throughFace(field);
    // startPointOnVertex(field);
    // startPointOnEdge(field);
    // startPointOnFace(field);
    // endPointOnVertex(field);
    // endPointOnEdge(field);
    // endPointOnFace(field);
    saveToFile(tetmesh, "testmesh.ovm");

    // throughVertex(tetmesh);
    
    // saveToFile(tetmesh, "oneTet.ovm");
    // deformedTestCase();


}

