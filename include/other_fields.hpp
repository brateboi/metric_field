#ifndef OTHERFIELD_H
#define OTHERFIELD_H

#include <Eigen/Dense>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>


#include <cmath>
#include <limits>

namespace MetricField {

using Quaternion = Eigen::Quaterniond;
using Vec3d = Eigen::Matrix<double, 3, 1>;
using Vec4d = Eigen::Matrix<double, 4, 1>;
using Vec9d = Eigen::Matrix<double, 9, 1>;
using Mat3d = Eigen::Matrix<double, 3, 3>;
using Mat4d = Eigen::Matrix<double, 4, 4>;
using Mat9d = Eigen::Matrix<double, 9, 9>;

namespace OVM = OpenVolumeMesh;
using TM = OVM::TetrahedralGeometryKernel<OVM::Vec3d,
                                          OVM::TetrahedralMeshTopologyKernel>;

void saveToFile(TM &tetmesh, std::string filename);
void saveToFileOVMB(TM &tetmesh, std::string filename);
void readFile(TM &mesh, std::string filename);

Mat3d rootMetricAnalytical(Vec3d &p);
Vec9d curlAnalytical(Vec3d p);
Mat3d eval_W_Analytical(Vec3d _p);
Vec3d integrate2PointAnalytical(Vec3d a, Vec3d b);
Mat3d recursiveDivideAnalytical(const Vec3d &a, const Vec3d &b, int depth);
Mat3d computeCoeffAnalytical(const Vec3d &q, const Vec3d &p, int maxdepth);
void attachConstantMetric(TM &tetmesh);
Mat3d rotationPushforward(Vec3d &p1, Vec3d &p2);
Mat3d pushforward(Vec3d &p);
Vec3d deformInverse(Vec3d &p);

// constant discretization
Mat3d computeCoeffConstant(const OVM::CellHandle a, const OVM::CellHandle b, TM &tetmesh);
OVM::FaceHandle commonFace(const OVM::CellHandle a, const OVM::CellHandle b, TM &tetmesh);
void calculate_curl_constant_and_W(const OVM::CellHandle a, const OVM::CellHandle b, TM &tetmesh, Mat3d &curl_constant, Mat3d &W_a, Mat3d &W_b);

} // end namespace MetricField

#endif