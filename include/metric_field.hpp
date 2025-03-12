#ifndef METRICFIELD_H
#define METRICFIELD_H

#include <Eigen/Dense>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>
#include <mesh_element.hpp>

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

namespace MetricField {


using Quaternion = Eigen::Quaterniond;
using Vec3d = Eigen::Matrix<double, 3, 1>;
using Vec9d = Eigen::Matrix<double, 9, 1>;
using Mat3d = Eigen::Matrix<double, 3, 3>;
using Mat4d = Eigen::Matrix<double, 4, 4>;
using Mat9d = Eigen::Matrix<double, 9, 9>;

namespace OVM = OpenVolumeMesh;

using TM = OVM::TetrahedralGeometryKernel<OVM::Vec3d,
                                          OVM::TetrahedralMeshTopologyKernel>;


class MetricField {

private:
  TM tetmesh;

  // persistent mesh properties
  OVM::CellPropertyT<std::vector<OVM::VertexHandle>> cell_vertices;
  OVM::VertexPropertyT<Mat3d> metric;
  OVM::CellPropertyT<Mat3d> curl;
  OVM::CellPropertyT<Mat4d> T_prop;

  void prepareTransformation();
  void attachMetric();
  void attachCurl();
  void attach_cell_vertices();
  void calculateCurl(const OVM::CellHandle &ch);
  Vec3d integrate2Point(const OVM::CellHandle &ch, const Vec3d &a, const Vec3d &b);
  Mat3d computeCoeff(const OVM::CellHandle &ch, const Vec3d &q, const Vec3d &p);
  Mat3d recursiveDivide(const OVM::CellHandle &ch, const Vec3d &q,
                        const Vec3d &p);
  
  Mat9d scaledInverse(const Mat3d &A);

  // tet finder
  void readInTetmesh(const std::string& filename);
  Vec3d intersection(const Vec3d &_q, const Vec3d &_p, OVM::VertexHandle &_u,
                     OVM::VertexHandle &_v, OVM::VertexHandle &_w);
  bool pointInTet(const OVM::CellHandle &ch, const Vec3d &q);
  double orient3dHelper(const OVM::VertexHandle a, const OVM::VertexHandle b,
                        const OVM::VertexHandle c, const Vec3d d);

  double orient3dHelper(const OVM::Vec3d a, const OVM::Vec3d b,
                        const OVM::Vec3d c, const OVM::Vec3d d);

  double orient3dHelper(const OVM::HalfFaceHandle hfh,
                                   const Vec3d d);

  double orient3dHelper(const OVM::VertexHandle a, const OVM::VertexHandle b,
                        const Vec3d c, const Vec3d d);
  double orient3dHelper(const OVM::VertexHandle a, const OVM::VertexHandle b,
                        const OVM::VertexHandle c, const OVM::VertexHandle d);

  
  OVM::HalfFaceHandle findOppositeNextHalfFace(const OVM::CellHandle& ch, const OVM::HalfFaceHandle& hfh, const Vec3d &q, const Vec3d &p);
  bool pointInTetWithElement(OVM::CellHandle ch, const OVM::Vec3d param, MeshElement& elem);
  OVM::CellHandle findCorrectStartTet(OVM::CellHandle start_tet, const Vec3d &q, const Vec3d &p);

  std::vector<OVM::HalfFaceHandle> getIncidentHalfFaces(const OVM::CellHandle& ch, const OVM::EdgeHandle& eh);
  std::vector<OVM::HalfFaceHandle> getIncidentHalfFaces(const OVM::CellHandle& ch, const OVM::VertexHandle& vh);

  public:
  
  MetricField();
  MetricField(std::string filename); // filename of mesh
  MetricField(TM &tetmesh); // supply tetmesh
  Mat3d metricAtPoint(const OVM::CellHandle &ch, const Vec3d &_p);
  Mat3d eval_W(const OVM::CellHandle &ch, const Vec3d &_p);
  Mat3d getCurl(const OVM::CellHandle &ch);

  Mat3d computeCoeff(const Vec3d &q, const Vec3d &p, std::vector<OVM::CellHandle> &cells);

  // computes the rotation between two points in the metric field, main function, start_tet is the cell that contains q
  Mat3d computeCoeffImproved(const OVM::CellHandle &start_tet, const Vec3d &q, const Vec3d &p);

  TM& get_tetmesh();

  OVM::CellHandle startCellFinder(OVM::CellHandle start_tet, const Vec3d &q, const Vec3d &p);
  std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> tetFinderTobias(OVM::CellHandle start_tet, const Vec3d &q, const Vec3d &p, bool compute_intersections = true);

};


void hello_world();


// helper finder rotation coefficient
Mat3d unstack(const Vec9d &_m);
Vec9d stack(const Mat3d &_v);
Mat3d Ai_x(const Vec3d &A_i);
Mat9d A_x(const Mat3d &A);
Mat3d lie_exp(Vec3d u);

// helper functions tet finder
OVM::Vec3d toOVMVec3d(const Vec3d &_q);
// OVM::Vec3d toOVMVec3d(Vec3d &_q);
Vec3d toVec3d(const OVM::Vec3d &_q);
bool contains(const std::vector<OVM::CellHandle>& v, OVM::CellHandle x);
double orient3dPaper(const double *a, const double *b, const double *c,
                     const double *d);
bool sameSign(double a, double b);

// visualizing
void markTetsAsIntersected(TM &tetmesh, std::vector<OVM::CellHandle> tets);
void showIntersectionPoints(TM &tetmesh, std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> ints);

}
#endif