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
#include <mutex>

using Quaternion = Eigen::Quaterniond;
using Vec3d = Eigen::Matrix<double, 3, 1>;
using Vec9d = Eigen::Matrix<double, 9, 1>;
using Mat3d = Eigen::Matrix<double, 3, 3>;
using Mat9d = Eigen::Matrix<double, 9, 9>;

namespace OVM = OpenVolumeMesh;

using TM = OVM::TetrahedralGeometryKernel<OVM::Vec3d,
                                          OVM::TetrahedralMeshTopologyKernel>;

struct cmpVec3d {
  bool operator()(const Vec3d &a, const Vec3d &b) const {
    return std::lexicographical_compare(a.data(), a.data() + a.size(), b.data(),
                                        b.data() + b.size());
  }
};

struct hashFunc{
    size_t operator()(const Vec3d &k) const{
    size_t h1 = std::hash<double>()(k(0));
    size_t h2 = std::hash<double>()(k(1));
    size_t h3 = std::hash<double>()(k(2));
    return (h1 ^ (h2 << 1)) ^ h3;
    }
};



struct xor_hash {
  template <class T1, class T2>
    std::size_t operator() (const std::tuple<T1, T2> &v) const
    {
        return std::hash<T1>()(std::get<0>(v)) ^ (std::hash<T2>()(std::get<1>(v)) << 1);
    }
};

class MetricField {

private:
  TM tetmesh;
  // std::map<Eigen::VectorXd, std::vector<OVM::CellHandle>, cmpVec3d> hash;
  std::unordered_map<Eigen::VectorXd, std::vector<OVM::CellHandle>, hashFunc> hash;
  std::mutex hash_mutex;
  OVM::CellHandle prevTetContainingPoint = OVM::CellHandle(0);

  OVM::CellPropertyT<std::vector<OVM::VertexHandle>> cell_vertices;

  int degenerate_counter = 0;

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
  OVM::VertexHandle linearNN(const Vec3d &q);
  //std::vector<OVM::CellHandle> locateTets(const Vec3d &q);
  bool initialization(OVM::CellHandle &t, const Vec3d &q, const Vec3d &p,
                      OVM::VertexHandle &u, OVM::VertexHandle &v,
                      OVM::VertexHandle &w);

  OVM::HalfFaceHandle
  find_halfface_in_cell(std::vector<OVM::VertexHandle> &verts,
                        const OVM::CellHandle &ch);

  bool validState(const Vec3d &q, const Vec3d &p, OVM::VertexHandle u,
                  OVM::VertexHandle v, OVM::VertexHandle w);

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


  void initializeFast(const OVM::VertexHandle q, const Vec3d p, OVM::VertexHandle &u, OVM::VertexHandle &v, OVM::VertexHandle &w, OVM::CellHandle &t);

  OVM::CellHandle neighbor(const OVM::VertexHandle &u, const OVM::VertexHandle &v, const OVM::VertexHandle &w, const OVM::CellHandle &t);
  OVM::VertexHandle vertex_of_t(const OVM::VertexHandle &u, const OVM::VertexHandle &v, const OVM::VertexHandle &w, const OVM::CellHandle &t);

  std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> tetFinderFast(const Vec3d &q, const Vec3d &p, bool &failed);
  std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> tetFinderRobust(const Vec3d &q, const Vec3d &p);


  
  OVM::HalfFaceHandle findOppositeNextHalfFace(const OVM::CellHandle& ch, const OVM::HalfFaceHandle& hfh, const Vec3d &q, const Vec3d &p);
  bool pointInTetWithElement(OVM::CellHandle ch, const OVM::Vec3d param, MeshElement& elem);
  OVM::CellHandle findCorrectStartTet(OVM::CellHandle start_tet, const Vec3d &q, const Vec3d &p);

  std::vector<OVM::HalfFaceHandle> getIncidentHalfFaces(const OVM::CellHandle& ch, const OVM::EdgeHandle& eh);
  std::vector<OVM::HalfFaceHandle> getIncidentHalfFaces(const OVM::CellHandle& ch, const OVM::VertexHandle& vh);

  public:
  // std::unordered_map<u_int64_t, Mat3d> metricHash; // hashes the metric at the midpoint of the node_code
  std::map<u_int64_t, Mat3d> metricHash; // hashes the metric at the midpoint of the node_code
  // std::unordered_map<std::tuple<u_int64_t, uint64_t>, Vec3d, xor_hash> rotationHash; // hashes the rotations as EulerZYX angles, always sorted in lower<higher
  std::map<std::tuple<u_int64_t, uint64_t>, Vec3d> rotationHash; // hashes the rotations as EulerZYX angles, always sorted in lower<higher
  // std::unordered_map<u_int64_t, std::map<u_int64_t, Vec3d>> rotationHashNested; // hashes the rotations as EulerZYX angles, always sorted in lower<higher


  Mat3d computeCoeff(const Vec3d& v_pos, const Vec3d& n_pos);
  MetricField();
  MetricField(std::string filename); // filename of mesh
  MetricField(TM &tetmesh); // supply tetmesh
  Mat3d metricAtPoint(const OVM::CellHandle &ch, const Vec3d &_p);
  Mat3d eval_W(const OVM::CellHandle &ch, const Vec3d &_p);
  Mat3d getCurl(const OVM::CellHandle &ch);
  OVM::CellHandle locateTetsFast(const OVM::VertexHandle &q, const Vec3d &p, bool &failed);
  std::vector<OVM::CellHandle> locateTets(const Vec3d &q);

  Mat3d computeCoeff(const Vec3d &q, const Vec3d &p, std::vector<OVM::CellHandle> &cells);

  Mat3d computeCoeff(const Vec3d &q, const Vec3d &p, OVM::CellHandle &cq, OVM::CellHandle &cp);

  // computes the rotation between two points in the metric field, main function, start_tet is the cell that contains q
  Mat3d computeCoeffImproved(const OVM::CellHandle &start_tet, const Vec3d &q, const Vec3d &p);

  int get_degenerate_counter();
  
  std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> tetFinder(const Vec3d &q, const Vec3d &p);

  TM& get_tetmesh();

  Mat3d metricAtPoint(const Vec3d &_p); // metric accessible globably

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
bool robustRayTriangle(const Vec3d& q, const Vec3d& p, const OVM::Vec3d& v1,
                       const OVM::Vec3d& v2, const OVM::Vec3d& v3);

// visualizing
void markTetsAsIntersected(TM &tetmesh, std::vector<OVM::CellHandle> tets);
void showIntersectionPoints(TM &tetmesh, std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> ints);
#endif