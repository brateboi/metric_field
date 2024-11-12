#include <Eigen/Dense>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

#include "metric_field.hpp"
#include "other_fields.hpp"
#include "wrapper.h"

#include <algorithm>

#include <cmath>
#include <limits>

using Quaternion = Eigen::Quaterniond;
using Vec3d = Eigen::Matrix<double, 3, 1>;
using Vec4d = Eigen::Matrix<double, 4, 1>;
using Vec9d = Eigen::Matrix<double, 9, 1>;
using Mat3d = Eigen::Matrix<double, 3, 3>;
using Mat4d = Eigen::Matrix<double, 4, 4>;
using Mat9d = Eigen::Matrix<double, 9, 9>;

using TM = OVM::TetrahedralGeometryKernel<OVM::Vec3d,
                                          OVM::TetrahedralMeshTopologyKernel>;

#include <iostream>



void hello_world(){
  std::cout << "HELLO WORLD, C++ call from C worked! hopefull" << std::endl;
  return;
} 



// exact predicates
extern "C" double orient3d(const double *, const double *, const double *,
                           const double *);
extern "C" double orient2d(const double *, const double *, const double *);
extern "C" void exactinit();

void MetricField::readInTetmesh(const std::string &filename)
{
  OVM::IO::FileManager fm;
  fm.readFile(filename, tetmesh, true, true);
}

TM& MetricField::get_tetmesh(){
  return tetmesh;
}

// int MetricField::get_degenerate_counter(){
//   return degenerate_counter;
// }

// computes the rotation between two points in the metric field, main function
// Mat3d MetricField::computeCoeff(const Vec3d &q, const Vec3d &p)
// {
//   // return Mat3d::Identity(); // no work

//   std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> lineSegments = tetFinder(q, p);

//   Mat3d coeff = Mat3d::Identity();
//   for (const auto &segment : lineSegments)
//   {
//     coeff = computeCoeff(std::get<0>(segment), std::get<1>(segment), std::get<2>(segment)) * coeff;
//   }

//   return coeff;
// }

// computes the rotation between two points in the metric field, main function, start_tet is the cell that contains q
Mat3d MetricField::computeCoeffImproved(const OVM::CellHandle &start_tet, const Vec3d &q, const Vec3d &p)
{
  assert(pointInTet(start_tet, q)); // q needs to be contained in startTet, otherwise tetFinderTobias doesnt work
  // return Mat3d::Identity(); // no work

  std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> lineSegments = tetFinderTobias(start_tet, q, p);

  Mat3d coeff = Mat3d::Identity();
  for (const auto &segment : lineSegments)
  {
    coeff = computeCoeff(std::get<0>(segment), std::get<1>(segment), std::get<2>(segment)) * coeff;
  }

  return coeff;
}

// std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> MetricField::tetFinder(const Vec3d &q, const Vec3d &p){
//   bool failed = false;

//   // std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> lineSegments = tetFinderFast(q, p, failed);
//   std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> lineSegments = tetFinderRobust(q, p);

//   // check if run was degenerate, use robust, but slow tet finder
//   if (failed){
//     degenerate_counter++;
//     lineSegments = tetFinderRobust(q, p);
//   }
//   return lineSegments;
// }

// computes the rotation between two points in the metric field, main function, pushes traversed cells to vector
// Mat3d MetricField::computeCoeff(const Vec3d &q, const Vec3d &p, std::vector<OVM::CellHandle> &cells)
// {
//   std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> lineSegments = tetFinder(q, p);
  
//   Mat3d coeff = Mat3d::Identity();
//   for (const auto &segment : lineSegments)
//   {
//     // for visualizing
//     // std::cout << "Cell handle " << std::get<0>(segment) << std::endl;
//     cells.push_back(std::get<0>(segment));

//     coeff = computeCoeff(std::get<0>(segment), std::get<1>(segment), std::get<2>(segment)) * coeff;
//   }

//   // for visualizing---------
//   // markTetsAsIntersected(this->tetmesh, cells);
//   // showIntersectionPoints(this->tetmesh, lineSegments);
//   // ------------------------

//   return coeff;
// }

// computes the rotation between two points in the metric field, main function, but with info in which CellHandle p/q are
// Mat3d MetricField::computeCoeff(const Vec3d &q, const Vec3d &p, OVM::CellHandle &cq, OVM::CellHandle &cp)
// {
//   // fill hash with tet info
//   std::vector<OVM::CellHandle> tets;
//   tets.push_back(cq);
//   #pragma omp critical
//   {
//     hash.try_emplace(q, tets);

//     tets.clear();
//     tets.push_back(cp);
//     hash.try_emplace(p, tets);
//   }

//   std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> lineSegments = tetFinder(q, p);
//   Mat3d coeff = Mat3d::Identity();
//   for (const auto &segment : lineSegments)
//   {
//     coeff = computeCoeff(std::get<0>(segment), std::get<1>(segment), std::get<2>(segment)) * coeff;
//   }

//   return coeff;
// }


MetricField::MetricField()
{
  exactinit();

  // std::string filename = "s01u_cube_aligned.vtk";
  // std::string filename = "s01u_cube_aligned.ovm";
  // std::string filename = "myMesh1.ovm";
  // std::string filename = "../myMesh1.ovm";
  std::string filename = "/Users/FloriGod/Development/cube_field/meshes/0_1_cube.ovm";
  // std::string filename = "../meshes/0_1_cube.ovm";
  // std::string filename = "bridge.ovm";
  // std::string filename = "s01c_cube.vtk";
  readInTetmesh(filename);
  attach_cell_vertices();

  prepareTransformation();
  attachMetric();
  attachCurl();

}

MetricField::MetricField(std::string filename)
{
  exactinit();

  // std::string filename = "s01u_cube_aligned.vtk";
  // std::string filename = "s01u_cube_aligned.ovm";
  // std::string filename = "myMesh1.ovm";
  // std::string filename = "s01c_cube.vtk";
  readInTetmesh(filename);
  std::cout << "read in tetmesh " << std::endl;
  attach_cell_vertices();
  std::cout << "attached cell vertices to property " << std::endl;
  prepareTransformation();
  std::cout << "prepare inverse matrices" << std::endl;
  attachMetric();
  std::cout << "attached metric " << std::endl;
  attachCurl();
  std::cout << "attached curl " << std::endl;
}

MetricField::MetricField(TM &tetmesh_) {
  exactinit();
  tetmesh = tetmesh_;
  attach_cell_vertices();

  prepareTransformation();
  attachMetric();
  attachCurl();
}

void MetricField::attach_cell_vertices(){
  cell_vertices = tetmesh.request_cell_property<std::vector<OVM::VertexHandle>>();

  for (auto cell : tetmesh.cells()){
    cell_vertices[cell] = tetmesh.get_cell_vertices(cell);
  }
  return;
}



// calculates the inverse matrix T^-1 for the barycentric coords
void MetricField::prepareTransformation()
{
  T_prop = tetmesh.request_cell_property<Mat4d>("T");

  for (const auto &ch : tetmesh.cells())
  {
    // const auto verts = tetmesh.get_cell_vertices(ch);
    const auto verts = cell_vertices[ch];
    Mat4d T_inverse;

    T_inverse.block<3, 1>(0, 0) = Vec3d(tetmesh.vertex(verts[0]).data());
    T_inverse.block<3, 1>(0, 1) = Vec3d(tetmesh.vertex(verts[1]).data());
    T_inverse.block<3, 1>(0, 2) = Vec3d(tetmesh.vertex(verts[2]).data());
    T_inverse.block<3, 1>(0, 3) = Vec3d(tetmesh.vertex(verts[3]).data());

    // boundary condition a+b+c+d=1
    Eigen::Matrix<double, 1, 4> ones = Eigen::Matrix<double, 1, 4>::Constant(1, 4, 1);

    T_inverse.block<1, 4>(3, 0) = ones;

    T_prop[ch] = T_inverse.inverse();
  }
}



// attaches linearly varying metric in the z-axis
void MetricField::attachMetric()
{
  Mat3d identity = Mat3d::Identity();
  // OVM::VertexPropertyT<Mat3d> metric = tetmesh.create_persistent_vertex_property<Mat3d>("metric").value();
  metric = tetmesh.request_vertex_property<Mat3d>("metric");
  for (OVM::VertexIter v_it = tetmesh.vertices_begin(); v_it != tetmesh.vertices_end(); ++v_it)
  {
    // metric[*v_it] = identity;
    //------------arc example or twisted rod, check rootMetricAnalytical implementation
    // Vec3d temp = toVec3d(tetmesh.vertex(*v_it));
    // metric[*v_it] = rootMetricAnalytical(temp);

    //--------------------------
    // double factor;
    // double k = 10.0;
    // double z = tetmesh.vertex(*v_it)[2];
    // factor = (k-1.0)*z+1.0;
    // metric[*v_it]= identity * factor;
    //--------------------------
    // double factor;
    // double k = 20.0;
    // double z = tetmesh.vertex(*v_it)[2];
    // factor = (4.5)*z+5.5;
    // metric[*v_it]= identity * factor;

    //-------- ex1 --------------- isotropic sizing, cube divided into 3 parts
    // double factor;
    // double k = 10;
    // double z = tetmesh.vertex(*v_it)[2];
    // if (z < 1.0/3.0) {
    //   factor = 1.0;
    // } else if (z > 1.0/3.0 && z < 2.0/3.0) {
    //   factor = 3.0*(k-1.0)*z-k+2.0;
    // } else {
    //   factor = k;
    // }
    // metric[*v_it]= identity * factor;

    //-------- ex2 --------------- anisotropic sizing, cube divided into 3 partss
    // double factor;
    // double k = 4;
    // double z = tetmesh.vertex(*v_it)[2];
    // if (z < 1.0/3.0) {
    //   factor = 1.0;
    // } else if (z > 1.0/3.0 && z < 2.0/3.0) {
    //   factor = 3*z*(k-1)-k+2;
    // } else {
    //   factor = k;
    // }
    // identity(0,0) = factor;
    // identity(2,2) = factor;
    metric[*v_it]= identity;

    //---------- anisotropic sizing, linearly increasing everywhere
    // double k = 4;
    // double z = tetmesh.vertex(*v_it)[2];

    // identity(0,0) = z;
    // identity(2,2) = z;
    // metric[*v_it]= identity;

    // ---------- isotropic sizing, larger at edge of the cube
    // Vec3d offset = (Vec3d() << 0.5, 0.5, 0.5).finished();
    // Vec3d point = toVec3d(tetmesh.vertex(*v_it)) - offset;

    // double k = 10.0;

    // // double factor = point.lpNorm<6>();
    // double factor = point.lpNorm<Eigen::Infinity>();

    // if (factor < 0.166666)
    // {
    //   metric[*v_it] = identity;
    // }
    // else if (factor > 0.166666 && factor < 0.166666 * 2.0)
    // {
    //   metric[*v_it] = identity * (6 * (k - 1) * factor + 2 - k);
    // }
    // else
    // {
    //   metric[*v_it] = identity * k;
    // }
  }
}

// calculates the curl in each tet and attaches it as a property
void MetricField::attachCurl()
{
  curl = tetmesh.request_cell_property<Mat3d>("curl");
  for (OVM::CellIter c_it = tetmesh.cells_begin(); c_it != tetmesh.cells_end(); ++c_it)
  {
    calculateCurl(*c_it);
  }
}

void MetricField::calculateCurl(const OVM::CellHandle &ch)
{
  // curl is constant over the tet

  const auto verts = cell_vertices[ch];
  // OVM::CellPropertyT<Mat3d> curl = tetmesh.request_cell_property<Mat3d>("curl");
  // Mat4d T = tetmesh.request_cell_property<Mat4d>("T")[ch];
  Mat4d T = T_prop[ch];
  for (int i = 0; i < 3; i++)
  {
    double value1 =
        T.col(1).dot((Vec4d() << metric[verts[0]](2, i), metric[verts[1]](2, i),
                      metric[verts[2]](2, i), metric[verts[3]](2, i))
                         .finished());
    double value2 =
        T.col(2).dot((Vec4d() << metric[verts[0]](1, i), metric[verts[1]](1, i),
                      metric[verts[2]](1, i), metric[verts[3]](1, i))
                         .finished());
    double value3 =
        T.col(2).dot((Vec4d() << metric[verts[0]](0, i), metric[verts[1]](0, i),
                      metric[verts[2]](0, i), metric[verts[3]](0, i))
                         .finished());
    double value4 =
        T.col(0).dot((Vec4d() << metric[verts[0]](2, i), metric[verts[1]](2, i),
                      metric[verts[2]](2, i), metric[verts[3]](2, i))
                         .finished());
    double value5 =
        T.col(0).dot((Vec4d() << metric[verts[0]](1, i), metric[verts[1]](1, i),
                      metric[verts[2]](1, i), metric[verts[3]](1, i))
                         .finished());
    double value6 =
        T.col(1).dot((Vec4d() << metric[verts[0]](0, i), metric[verts[1]](0, i),
                      metric[verts[2]](0, i), metric[verts[3]](0, i))
                         .finished());
    Vec3d curl_i;
    curl_i << value1 - value2, value3 - value4, value5 - value6;
    curl[ch].block<3, 1>(0, i) = curl_i; // assign to curl
  }
}

/**
 * Calculates metric at point within tet by barycentric interpolation
 *
 * @param Vec3d
 * @returns Matrix3d interpolatedMetric
 */
Mat3d MetricField::metricAtPoint(const OVM::CellHandle &ch, const Vec3d &_p)
{
  // assert(pointInTet(ch, _p));
  Mat4d T = T_prop[ch];

  Vec4d p;
  p << _p(0), _p(1), _p(2), 1;
  const Vec4d lambda = T * p;
  // const auto verts = tetmesh.get_cell_vertices(ch);
  const auto verts = cell_vertices[ch];
  Mat3d res =
      lambda.coeff(0) * metric[verts[0]] + lambda.coeff(1) * metric[verts[1]] +
      lambda.coeff(2) * metric[verts[2]] + lambda.coeff(3) * metric[verts[3]];

  return res;
}


// frobenius norm
Mat3d MetricField::recursiveDivide(const OVM::CellHandle &ch, const Vec3d &a, const Vec3d &b)
{
  Mat3d r1 = lie_exp(integrate2Point(ch, a, b));
  Vec3d midpoint = (a + b) / 2.0;
  double length = (a - b).squaredNorm();
  std::cout << "length " << length << std::endl;
  if (((r1 - lie_exp(integrate2Point(ch, a, midpoint)) * lie_exp(integrate2Point(ch, midpoint, b))).squaredNorm()) /
          length <
      10e-6)
  {
    return r1;
  }
  else
  {
    return recursiveDivide(ch, a, midpoint) *
           recursiveDivide(ch, midpoint, b);
  }
}

/**
 * calculates matrix exponential with rodrigues formula for antisymmetric
 * matrices
 *
 * @param Vec3d u rotation axis
 * @returns rotation around u by angle norm(u)
 */
Mat3d lie_exp(Vec3d u)
{
  double angle = u.norm();
  u.normalize();
  Mat3d mat = Ai_x(u);
  Mat3d result = Mat3d().Identity() + sin(angle) * mat + (1 - cos(angle)) * mat * mat;
  return result;
}


/**
 * The precomputed inverse of A, to solve the linear system
 * 
 * @param Mat9d A
 * @returns Mat9d A^{-1}
*/
Mat9d MetricField::scaledInverse(const Mat3d &A)
{
  double factor = A.determinant() * -2.;
  Mat9d inv;
  double a_11 = A(0, 0);
  double a_12 = A(0, 1);
  double a_13 = A(0, 2);
  double a_22 = A(1, 1);
  double a_23 = A(1, 2);
  double a_33 = A(2, 2);

  inv << a_11 * a_11, a_11 * a_12, a_11 * a_13, a_11 * a_12, 2 * a_12 * a_12 - a_11 * a_22, 2 * a_12 * a_13 - a_11 * a_23, a_11 * a_13, 2 * a_12 * a_13 - a_11 * a_23, 2 * a_13 * a_13 - a_11 * a_33,
      a_11 * a_12, a_12 * a_12, a_12 * a_13, -a_12 * a_12 + 2 * a_11 * a_22, a_12 * a_22, 2 * a_13 * a_22 - a_12 * a_23, -a_12 * a_13 + 2 * a_11 * a_23, a_12 * a_23, 2 * a_13 * a_23 - a_12 * a_33,
      a_11 * a_13, a_12 * a_13, a_13 * a_13, -a_12 * a_13 + 2 * a_11 * a_23, -a_13 * a_22 + 2 * a_12 * a_23, a_13 * a_23, -a_13 * a_13 + 2 * a_11 * a_33, -a_13 * a_23 + 2 * a_12 * a_33, a_13 * a_33,
      a_11 * a_12, -a_12 * a_12 + 2 * a_11 * a_22, -a_12 * a_13 + 2 * a_11 * a_23, a_12 * a_12, a_12 * a_22, a_12 * a_23, a_12 * a_13, 2 * a_13 * a_22 - a_12 * a_23, 2 * a_13 * a_23 - a_12 * a_33,
      2 * a_12 * a_12 - a_11 * a_22, a_12 * a_22, -a_13 * a_22 + 2 * a_12 * a_23, a_12 * a_22, a_22 * a_22, a_22 * a_23, -a_13 * a_22 + 2 * a_12 * a_23, a_22 * a_23, 2 * a_23 * a_23 - a_22 * a_33,
      2 * a_12 * a_13 - a_11 * a_23, 2 * a_13 * a_22 - a_12 * a_23, a_13 * a_23, a_12 * a_23, a_22 * a_23, a_23 * a_23, -a_13 * a_23 + 2 * a_12 * a_33, -a_23 * a_23 + 2 * a_22 * a_33, a_23 * a_33,
      a_11 * a_13, -a_12 * a_13 + 2 * a_11 * a_23, -a_13 * a_13 + 2 * a_11 * a_33, a_12 * a_13, -a_13 * a_22 + 2 * a_12 * a_23, -a_13 * a_23 + 2 * a_12 * a_33, a_13 * a_13, a_13 * a_23, a_13 * a_33,
      2 * a_12 * a_13 - a_11 * a_23, a_12 * a_23, -a_13 * a_23 + 2 * a_12 * a_33, 2 * a_13 * a_22 - a_12 * a_23, a_22 * a_23, -a_23 * a_23 + 2 * a_22 * a_33, a_13 * a_23, a_23 * a_23, a_23 * a_33,
      2 * a_13 * a_13 - a_11 * a_33, 2 * a_13 * a_23 - a_12 * a_33, a_13 * a_33, 2 * a_13 * a_23 - a_12 * a_33, 2 * a_23 * a_23 - a_22 * a_33, a_23 * a_33, a_13 * a_33, a_23 * a_33, a_33 * a_33;

  return (1./factor) * inv;
}

/**
 * Evaluates the coefficients of the one-forms W through a 9x9 linear system
 *
 * @param Vec3d point
 * @returns Matrix3d W
 */
Mat3d MetricField::eval_W(const OVM::CellHandle &ch, const Vec3d &_p)
{
  // OVM::CellPropertyT<Mat3d> curl = tetmesh.request_cell_property<Mat3d>("curl");
  
  return unstack(scaledInverse(metricAtPoint(ch, _p)) * stack(curl[ch]));
  // return unstack(A_x(metricAtPoint(ch, _p)).partialPivLu().solve(stack(curl[ch])));;
  // return unstack(A_x(metricAtPoint(ch, _p)).fullPivLu().solve(stack(curl[ch])));
}


/**
 * @returns curl of certain tet
*/
Mat3d MetricField::getCurl(const OVM::CellHandle &ch){
  
  // OVM::CellPropertyT<Mat3d> curl = tetmesh.request_cell_property<Mat3d>("curl");
  return curl[ch];
}

/**
 * integrates with trapezoidal rule between point a and b within a tet
 * not weighted since edge lies within tet
 *
 * @param Vec3d a
 * @param Vec3d b
 * @returns Vec3d integrated vector
 */
Vec3d MetricField::integrate2Point(const OVM::CellHandle &ch, const Vec3d &a, const Vec3d &b)
{
  // std::cout << "eval W linear \n" << eval_W(ch, a) << std::endl;
  Vec3d res = 0.5 * (eval_W(ch, a) + eval_W(ch, b)).transpose() * (b - a);
  return res;
}


/**
 * computes the rotation coefficient between points p and q where p/q lie in the same cell
 * 
 * @param Vec3d q
 * @param Vec3d p
 * @returns Mat3d rotation
*/
Mat3d MetricField::computeCoeff(const OVM::CellHandle &ch, const Vec3d &q, const Vec3d &p)
{
  return lie_exp(integrate2Point(ch, q, p)); // go only to boundary of cell
  //return recursiveDivide(ch, q, p); // do normal recursive division
}

/**
 * arranges 3x3 matrix A to 9x9 cross matrix
 *
 * @param Mat3d A
 * @returns MatrixXd 9x9 cross matrix
 */
Mat9d A_x(const Mat3d &A)
{
  Mat9d res;

  Mat3d A1_x = Ai_x(A.col(0));
  Mat3d A2_x = Ai_x(A.col(1));
  Mat3d A3_x = Ai_x(A.col(2));
  Mat3d zero = Mat3d::Zero();

  res.block<3, 3>(0, 0) = zero;
  res.block<3, 3>(0, 3) = -A3_x;
  res.block<3, 3>(0, 6) = A2_x;
  res.block<3, 3>(3, 0) = A3_x;
  res.block<3, 3>(3, 3) = zero;
  res.block<3, 3>(3, 6) = -A1_x;
  res.block<3, 3>(6, 0) = -A2_x;
  res.block<3, 3>(6, 3) = A1_x;
  res.block<3, 3>(6, 6) = zero;

  return res;
}

/**
 * arranges column to 3x3 cross matrix
 *
 * @param Vec3d column to be arranged
 * @returns Matrix3d cross matrix
 */
Mat3d Ai_x(const Vec3d &A_i)
{
  Mat3d res;

  res << 0, -A_i(2), A_i(1), A_i(2), 0, -A_i(0), -A_i(1), A_i(0), 0;

  return res;
}

/**
 * stacks columns of 3x3 Matrix to 9x1 Vector
 *
 * @param Mat3d mat
 * @returns VectorXd stacked vector
 */
Vec9d stack(const Mat3d &_v) { return _v.reshaped(); }

/**
 * unstacks columns of 9x1 Vector to 3x3 Matrix
 *
 * @param Vec9d vector
 * @returns Matrix3d unstacked vector
 */
Mat3d unstack(const Vec9d &_m)
{
  Mat3d res = _m.reshaped(3, 3);
  return res;
}

//-------------------------------------------TET_FINDER_CODE------------------------------------

/**
 * TOBIAS VERY SUPER SMART TRACING CODE
 */

OVM::HalfFaceHandle MetricField::findOppositeNextHalfFace(const OVM::CellHandle& ch, const OVM::HalfFaceHandle& hfh, const Vec3d &q, const Vec3d &p)
    {
    assert(!hfh.is_valid() || tetmesh.incident_cell(hfh) == ch);

    for (const auto& checkHfh : tetmesh.cell(ch).halffaces()) {
        if (checkHfh == hfh) continue;

        const auto& vertices = tetmesh.get_halfface_vertices(checkHfh);
        auto u = vertices[0]; // params[0]
        auto v = vertices[1];// params[1]
        auto w = vertices[2];// params[2]


        // Primary direction must cut through the plane given by the face
        // if (ori3d(params[0], params[1], params[2], u) != ORI_ABOVE || ori3d(params[0], params[1], params[2], u + d1) != ORI_BELOW)
        if (orient3dHelper(u, v, w, q) != 1 || orient3dHelper(u, v, w, p) != -1)
            continue;

        // Get the orientations of the 3 tets formed around u, u+d1
        std::array<int, 3> oris;
        oris[0] = orient3dHelper(u, v, q, p);
        oris[1] = orient3dHelper(v, w, q, p);
        oris[2] = orient3dHelper(w, u, q, p);

        // Check if primary direction cuts through interior/center of triangle
        if (oris[0] == oris[1] && oris[0] == oris[2]) {
            assert(oris[0] == -1);
            return checkHfh; // case 1
        }

        // Check if primary direction does not cut through triangle at all
        if (std::any_of(oris.begin(), oris.end(), [](int ori){return ori == -1;}) && std::any_of(oris.begin(), oris.end(), [](int ori){return ori == 1;})) {
            continue;
        }

        // Otherwise the primary direction cuts through the triangles boundary (edge or vertex). First get all edge intersection (amount must be 1 or 2)
        std::vector<int> edgeIntersections;
        for (unsigned char i = 0; i < oris.size(); ++i)
        { 
          if (oris[i] == 0) edgeIntersections.push_back(i);
        }

        if (edgeIntersections.size() == 1) {  // Pointing through edge - case 2
            return checkHfh;

        } else { // Pointing through vertex - case 3
            assert(edgeIntersections.size() == 2);

            return checkHfh;
        }
    }


    }


/**
 * Given a start_tet where q is containted, returns tet containing p, starting from q, without calculating all intersection points.
 */
OVM::CellHandle MetricField::startCellFinder(OVM::CellHandle start_tet, const Vec3d &q, const Vec3d &p) {
  auto result = tetFinderTobias(start_tet, q, p, false);
  auto last_segment = result.back();

  auto contained_cell = std::get<0>(last_segment);
  assert(pointInTet(contained_cell, p));
  return contained_cell;
}

/**
 * assumes q is contained within start_tet
 */
std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> MetricField::tetFinderTobias(OVM::CellHandle start_tet, const Vec3d &q, const Vec3d &p, bool compute_intersections) {
  assert(pointInTet(start_tet, q));

  std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> ints;

  if (pointInTet(start_tet, p)){ // if p is already contained, we are done here
    ints.push_back(std::make_tuple(start_tet, q, p));
    return ints;
  }

  
  // get correct start_tet, initial tet may not be valid because the ray from qp does not go through a correct halfface
  // std::cout << "curr tet " << start_tet << std::endl;
  start_tet = findCorrectStartTet(start_tet, q, p);
  OVM::CellHandle curr_tet = start_tet;
  OVM::HalfFaceHandle curr_face = OVM::HalfFaceHandle(-1);

  Vec3d prev = q;


  while(!pointInTet(curr_tet, p)){
    curr_face = findOppositeNextHalfFace(curr_tet, curr_face, q, p).opposite_handle();
    curr_tet = tetmesh.incident_cell(curr_face);

    if (compute_intersections){
      auto verts =  tetmesh.get_halfface_vertices(curr_face);

      Vec3d curr = intersection(q, p, verts[0], verts[1], verts[2]);
      ints.push_back(std::make_tuple(curr_tet, prev, curr));

      prev = curr;
        // std::cout << "curr tet " << curr_tet << std::endl;

    } else {
      ints.push_back(std::make_tuple(curr_tet, Vec3d(0,0,0), Vec3d(0,0,0)));
        // std::cout << "starting curr tet " << curr_tet << std::endl;

    }
    

  }

  assert(ints.size() > 0);

  if (compute_intersections){
    // also make sure to get element from q to first intersection point
    auto first_segment = ints.at(0);

    ints.push_back(std::make_tuple(start_tet, q, std::get<1>(first_segment)));
    std::rotate(ints.rbegin(), ints.rbegin() + 1, ints.rend());

    // append from last intersection point to p
    ints.push_back(std::make_tuple(curr_tet, prev, p));
  }
  

  return ints;
}

OVM::CellHandle MetricField::findCorrectStartTet(OVM::CellHandle start_tet, const Vec3d &q, const Vec3d &p){
  MeshElement start_element = MeshElement();

  bool valid_start = pointInTetWithElement(start_tet, toOVMVec3d(q), start_element);
  assert(valid_start);

  if (start_element.is_cell()){
    // all good

  } else if (start_element.is_face()){
    auto verts = tetmesh.get_halfface_vertices(start_element.hfh());

    // take cell that has remaining vertex on same side as p
    if (orient3dHelper(verts[0], verts[1], verts[2], p) == -1){ // not on same side
      start_tet = tetmesh.incident_cell(start_element.hfh().opposite_handle());
    }
  } else if (start_element.is_edge()){
    for (auto ec_it = tetmesh.ec_iter(start_element.eh()); ec_it.is_valid(); ++ec_it){
      const auto& hfhs = getIncidentHalfFaces(*ec_it, start_element.eh());
      assert(hfhs.size()==2);
      if (orient3dHelper(hfhs[0], p) == -1) continue;
      if (orient3dHelper(hfhs[1], p) == -1) continue;
      start_tet = *ec_it;
      break;
    }
  } else if (start_element.is_vertex()){
    for (auto vc_it = tetmesh.vc_iter(start_element.vh()); vc_it.is_valid(); ++vc_it){
      const auto& hfhs = getIncidentHalfFaces(*vc_it, start_element.vh());
      assert(hfhs.size()==3);
      if (orient3dHelper(hfhs[0], p) == -1) continue;
      if (orient3dHelper(hfhs[1], p) == -1) continue;
      if (orient3dHelper(hfhs[2], p) == -1) continue;
      start_tet = *vc_it;
      break;
    }
  }
  return start_tet;
}

bool MetricField::pointInTetWithElement(OVM::CellHandle ch, const OVM::Vec3d param, MeshElement& elem)
{
    // Get the 4 vertices of the tet and their parameters in the chart of the tet
    const auto& vhs = cell_vertices[ch];



    std::vector<OVM::Vec3d> params = {tetmesh.vertex(vhs[0]), tetmesh.vertex(vhs[1]), tetmesh.vertex(vhs[2]), tetmesh.vertex(vhs[3])};
    for (char i = 0; i <= 3; ++i) {if (params[i] == param) {elem.set(vhs[i]); return true;}}

    auto oris = std::array<int, 4>();
    oris[0] = orient3dHelper(params[0], params[1], params[2], param); if (oris[0] == -1) return false;
    oris[1] = orient3dHelper(params[0], params[2], params[3], param); if (oris[1] == -1) return false;
    oris[2] = orient3dHelper(params[0], params[3], params[1], param); if (oris[2] == -1) return false;
    oris[3] = orient3dHelper(params[1], params[3], params[2], param); if (oris[3] == -1) return false;

    auto zeros = std::vector<char>();
    zeros.reserve(2);
    for (char i = 0; i <= 3; ++i) {if (oris[i] == 0) {zeros.push_back(i);}}

    if (zeros.size() == 2) {
        if (zeros[0] == 0 && zeros[1] == 1) {elem.set(tetmesh.find_halfedge(vhs[0], vhs[2]));}
        else if (zeros[0] == 0 && zeros[1] == 2) {elem.set(tetmesh.find_halfedge(vhs[0], vhs[1]));}
        else if (zeros[0] == 0 && zeros[1] == 3) {elem.set(tetmesh.find_halfedge(vhs[1], vhs[2]));}
        else if (zeros[0] == 1 && zeros[1] == 2) {elem.set(tetmesh.find_halfedge(vhs[0], vhs[3]));}
        else if (zeros[0] == 1 && zeros[1] == 3) {elem.set(tetmesh.find_halfedge(vhs[2], vhs[3]));}
        else if (zeros[0] == 2 && zeros[1] == 3) {elem.set(tetmesh.find_halfedge(vhs[1], vhs[3]));}
        return true;
    }
    if (zeros.size() == 1) {
        if (zeros[0] == 0) {elem.set(tetmesh.find_halfface_in_cell({vhs[0], vhs[1], vhs[2]}, ch));}
        else if (zeros[0] == 1) {elem.set(tetmesh.find_halfface_in_cell({vhs[0], vhs[2], vhs[3]}, ch));}
        else if (zeros[0] == 2) {elem.set(tetmesh.find_halfface_in_cell({vhs[0], vhs[3], vhs[1]}, ch));}
        else if (zeros[0] == 3) {elem.set(tetmesh.find_halfface_in_cell({vhs[1], vhs[3], vhs[2]}, ch));}
        return true;
    }
    if (zeros.size() == 0) {
        elem.set(ch);
        return true;
    }
    return false;
}


// returns intersection point with triangle uvw, Müller-Trombone Algorithm
Vec3d MetricField::intersection(const Vec3d &_q, const Vec3d &_p,
                                OVM::VertexHandle &_u, OVM::VertexHandle &_v,
                                OVM::VertexHandle &_w)
{
  OVM::Vec3d v0, v1, v2;
  v0 = tetmesh.vertex(_u);
  v1 = tetmesh.vertex(_v);
  v2 = tetmesh.vertex(_w);
  OVM::Vec3d edge1, edge2, h, s, q, rayOrigin, rayVector;
  rayOrigin = OVM::Vec3d(_q.data());
  rayVector = toOVMVec3d(_p - _q);
  double a, f, u, v;
  edge1 = v1 - v0;
  edge2 = v2 - v0;
  h = rayVector.cross(edge2);
  a = edge1.dot(h);

  double orientQ = orient3d(_q.data(), v0.data(), v1.data(), v2.data());
  double orientP = orient3d(_p.data(), v0.data(), v1.data(), v2.data());

  // coplanar ray triangle, 2D case
  if (orientQ == 0.0 && orientP == 0.0)
  {
    std::cout << "coplanar" << std::endl;
    return Vec3d(0.0, 0.0, 0.0);
  }

  f = 1.0 / a;
  s = rayOrigin - v0;
  u = f * s.dot(h);


  q = s.cross(edge1);
  v = f * rayVector.dot(q);  

  // Compute t to find out where intersection point is
  double t = f * edge2.dot(q);
  return toVec3d(rayOrigin + rayVector * t);
}

OVM::Vec3d toOVMVec3d(const Vec3d &_q) { return OVM::Vec3d(_q.data()); }
// OVM::Vec3d toOVMVec3d(Vec3d &_q) { return OVM::Vec3d(_q.data()); }

Vec3d toVec3d(const OVM::Vec3d &_q) { return Vec3d(_q.data()); }

// check if two doubles have same sign, if one value is zero, it returns true
bool sameSign(double a, double b)
{
  return a == 0.0 || b == 0.0 || std::signbit(a) == std::signbit(b);
}

// does not handle the case where q is exactly on a vertex
bool MetricField::pointInTet(const OVM::CellHandle &ch, const Vec3d &q)
{
  // std::vector<OVM::VertexHandle> vertices = tetmesh.get_cell_vertices(ch);
  std::vector<OVM::VertexHandle> vertices = cell_vertices[ch];
  const double *v1 = tetmesh.vertex(vertices[0]).data();
  const double *v2 = tetmesh.vertex(vertices[1]).data();
  const double *v3 = tetmesh.vertex(vertices[2]).data();
  const double *v4 = tetmesh.vertex(vertices[3]).data();

  return sameSign(orient3d(v1, v2, v3, q.data()), orient3d(v1, v2, v3, v4)) &&
         sameSign(orient3d(v2, v3, v4, q.data()), orient3d(v2, v3, v4, v1)) &&
         sameSign(orient3d(v3, v4, v1, q.data()), orient3d(v3, v4, v1, v2)) &&
         sameSign(orient3d(v4, v1, v2, q.data()), orient3d(v4, v1, v2, v3));
}

double orient3dPaper(const double *a, const double *b, const double *c,
                     const double *d)
{
  const double result = orient3d(b, c, d, a);
  return (result > 0.0) - (result < 0.0);
}

double MetricField::orient3dHelper(const OVM::VertexHandle a,
                                   const OVM::VertexHandle b,
                                   const OVM::VertexHandle c, const Vec3d d)
{
  return orient3dPaper(tetmesh.vertex(a).data(), tetmesh.vertex(b).data(),
                       tetmesh.vertex(c).data(), d.data());
}

double MetricField::orient3dHelper(const OVM::HalfFaceHandle hfh,
                                   const Vec3d d)
{
  auto verts = tetmesh.get_halfface_vertices(hfh);


  return orient3dPaper(tetmesh.vertex(verts[0]).data(), tetmesh.vertex(verts[1]).data(),
                       tetmesh.vertex(verts[2]).data(), d.data());
}

double MetricField::orient3dHelper(const OVM::VertexHandle a,
                                   const OVM::VertexHandle b, const Vec3d c,
                                   const Vec3d d)
{
  return orient3dPaper(tetmesh.vertex(a).data(), tetmesh.vertex(b).data(),
                       c.data(), d.data());
}

double MetricField::orient3dHelper(const OVM::VertexHandle a,
                                   const OVM::VertexHandle b,
                                   const OVM::VertexHandle c,
                                   const OVM::VertexHandle d)
{
  return orient3dPaper(tetmesh.vertex(a).data(), tetmesh.vertex(b).data(),
                       tetmesh.vertex(c).data(), tetmesh.vertex(d).data());
}

double MetricField::orient3dHelper(const OVM::Vec3d a, const OVM::Vec3d b,
                        const OVM::Vec3d c, const OVM::Vec3d d) {
  return orient3dPaper(a.data(), b.data(),
                       c.data(), d.data());
                        }

// returns true if v contains x
bool contains(const std::vector<OVM::CellHandle> &v, OVM::CellHandle x)
{
  return std::find(v.begin(), v.end(), x) != v.end();
}

// removes tets from tetmesh that are not contained in tets
void markTetsAsIntersected(TM &tetmesh, std::vector<OVM::CellHandle> tets) {
    OVM::CellPropertyT<bool> intersected = tetmesh.create_persistent_cell_property<bool>("intersected").value();
    for (const OVM::CellHandle &ch : tetmesh.cells()) {
        intersected[ch] = false;
    }
    for (const OVM::CellHandle &ch : tets){
        intersected[ch] = true;
        // std::cout << "intersected " << ch << std::endl;
    }
    // to visualize in OpenFlipper
    for (const OVM::CellHandle &ch : tetmesh.cells()) {
        if (!intersected[ch]) {
            tetmesh.delete_cell(ch);
        }
    }
    
    saveToFile(tetmesh, "intersected.ovm");
    
}

void showIntersectionPoints(TM &tetmesh, std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> ints) {
    TM points;
    for (const auto &el : ints) {
        points.add_vertex(toOVMVec3d(std::get<1>(el)));
        points.add_vertex(toOVMVec3d(std::get<2>(el)));
    }
    saveToFile(points, "intersections.ovm");
}


std::vector<OVM::HalfFaceHandle> MetricField::getIncidentHalfFaces(const OVM::CellHandle& ch, const OVM::EdgeHandle& eh)
{
    OVM::OpenVolumeMeshCell c = tetmesh.cell(ch);
    OVM::OpenVolumeMeshEdge e = tetmesh.edge(eh);
    std::vector<OVM::HalfFaceHandle> res;
    res.reserve(2);
    for (const auto& hfh : c.halffaces())
        for (const auto& heh : tetmesh.halfface(hfh).halfedges())
            if (heh.edge_handle() == eh) {
                res.push_back(hfh);
                if (res.size()==2) return res;
            }
    return res;
}

std::vector<OVM::HalfFaceHandle> MetricField::getIncidentHalfFaces(const OVM::CellHandle& ch, const OVM::VertexHandle& vh)
{
    OVM::OpenVolumeMeshCell c = tetmesh.cell(ch);
    std::vector<OVM::HalfFaceHandle> res;
    res.reserve(3);
    for (const auto& hfh : c.halffaces())
        for (const auto& heh : tetmesh.halfface(hfh).halfedges())
            if (tetmesh.from_vertex_handle(heh) == vh) {
                res.push_back(hfh);
                if (res.size()==3) return res;
            }
    return res;
}