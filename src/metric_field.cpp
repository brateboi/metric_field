#include <Eigen/Dense>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

#include "metric_field.hpp"

#include <cmath>
#include <limits>
#include <other_fields.hpp>

using Quaternion = Eigen::Quaterniond;
using Vec3d = Eigen::Matrix<double, 3, 1>;
using Vec4d = Eigen::Matrix<double, 4, 1>;
using Vec9d = Eigen::Matrix<double, 9, 1>;
using Mat3d = Eigen::Matrix<double, 3, 3>;
using Mat4d = Eigen::Matrix<double, 4, 4>;
using Mat9d = Eigen::Matrix<double, 9, 9>;

using TM = OVM::TetrahedralGeometryKernel<OVM::Vec3d,
                                          OVM::TetrahedralMeshTopologyKernel>;



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


// computes the rotation between two points in the metric field, main function
Mat3d MetricField::computeCoeff(const Vec3d &q, const Vec3d &p)
{

  std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> lineSegments = tetFinder(q, p);
  Mat3d coeff = Mat3d::Identity();
  for (const auto &segment : lineSegments)
  {
    coeff = computeCoeff(std::get<0>(segment), std::get<1>(segment), std::get<2>(segment)) * coeff;
  }

  return coeff;
}

// computes the rotation between two points in the metric field, main function, but with info in which CellHandle p/q are
Mat3d MetricField::computeCoeff(const Vec3d &q, const Vec3d &p, OVM::CellHandle &cq, OVM::CellHandle &cp)
{
  // fill hash with tet info
  std::vector<OVM::CellHandle> tets;
  tets.push_back(cq);
  hash.try_emplace(q, tets);

  tets.clear();
  tets.push_back(cp);
  hash.try_emplace(p, tets);


  std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> lineSegments = tetFinder(q, p);
  Mat3d coeff = Mat3d::Identity();
  for (const auto &segment : lineSegments)
  {
    coeff = computeCoeff(std::get<0>(segment), std::get<1>(segment), std::get<2>(segment)) * coeff;
  }

  return coeff;
}


// computes the rotation between two points in the metric field, main function, pushes traversed cells to vector
Mat3d MetricField::computeCoeff(const Vec3d &q, const Vec3d &p, std::vector<OVM::CellHandle> &cells)
{


  std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> lineSegments = tetFinder(q, p);
  Mat3d coeff = Mat3d::Identity();
  for (const auto &segment : lineSegments)
  {
    // for visualizing
    // std::cout << "Cell handle " << std::get<0>(segment) << std::endl;
    cells.push_back(std::get<0>(segment));


    coeff = computeCoeff(std::get<0>(segment), std::get<1>(segment), std::get<2>(segment)) * coeff;
  }


  // for visualizing---------
  // markTetsAsIntersected(this->tetmesh, cells);
  // showIntersectionPoints(this->tetmesh, lineSegments);
  // ------------------------

  return coeff;
}

MetricField::MetricField()
{
  exactinit();

  // std::string filename = "s01u_cube_aligned.vtk";
  // std::string filename = "s01u_cube_aligned.ovm";
  // std::string filename = "myMesh1.ovm";
  std::string filename = "../myMesh1.ovm";
  // std::string filename = "bridge.ovm";
  // std::string filename = "s01c_cube.vtk";
  readInTetmesh(filename);
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
  prepareTransformation();
  attachMetric();
  attachCurl();
}

MetricField::MetricField(TM &tetmesh_) {
  exactinit();
  tetmesh = tetmesh_;

  prepareTransformation();
  attachMetric();
  attachCurl();
}



// calculates the inverse matrix T^-1 for the barycentric coords
void MetricField::prepareTransformation()
{
  OVM::CellPropertyT<Mat4d> T = tetmesh.create_persistent_cell_property<Mat4d>("T").value();
  for (const auto &ch : tetmesh.cells())
  {
    const auto verts = tetmesh.get_cell_vertices(ch);
    Mat4d T_inverse;

    T_inverse.block<3, 1>(0, 0) = Vec3d(tetmesh.vertex(verts[0]).data());
    T_inverse.block<3, 1>(0, 1) = Vec3d(tetmesh.vertex(verts[1]).data());
    T_inverse.block<3, 1>(0, 2) = Vec3d(tetmesh.vertex(verts[2]).data());
    T_inverse.block<3, 1>(0, 3) = Vec3d(tetmesh.vertex(verts[3]).data());

    // boundary condition a+b+c+d=1
    Eigen::Matrix<double, 1, 4> ones = Eigen::Matrix<double, 1, 4>::Constant(1, 4, 1);

    T_inverse.block<1, 4>(3, 0) = ones;

    T[ch] = T_inverse.inverse();
  }
}



// attaches linearly varying metric in the z-axis
void MetricField::attachMetric()
{

  Mat3d identity = Mat3d::Identity();
  OVM::VertexPropertyT<Mat3d> metric = tetmesh.create_persistent_vertex_property<Mat3d>("metric").value();
  for (OVM::VertexIter v_it = tetmesh.vertices_begin(); v_it != tetmesh.vertices_end(); ++v_it)
  {
    //------------arc example or twisted rod, check rootMetricAnalytical implementation
    // Vec3d temp = toVec3d(tetmesh.vertex(*v_it));
    // metric[*v_it] = rootMetricAnalytical(temp);

    //--------------------------
    // double factor;
    // double k = 10.0;
    // double z = tetmesh.vertex(*v_it)[2];
    // factor = (k-1.0)*z+1.0;
    // metric[*v_it]= identity * factor;

    //-------- ex1 --------------- isotropic sizing, cube divided into 3 parts
    double factor;
    double k = 10;
    double z = tetmesh.vertex(*v_it)[2];
    if (z < 1.0/3.0) {
      factor = 1.0;
    } else if (z > 1.0/3.0 && z < 2.0/3.0) {
      factor = 3.0*(k-1.0)*z-k+2.0;
    } else {
      factor = k;
    }
    metric[*v_it]= identity * factor;

    //-------- ex2 --------------- anisotropic sizing, cube divided into 3 partss
    // double factor;
    // double k = 4;
    // double z = tetmesh.vertex(*v_it)[2];
    // if (z < 1.0/3.0) {
    //   factor = 1.0;
    // } else if (z > 1.0/3.0 && z < 2.0/3.0) {
    //   factor = 3.0*(k-1.0)*z-k+2.0;
    // } else {
    //   factor = k;
    // }
    // identity(0,0) = factor;
    // identity(2,2) = factor;
    // metric[*v_it]= identity;

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
  OVM::CellPropertyT<Mat3d> curl = tetmesh.create_persistent_cell_property<Mat3d>("curl").value();
  for (OVM::CellIter c_it = tetmesh.cells_begin(); c_it != tetmesh.cells_end(); ++c_it)
  {
    calculateCurl(*c_it);
  }
}

void MetricField::calculateCurl(const OVM::CellHandle &ch)
{
  // curl is constant over the tet
  OVM::VertexPropertyT<Mat3d> metric = tetmesh.request_vertex_property<Mat3d>("metric");
  const auto verts = tetmesh.get_cell_vertices(ch);
  OVM::CellPropertyT<Mat3d> curl = tetmesh.request_cell_property<Mat3d>("curl");
  Mat4d T = tetmesh.request_cell_property<Mat4d>("T")[ch];
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
  Mat4d T = tetmesh.request_cell_property<Mat4d>("T")[ch];
  OVM::VertexPropertyT<Mat3d> metric = tetmesh.request_vertex_property<Mat3d>("metric");

  Vec4d p;
  p << _p(0), _p(1), _p(2), 1;
  Vec4d lambda = T * p;
  const auto verts = tetmesh.get_cell_vertices(ch);
  Mat3d res =
      lambda.coeff(0) * metric[verts[0]] + lambda.coeff(1) * metric[verts[1]] +
      lambda.coeff(2) * metric[verts[2]] + lambda.coeff(3) * metric[verts[3]];

  return res;
}

// frobenius norm
Mat3d MetricField::recursiveDivide(const OVM::CellHandle &ch, const Vec3d &a, const Vec3d &b, int depth)
{
  

  Mat3d r1 = lie_exp(integrate2Point(ch, a, b));
  Vec3d midpoint = (a + b) / 2.0;
  double length = (a - b).squaredNorm();
  if (((r1 - lie_exp(integrate2Point(ch, a, midpoint)) * lie_exp(integrate2Point(ch, midpoint, b))).squaredNorm()) /
          length <
      10e-6)
  {
    return r1;
  }
  else
  {
    return recursiveDivide(ch, a, midpoint, ++depth) *
           recursiveDivide(ch, midpoint, b, ++depth);
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
  OVM::CellPropertyT<Mat3d> curl = tetmesh.request_cell_property<Mat3d>("curl");
  
  return unstack(scaledInverse(metricAtPoint(ch, _p)) * stack(curl[ch]));
  // return unstack(A_x(metricAtPoint(ch, _p)).partialPivLu().solve(stack(curl[ch])));;
  // return unstack(A_x(metricAtPoint(ch, _p)).fullPivLu().solve(stack(curl[ch])));
}

Mat3d MetricField::getCurl(const OVM::CellHandle &ch){
  
  OVM::CellPropertyT<Mat3d> curl = tetmesh.request_cell_property<Mat3d>("curl");
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

Mat3d MetricField::computeCoeff(const OVM::CellHandle &ch, const Vec3d &q, const Vec3d &p)
{
  // return lie_exp(integrate2Point(ch, q, p)); // go only to boundary of cell
  return recursiveDivide(ch, q, p, 0);
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
 * @brief walks through the mesh starting from q going to p
 * marks a path of tets and their intersection points with the ray inbetween
 * them
 *
 * @param q
 * @param p
 * @return std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>>
 */
std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>>
MetricField::tetFinder(const Vec3d &q, const Vec3d &p)
{
  // we work with OVM::Vec3d here, copy them
  // const OVM::Vec3d q = OVM::Vec3d(_q.data());
  // const OVM::Vec3d p = OVM::Vec3d(_p.data());

  auto starts = locateTets(q);
  OVM::CellHandle t;
  OVM::VertexHandle u, v, w, s;
  OVM::HalfFaceHandle triangle;
  std::vector<OVM::VertexHandle> verts;
  std::vector<std::tuple<OVM::CellHandle, Vec3d, Vec3d>> ints;
  bool success = false;
  for (OVM::CellHandle _t : starts)
  {
    // if p is in same tet as q, then we are done here
    if (pointInTet(_t, p))
    {
      ints.push_back(std::make_tuple(_t, q, p));
      return ints;
    }
    if (success = initialization(_t, q, p, u, v, w))
    {
      t = _t;
      break;
    };
  }

  // where p is actually located, to know if we are finished
  std::vector<OVM::CellHandle> t_others = locateTets(p);

  // start walk
  Vec3d prev = q;
  if (!success)
  {
    std::cout << "failed to initialize correctly" << std::endl;
    std::cout << "q" << q.transpose() << std::endl;
    std::cout << "p" << p.transpose() << std::endl;
    for (const auto _t : starts)
    {
      std::cout << "\t tried tet " << _t << std::endl;
      const auto temp = tetmesh.get_cell_vertices(_t);
      for (const auto &vh : temp)
      {
        std::cout << "\t \t vertices " << tetmesh.vertex(vh) << std::endl;
      }
    }
  }
  Vec3d curr = intersection(q, p, u, v, w);
  ints.push_back(std::make_tuple(t, prev, curr));

  // qp intersects triangle uvw
  // wvqp, vuqp, uwqp are positively oriented
  while (orient3dHelper(u, w, v, p) > 0 || !contains(t_others, t))
  {
    // degenerate cases can lead to invalid states, need reinit
    if (/*orient3dHelper(u, w, v, p) < 0 ||*/ !validState(q, p, u, v, w))
    {
      // std::cout << "p " << p.transpose() << std::endl; 
      // std::cout << "q " << q.transpose() << std::endl;
      std::cout << "reinitialization" << std::endl;
      
      // after a reinitialization, the next triangle traversed must not
      // necessarily be a triangle that gets intersected
      initialization(t, q, p, u, v, w);
    }
    // std::cout << "validstate before reassignment: -----------" << validState(q,p,u,v,w) << std::endl;
    // std::cout << "u " << tetmesh.vertex(u) << std::endl; 
    // std::cout << "v " << tetmesh.vertex(v) << std::endl; 
    // std::cout << "w " << tetmesh.vertex(w) << std::endl; 

    verts.clear();
    verts.push_back(u);
    verts.push_back(v);
    verts.push_back(w);
    triangle = find_halfface_in_cell(verts, t);

    // jump to neighbor
    triangle = tetmesh.opposite_halfface_handle(triangle);
    t = tetmesh.incident_cell(triangle);


    // reassign s
    const auto candidates = tetmesh.get_cell_vertices(t);
    for (const OVM::VertexHandle &vh : candidates)
    {
      if (vh != u && vh != v && vh != w)
      {
        s = vh;
        break;
      }
    }
    prev = curr;
    if (orient3dHelper(u, s, q, p) > 0)
    { // qp does not intersect triangle usw
      if (orient3dHelper(v, s, q, p) > 0)
      { // qp intersects triangle vsw
        u = s;
      }
      else
      { // qp intersects triangle usv
        w = s;
      }
    }
    else
    { // qp doesnt not intersect triangle usv
      if (orient3dHelper(w, s, q, p) > 0)
      { // qp intersects triangle usw
        v = s;
      }
      else
      { // qp intersects triangle vsw
        u = s;
      }
    }

    // if we have to make reinitialization, triangle traversed isn't necessarily
    // a triangle that gets intersected check before each intersection
    if (robustRayTriangle(q, p, tetmesh.vertex(u), tetmesh.vertex(v),
                          tetmesh.vertex(w)))
    {
      curr = intersection(q, p, u, v, w);

      ints.push_back(std::make_tuple(t, prev, curr));
    }

  } // t contains p
  // add last element because loop terminated
  ints.push_back(std::make_tuple(t, prev, p));
  return ints;
}

OVM::VertexHandle MetricField::linearNN(const Vec3d &_q)
{
  OVM::Vec3d q = toOVMVec3d(_q);
  OVM::VertexHandle nn;
  double prev_dist = std::numeric_limits<double>::max();

  for (OVM::VertexIter v_it = tetmesh.vertices_begin();
       v_it != tetmesh.vertices_end(); ++v_it)
  {
    OVM::VertexHandle curr = *v_it;
    OVM::Vec3d curr_pos = tetmesh.vertex(curr);
    double dist = (curr_pos - q).dot(curr_pos - q);
    if (dist < prev_dist)
    {
      prev_dist = dist;
      nn = curr;
    }
  }
  return nn;
}

// returns a list of all tets containing q
std::vector<OVM::CellHandle> MetricField::locateTets(const Vec3d &q)
{
  // check if point was hashed already
  if (hash.contains(q))
  {
    // std::cout << "contained in hash" << std::endl;
    return hash[q];
  }

  std::vector<OVM::CellHandle> tets;

  // check if previous cell contains point as good guess
  if (prevTetContainingPoint != OVM::TopologyKernel::InvalidCellHandle) {
    if (pointInTet(prevTetContainingPoint, q)) {
      tets.push_back(prevTetContainingPoint);
      // std::cout << " found in prev tet BIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIG " << std::endl;
      return tets;
    }
  }

  // std::cout << "searching q " << q << std::endl;

  // Locate which tet contains q
  OVM::VertexHandle s = linearNN(q);
  // std::cout << "Nearest neighbor" << tetmesh.vertex(s) << std::endl;

  // Test all incident cells to s if they contain q
  for (const OVM::CellHandle &ch : tetmesh.vertex_cells(s))
  {
    // std::cout << "inspecting cell " << ch << std::endl;

    if (pointInTet(ch, q))
    {
      // std::cout << "contained " << std::endl;
      prevTetContainingPoint = ch;
      tets.push_back(ch);
    }
  }

  if (tets.empty())
  {
    // std::cout << "did not find " << q.transpose() << " -------- fallback"
              // << std::endl;
    for (const auto &ch : tetmesh.cells())
    {
      if (pointInTet(ch, q))
      {
        prevTetContainingPoint = ch;
        tets.push_back(ch);
        // std::cout << "in tet " << ch << std::endl;
      }
    }

    // if still empty, then mesh doesn't cover the domain
    if (tets.empty())
    {
      std::cout << "something is off, mesh doesn't cover the domain" << std::endl;
      tets.push_back(OVM::CellHandle(-1));
    }
  }
  hash.try_emplace(q, tets);

  return tets;
}

// returns true if initialization was successful
bool MetricField::initialization(OVM::CellHandle &t, const Vec3d &q,
                                 const Vec3d &p, OVM::VertexHandle &u,
                                 OVM::VertexHandle &v, OVM::VertexHandle &w)
{
  bool continueSearch = true;
  // cycle through every halfface, check if intersects with ray qp
  for (const OVM::HalfFaceHandle &hfh : tetmesh.cell_halffaces(t))
  {
    if (!continueSearch)
    {
      break;
    }

    const auto vertices = tetmesh.get_halfface_vertices(hfh);
    if (robustRayTriangle(q, p, tetmesh.vertex(vertices[0]),
                          tetmesh.vertex(vertices[1]),
                          tetmesh.vertex(vertices[2])))
    {
      // ray intersects with triangle,
      // brute force combination until valid state is found, max 6 combinations
      int combinations[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
      for (int i = 0; i < 6; i++)
      {
        u = vertices[combinations[i][0]];
        v = vertices[combinations[i][1]];
        w = vertices[combinations[i][2]];
        if (validState(q, p, u, v, w))
        {
          continueSearch = false;
          break;
        }
      }
    }
  }
  // if search was broken, init was successful
  return !continueSearch;
}

bool MetricField::validState(const Vec3d &q, const Vec3d &p,
                             OVM::VertexHandle u, OVM::VertexHandle v,
                             OVM::VertexHandle w)
{
  // std::cout << "wvqp " << orient3dHelper(tetmesh, w,v,q,p) << std::endl;
  // std::cout << "vuqp " << orient3dHelper(tetmesh, v,u,q,p) << std::endl;
  // std::cout << "uwqp " << orient3dHelper(tetmesh, u,w,q,p) << std::endl;
  // std::cout << "intersects uvw " << robustRayTriangle(q, p,
  // tetmesh.vertex(u), tetmesh.vertex(v), tetmesh.vertex(w)) << std::endl;
  return orient3dHelper(w, v, q, p) >= 0 && orient3dHelper(v, u, q, p) >= 0 &&
         orient3dHelper(u, w, q, p) >= 0 &&
         robustRayTriangle(q, p, tetmesh.vertex(u), tetmesh.vertex(v),
                           tetmesh.vertex(w));
}

bool robustRayTriangle(const Vec3d &q, const Vec3d &p, const OVM::Vec3d &v1,
                       const OVM::Vec3d &v2, const OVM::Vec3d &v3)
{
  double orientQ = orient3d(q.data(), v1.data(), v2.data(), v3.data());
  double orientP = orient3d(p.data(), v1.data(), v2.data(), v3.data());
  double case1 = orient3d(p.data(), v1.data(), q.data(), v2.data());
  double case2 = orient3d(p.data(), v3.data(), v2.data(), q.data());
  double case3 = orient3d(p.data(), v1.data(), v3.data(), q.data());
  // std::cout << "ray triangle q" << orientQ << std::endl;
  // std::cout << "ray triangle p" << orientP << std::endl;
  // std::cout <<  "\t case1 " << case1 << std::endl;
  // std::cout <<  "\t case2 " << case2 << std::endl;
  // std::cout <<  "\t case3 " << case3 << std::endl;
  if ((orientQ > 0 && orientP > 0) ||
      (orientQ < 0 &&
       orientP < 0))
  { // p and q lie on same side of plane -> no intersection
    return false;
  }
  if (orientQ == 0.0 &&
      orientP ==
          0.0)
  { // p, q lie within the plane of the triange, reduced to 2d case
    // case1 = orient2d(v1.data(), q.data(), p.data());
    // case2 = orient2d(v2.data(), q.data(), p.data());
    // case3 = orient2d(v3.data(), q.data(), p.data());
    // TODO, need to project 3 dimensional vectors onto 2d plane,
    // intersection if atleast one orientation test is positive and one
    // negative. std::cout << "coplanar points, p q may lie within triangle
    // face" << std::endl;
    return false;
  }

  if (orientQ >= 0.0 && orientP <= 0.0)
  {
    if (case1 >= 0.0 && case2 >= 0.0 && case3 >= 0.0)
    {
      return true;
    } // intersects triangle
  }
  if (orientQ <= 0.0 && orientP >= 0.0)
  {
    if (case1 <= 0.0 && case2 <= 0.0 && case3 <= 0.0)
    {
      return true;
    } // intersects triangle
  }
  return false;
}

// finds the halfface with the 3 vertices in a given cell
OVM::HalfFaceHandle
MetricField::find_halfface_in_cell(std::vector<OVM::VertexHandle> &verts,
                                   OVM::CellHandle &ch)
{
  // compares every halfface to the list of vertices
  for (const OVM::HalfFaceHandle &hfh : tetmesh.cell_halffaces(ch))
  {
    if (std::is_permutation(verts.begin(), verts.end(),
                            tetmesh.get_halfface_vertices(hfh).begin()))
    {
      return hfh;
    }
  }
  return OVM::TopologyKernel::InvalidHalfFaceHandle;
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

  if (u < 0.0 || u > 1.0)
  {
    std::cout << "Did not intersect with this" << std::endl;
    std::cout << "\t u " << tetmesh.vertex(_u) << std::endl;
    std::cout << "\t v " << tetmesh.vertex(_v) << std::endl;
    std::cout << "\t w " << tetmesh.vertex(_w) << std::endl;
    std::cout << "1st, no intersection, this should not have been the case"
              << std::endl;
    return Vec3d(0.0, 0.0, 0.0);
  }

  q = s.cross(edge1);
  v = f * rayVector.dot(q);

  if (v < 0.0 || u + v > 1.0)
  {
    std::cout << "v " << v << " u " << u << std::endl;
    std::cout << "2nd, no intersection, this should not have been the case"
              << std::endl;
    return Vec3d(0.0, 0.0, 0.0);
  }

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
  std::vector<OVM::VertexHandle> vertices = tetmesh.get_cell_vertices(ch);
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
/* TODO

    - covering ovm mesh
    - attach metric to each vertex
    - extend tet class by metric field functionality
    - port tetFinder code
*/



/**
 * find tets containting p but with walking strategy from given vertex
*/
std::vector<OVM::CellHandle> MetricField::locateTetsFast(const OVM::VertexHandle &q, const Vec3d &p)
{
  std::cout << "start vertex q " << tetmesh.vertex(q) << std::endl;

  OVM::CellHandle t;
  OVM::VertexHandle u, v, w, s;
  OVM::HalfFaceHandle triangle;
  std::vector<OVM::VertexHandle> verts;

  /*
  // do initial assignments
  t = *tetmesh.vc_iter(q); // get start tet
  verts = tetmesh.get_cell_vertices(t, q);

  

  u = verts[1];
  v = verts[2];
  w = verts[3];
  assert (u != q && v != q && w != q);

  // initialization step
  if (orient3dHelper(v,u,q,p) > 0){

    while (orient3dHelper(w,u,q,p) > 0) {

      v = w;

      // jump to neighbor
      verts.clear();
      verts[0] = q;
      verts[1] = u;
      verts[2] = w;

      triangle = find_halfface_in_cell(verts, t);
      triangle = tetmesh.opposite_halfface_handle(triangle);
      t = tetmesh.incident_cell(triangle);

      // reassign w
      const auto candidates = tetmesh.get_cell_vertices(t);
      for (const OVM::VertexHandle &vh : candidates)
      {
        if (vh != u && vh != v && vh != q)
        {
          w = vh;
          break;
        }
      }
    }    
  } else {

    do
    {
      w = v;

      // jump to neighbor
      verts.clear();
      verts[0] = q;
      verts[1] = u;
      verts[2] = v;

      triangle = find_halfface_in_cell(verts, t);
      triangle = tetmesh.opposite_halfface_handle(triangle);
      t = tetmesh.incident_cell(triangle);

      // reassign v
      const auto candidates = tetmesh.get_cell_vertices(t);
      for (const OVM::VertexHandle &vh : candidates)
      {
        if (vh != u && vh != w && vh != q)
        {
          v = vh;
          break;
        }
      }
      
    } while (orient3dHelper(v, u, q, p) < 0);
  }
  
  // now v and w lie on opposite side of plane uqp
  // vuqp is positively oriented and wuqp negatively
  while (orient3dHelper(v,w,q,p) > 0) {
    // jump to neighbor
    verts.clear();
    verts[0] = q;
    verts[1] = v;
    verts[2] = w;

    triangle = find_halfface_in_cell(verts, t);
    triangle = tetmesh.opposite_halfface_handle(triangle);
    t = tetmesh.incident_cell(triangle);

    // reassign s
    const auto candidates = tetmesh.get_cell_vertices(t);
    for (const OVM::VertexHandle &vh : candidates)
    {
      if (vh != v && vh != w && vh != q)
      {
        s = vh;
        break;
      }
    }
    if (orient3dHelper(s,u,q,p) > 0) v = s;
    else w = s;
  }

  // reassign u
  const auto candidates = tetmesh.get_cell_vertices(t);
  for (const OVM::VertexHandle &vh : candidates)
  {
    if (vh != v && vh != w && vh != q)
    {
      u = vh;
      break;
    }
  }
  // end of initialization step
  */

  initializeFast(q,p,u,v,w, t);
  {// debug 
    TM DEBUG_MESH;
    std::vector<OVM::VertexHandle> faceseseses;
    faceseseses.push_back(DEBUG_MESH.add_vertex(tetmesh.vertex(u)));
    faceseseses.push_back(DEBUG_MESH.add_vertex(tetmesh.vertex(v)));
    faceseseses.push_back(DEBUG_MESH.add_vertex(tetmesh.vertex(w)));
    DEBUG_MESH.add_face(faceseseses);
    saveToFile(DEBUG_MESH, "DEBUG_TRIANGLE.ovm");
  }


  // qp intersects triangle u,v,w
  // wvqp, vuqp, u,w,q,p are positively oriented
  while (orient3dHelper(u, w, v, p) > 0)
  {
    verts.clear();
    verts.push_back(u);
    verts.push_back(v);
    verts.push_back(w);
    triangle = find_halfface_in_cell(verts, t);

    // jump to neighbor
    triangle = tetmesh.opposite_halfface_handle(triangle);
    t = tetmesh.incident_cell(triangle);

    // reassign s
    const auto candidates = tetmesh.get_cell_vertices(t);
    for (const OVM::VertexHandle &vh : candidates)
    {
      if (vh != u && vh != v && vh != w)
      {
        s = vh;
        break;
      }
    }
    if (orient3dHelper(u, s, q, p) > 0)
    { // qp does not intersect triangle usw
      if (orient3dHelper(v, s, q, p) > 0)
      { // qp intersects triangle vsw
        u = s;
      }
      else
      { // qp intersects triangle usv
        w = s;
      }
    }
    else
    { // qp doesnt not intersect triangle usv
      if (orient3dHelper(w, s, q, p) > 0)
      { // qp intersects triangle usw
        v = s;
      }
      else
      { // qp intersects triangle vsw
        u = s;
      }
    }

  } // t contains p
  std::vector<OVM::CellHandle> tets;
  tets.push_back(t);

  std::cout << "point is contained in " << t << std::endl;

  OVM::CellPropertyT<bool> T = tetmesh.create_persistent_cell_property<bool>("contained").value();
  T[t] = true;
  saveToFile(tetmesh, "contained.ovm");
  return tets;

} 


/**
 * takes in Vertex q, assigns u,v,w and CellHandle t
 * such that it is a valid configuration
*/
void MetricField::initializeFast(const OVM::VertexHandle q, const Vec3d p, OVM::VertexHandle &u, OVM::VertexHandle &v, OVM::VertexHandle &w, OVM::CellHandle &t){
  for (auto vc_it = tetmesh.vc_iter(q); vc_it.is_valid(); ++vc_it){ 
    auto ch = *vc_it;

    auto verts = tetmesh.get_cell_vertices(ch, q);
    u = verts[1];
    v = verts[2];
    w = verts[3];

    // check if ray qp points into cell
    bool pointsIn = orient3dHelper(v, u, q, p) <= 0.0 &&
    orient3dHelper(u, w, q, p) <= 0.0 &&
    orient3dHelper(w, v, q, p) <= 0.0;

    // qp intersects triangle u,v,w
    // wvqp, vuqp, u,w,q,p are positively oriented
    if (pointsIn){
      t = ch;
      std::swap(u,w); // paper uses different convention 
      return;
    }
  }
  assert(false);

  std::cerr << " no valid configuration found, wtf ?" << std::endl;
  return;
}


