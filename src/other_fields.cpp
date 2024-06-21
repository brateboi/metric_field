#include <Eigen/Dense>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

#include "OpenVolumeMesh/IO/WriteOptions.hh"
#include "OpenVolumeMesh/IO/ReadOptions.hh"
#include "OpenVolumeMesh/IO/PropertyCodecs.hh"

#include "metric_field.hpp"
#include "other_fields.hpp"

#include <cmath>
#include <limits>
#include "OpenVolumeMesh/IO/PropertyCodecsEigen.hh"
#include <OpenVolumeMesh/IO/ovmb_write.hh>
#include <OpenVolumeMesh/IO/ovmb_read.hh>

using Quaternion = Eigen::Quaterniond;
using Vec3d = Eigen::Matrix<double, 3, 1>;
using Vec4d = Eigen::Matrix<double, 4, 1>;
using Vec9d = Eigen::Matrix<double, 9, 1>;
using Mat3d = Eigen::Matrix<double, 3, 3>;
using Mat4d = Eigen::Matrix<double, 4, 4>;
using Mat9d = Eigen::Matrix<double, 9, 9>;

using TM = OVM::TetrahedralGeometryKernel<OVM::Vec3d,
                                          OVM::TetrahedralMeshTopologyKernel>;


void saveToFile(TM &tetmesh, std::string filename) {
    // store mesh to .ovm to view it in OpenFlipper
    OVM::IO::FileManager fileManager;
    // Store mesh to file "myMesh.ovm" in the current directory
    if (!fileManager.writeFile(filename, tetmesh)) {
        std::cout << "somethings wrong" << std::endl;
    };
}

void saveToFileOVMB(TM &tetmesh, std::string filename){
    auto wo = (OVM::IO::WriteOptions());
    auto propCodecs = (OVM::IO::g_default_property_codecs);
    OVM::IO::register_eigen_codecs(propCodecs);
    std::ofstream off(filename.c_str(), std::ios::binary);
    OVM::IO::ovmb_write(off, tetmesh, wo, propCodecs);
    off.close();
}

void readFile(TM &mesh, std::string filename) {
    auto ro = (OVM::IO::ReadOptions());
    auto propCodecs = (OVM::IO::g_default_property_codecs);
    OVM::IO::register_eigen_codecs(propCodecs);
    std::ifstream iff(filename.c_str(), std::ios::binary);
    OVM::IO::ovmb_read(iff, mesh, ro, propCodecs);
    iff.close();
}


// analytical root metric
Mat3d rootMetricAnalytical(Vec3d &p)
{
  //arc example --------
  double x = p(0);
  double y = p(1);
  double z = p(2);


  double x2y2 = x * x + y * y;
  double temp = sqrt(1 / x2y2);
  
  return (Mat3d() << (x * x + y * y * temp) / x2y2, x * y * (1 - temp) / x2y2, 0,
          x * y * (1 - temp) / x2y2, (x * x * temp + y * y) / x2y2, 0,
          0, 0, 1)
      .finished();

  // twisted rod example------------
  // double x = p(0);
  // double y = p(1);
  // double z = p(2);

  // Mat3d metric;
  // metric << 1, -z, 0,
  //   -z, x*x+z*z+1, x,
  //   0, x, 1;
  
  // Eigen::SelfAdjointEigenSolver<Mat3d> es(metric);
  // Mat3d rootmetric = es.operatorSqrt();
  // return rootmetric;

}


Vec3d deformInverse(Vec3d &p) {
  Vec3d result;
  result << sqrt(p(0)*p(0)+p(1)*p(1)), atan2(p(1), p(0)), p(2); // arc example



  // twisted rod example
  // double x = p[0], y = p[1], z = p[2];
  // result << (cos(y)*x-sin(y)*z) , y, (sin(y)*x + cos(y)*z); 

  return result;
}

Mat3d pushforward(Vec3d &p_){
  Vec3d p = deformInverse(p_);

  Mat3d pf;

  // arc example pushforward
  double r = p(0), theta = p(1), z = p(2);
  double cosT = cos(theta);
  double sinT = sin(theta);
  pf.setZero();
  pf(0,0) = cosT;
  pf(0,1) = -r*sinT;
  pf(1,0) = sinT;
  pf(1,1) = r*cosT;
  pf(2,2) = 1.0;

  // twisted rod example ---------------
  // pf.setZero();
  // double x = p(0),y = p(1),z = p(2);
  // double sinY = sin(y), cosY = cos(y);
  // pf << cosY, z * cosY - x * sinY, sinY,
  //   0, 1, 0,
  //   -sinY, -x* cosY - z * sinY, cosY;



  return pf;
}

Mat3d rotationPushforward(Vec3d &p1, Vec3d &p2) {
  return pushforward(p1).inverse() * rootMetricAnalytical(p1).inverse() * rootMetricAnalytical(p2) * pushforward(p2);
}

// curl analytical for arc example
Vec9d curlAnalytical(Vec3d p){
  Vec9d curl;
  curl.setZero();
  double x = p(0);
  double y = p(1);
  double z = p(2);


  curl[2] = y/(x*x+y*y);
  curl[5] = -x/(x*x+y*y);

  return curl;
}

/**
 * Evaluates the coefficients of the one-forms W through a 9x9 linear system
 *
 * @param Vec3d point
 * @returns Matrix3d W
 */
Mat3d eval_W_Analytical(Vec3d _p)
{  

  Vec9d curl = curlAnalytical(_p);
  Mat3d metric = rootMetricAnalytical(_p);

  return unstack(A_x(metric).partialPivLu().solve(curl));
}

/**
 * integrates with trapezoidal rule between point a and b, analytical
 * not weighted as we don't have path lengths
 *
 * @param Vec3d a
 * @param Vec3d b
 * @returns Vec3d integrated vector
 */
Vec3d integrate2PointAnalytical(Vec3d a, Vec3d b)
{
  Vec3d res = 0.5 * (eval_W_Analytical(a) + eval_W_Analytical(b)).transpose() * (b - a);
  return res;
}

// frobenius norm
Mat3d recursiveDivideAnalytical(const Vec3d &a, const Vec3d &b, int depth )
{
  Mat3d r1 = lie_exp(integrate2PointAnalytical(a, b));
  Vec3d midpoint = (a + b) / 2.0;
  double length = (a - b).squaredNorm();
  // double length = (a - b).norm();
  // double length = 1.0;
  if (((r1 - lie_exp(integrate2PointAnalytical(a, midpoint)) * lie_exp(integrate2PointAnalytical(midpoint, b))).squaredNorm()) /
          length < 10e-6 || depth == 0)
  {
    return r1;
  }
  else
  {
    return recursiveDivideAnalytical(a, midpoint, --depth) *
           recursiveDivideAnalytical(midpoint, b, --depth);
  }
}

Mat3d computeCoeffAnalytical(const Vec3d &q, const Vec3d &p, int maxdepth = 10)
{
  return recursiveDivideAnalytical(q, p, maxdepth);
}


/**
 * calculates the curl and W_a, W_b with the constant discretization and writes it in
 * @param curl
 * @param W_a
 * @param W_b
*/
void calculate_curl_constant_and_W(const OVM::CellHandle a, const OVM::CellHandle b, TM &tetmesh, Mat3d &curl_constant, Mat3d &W_a, Mat3d &W_b) {
  auto face = commonFace(a,b,tetmesh);

  // assert that A and B are actually neighbors
  assert(face.is_valid());

  Vec3d n_f = toVec3d(tetmesh.normal(face.halfface_handle(0)));

  Vec3d c_a = toVec3d(tetmesh.barycenter(a));
  Vec3d c_b = toVec3d(tetmesh.barycenter(b));
  Vec3d c_f = toVec3d(tetmesh.barycenter(face));

  double s_a = ((c_f - c_a).dot(n_f)) / ((c_b - c_a).dot(n_f));
  double s_b = 1.0 - s_a;

  // eval_W_constant
  OVM::CellPropertyT<Mat3d> metric = tetmesh.request_cell_property<Mat3d>("constant_metric");

  Mat3d grad_A_scalars = (metric[b] - metric[a]) / (c_b - c_a).dot(n_f);
  
  Mat3d grad_A_1 = grad_A_scalars.col(0) * n_f.transpose();
  Mat3d grad_A_2 = grad_A_scalars.col(1) * n_f.transpose();
  Mat3d grad_A_3 = grad_A_scalars.col(2) * n_f.transpose();

  curl_constant(0,0) = grad_A_1(2,1) - grad_A_1(1,2);
  curl_constant(1,0) = grad_A_1(0,2) - grad_A_1(2,0);
  curl_constant(2,0) = grad_A_1(1,0) - grad_A_1(0,1);

  curl_constant(0,1) = grad_A_2(2,1) - grad_A_2(1,2);
  curl_constant(1,1) = grad_A_2(0,2) - grad_A_2(2,0);
  curl_constant(2,1) = grad_A_2(1,0) - grad_A_2(0,1);

  curl_constant(0,2) = grad_A_3(2,1) - grad_A_3(1,2);
  curl_constant(1,2) = grad_A_3(0,2) - grad_A_3(2,0);
  curl_constant(2,2) = grad_A_3(1,0) - grad_A_3(0,1);

  // W_a = unstack(A_x(metric[a]).partialPivLu().solve(stack(curl_constant)));
  W_a = unstack(A_x(metric[a]).partialPivLu().solve(stack(4.5*curl_constant))); // 5 for rod, 4.5 for arc
  // W_b = unstack(A_x(metric[b]).partialPivLu().solve(stack(curl_constant)));
  W_b = unstack(A_x(metric[b]).partialPivLu().solve(stack(4.5*curl_constant)));
}

// computesCoeff constant between two adjacent cells
Mat3d computeCoeffConstant(const OVM::CellHandle a, const OVM::CellHandle b, TM &tetmesh) {

  auto face = commonFace(a,b,tetmesh);

  // assert that A and B are actually neighbors
  assert(face.is_valid());

  Vec3d n_f = toVec3d(tetmesh.normal(face.halfface_handle(0)));

  Vec3d c_a = toVec3d(tetmesh.barycenter(a));
  Vec3d c_b = toVec3d(tetmesh.barycenter(b));
  Vec3d c_f = toVec3d(tetmesh.barycenter(face));

  double s_a = ((c_f - c_a).dot(n_f)) / ((c_b - c_a).dot(n_f));
  double s_b = 1.0 - s_a;
  Mat3d W_a, W_b;

  // eval_W_constant
  OVM::CellPropertyT<Mat3d> metric = tetmesh.request_cell_property<Mat3d>("constant_metric");

  Mat3d grad_A_scalars = (metric[b] - metric[a]) / (c_b - c_a).dot(n_f);
  Mat3d curl_constant;
  
  Mat3d grad_A_1 = grad_A_scalars.col(0) * n_f.transpose();
  Mat3d grad_A_2 = grad_A_scalars.col(1) * n_f.transpose();
  Mat3d grad_A_3 = grad_A_scalars.col(2) * n_f.transpose();

  curl_constant(0,0) = grad_A_1(2,1) - grad_A_1(1,2);
  curl_constant(1,0) = grad_A_1(0,2) - grad_A_1(2,0);
  curl_constant(2,0) = grad_A_1(1,0) - grad_A_1(0,1);

  curl_constant(0,1) = grad_A_2(2,1) - grad_A_2(1,2);
  curl_constant(1,1) = grad_A_2(0,2) - grad_A_2(2,0);
  curl_constant(2,1) = grad_A_2(1,0) - grad_A_2(0,1);

  curl_constant(0,2) = grad_A_3(2,1) - grad_A_3(1,2);
  curl_constant(1,2) = grad_A_3(0,2) - grad_A_3(2,0);
  curl_constant(2,2) = grad_A_3(1,0) - grad_A_3(0,1);

  W_a = unstack(A_x(metric[a]).partialPivLu().solve(stack(curl_constant)));
  // W_a = unstack(A_x(metric[a]).partialPivLu().solve(stack(4.5*curl_constant))); // 5 for rod, 4.5 for arc
  W_b = unstack(A_x(metric[b]).partialPivLu().solve(stack(curl_constant)));
  // W_b = unstack(A_x(metric[b]).partialPivLu().solve(stack(4.5*curl_constant)));


  // weighted trapezoidal rule
  Mat3d rot = lie_exp((s_a * W_a + s_b * W_b).transpose() * (c_b - c_a));
  
  return rot;
  
}





void attachConstantMetric(TM &tetmesh) {
  // Mat3d identity = Mat3d::Identity();
  OVM::CellPropertyT<Mat3d> metric = tetmesh.create_persistent_cell_property<Mat3d>("constant_metric").value();
  for (OVM::CellIter c_it = tetmesh.cells_begin(); c_it != tetmesh.cells_end(); ++c_it)
  {
    //------------arc example
    // sample metric at barycenter of cell
    Vec3d temp = toVec3d(tetmesh.barycenter(*c_it));
    metric[*c_it] = rootMetricAnalytical(temp);
  }
}




OVM::FaceHandle commonFace(const OVM::CellHandle a, const OVM::CellHandle b, TM &tetmesh){
  for (const auto hfh1 : tetmesh.cell_halffaces(a)) {
    for (const auto hfh2 : tetmesh.cell_halffaces(b)) {
      if (hfh1.face_handle() == hfh2.face_handle()) {
        return hfh1.face_handle();
      }
    }
  }
  std::cout << " weird " << a << std::endl;
  std::cout << "aswell" << b << std::endl;
  return OVM::FaceHandle(-1);
}