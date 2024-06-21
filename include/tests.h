#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>

#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh> 
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>
#include <OpenVolumeMesh/FileManager/FileManager.hh>

#include <Eigen/Dense>
#include <typeinfo>
#include "math.h"
#include "metric_field.hpp"
#include "other_fields.hpp"

#include <iostream>
using Eigen::Matrix;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::MatrixXd;

namespace OVM = OpenVolumeMesh;

using TM = OVM::TetrahedralGeometryKernel<OVM::Vec3d, OVM::TetrahedralMeshTopologyKernel>;

void addTet(TM &tetmesh, OVM::VertexHandle &v0, OVM::VertexHandle &v1, OVM::VertexHandle &v2, OVM::VertexHandle &v3);
void populateMesh(TM &tetmesh);
void normalCase(MetricField &field);
void throughVertex(MetricField &field);
void throughEdge(MetricField &field);
void throughFace(MetricField &field);
void startPointOnVertex(MetricField &field);
void startPointOnEdge(MetricField &field);
void startPointOnFace(MetricField &field);
void endPointOnVertex(MetricField &field);
void endPointOnEdge(MetricField &field);
void endPointOnFace(MetricField &field);
void showPoint();
