#include <Eigen/Dense>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>
#include "OpenVolumeMesh/IO/WriteOptions.hh"
#include "OpenVolumeMesh/IO/ReadOptions.hh"
#include "OpenVolumeMesh/IO/PropertyCodecs.hh"
#include "OpenVolumeMesh/IO/PropertyCodecsEigen.hh"
#include <OpenVolumeMesh/IO/ovmb_write.hh>
#include <OpenVolumeMesh/IO/ovmb_read.hh>
 
#include <fstream>
#include <cassert>
#include <list>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>

#include <vector>
#include <iostream>


// Typedefs for CGAL
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
// typedef CGAL::Mesh_triangulation_3<K>::type Delaunay; // Mesh_3-compatible triangulation

typedef Delaunay::Cell_handle Cell_handle;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Delaunay> C3T3;


using Vec3d = Eigen::Matrix<double, 3, 1>;
using Mat3d = Eigen::Matrix<double, 3, 3>;

namespace OVM = OpenVolumeMesh;

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

void readInTetmesh(TM &tetmesh, const std::string &filename)
{
  if (filename.ends_with(".ovmb")){
    readFile(tetmesh, filename); // read ovmb mesh
  } else if (filename.ends_with(".ovm")){
    OVM::IO::FileManager fm;
    fm.readFile(filename, tetmesh, true, true);
  } else {
    std::cerr << "unknown file extension while reading" << std::endl;
    exit(0);
  }
}


void populateTestMesh(TM &tetmesh){
    auto v1 = tetmesh.add_vertex(OVM::Vec3d(0,0,0));
    auto v2 = tetmesh.add_vertex(OVM::Vec3d(1,0,0));
    auto v3 = tetmesh.add_vertex(OVM::Vec3d(0,1,0));
    auto v4 = tetmesh.add_vertex(OVM::Vec3d(1,1,0));
    auto v5 = tetmesh.add_vertex(OVM::Vec3d(0.5,0.5,0.5));

    tetmesh.add_cell(v2,v1,v5,v4);
    tetmesh.add_cell(v3,v1,v4,v5);
}

Point insert_from_degenerate(const OVM::Vec3d &v1, const OVM::Vec3d &v2, const OVM::Vec3d &v3 ){
    auto dir1 = v1 - v2;
    auto dir2 = v1 - v3;
    auto normal = dir1.cross(v2).normalize();

    auto barycenter = (v1 + v2 + v3)/3.0;

    double dist = (v1 - barycenter).norm();

    auto new_p = barycenter + dist * normal;
    Point p = Point(new_p[0], new_p[1], new_p[2]);
    return p;

}

bool mark_degenerate(TM &tetmesh, std::vector<Point> &additional ){
    // first, clear previous points
    additional.clear();

    // marks degenrate cells with a true based on volume being near zero

    const auto volume = [&tetmesh](OVM::VertexHandle v1, OVM::VertexHandle v2, OVM::VertexHandle v3, OVM::VertexHandle v4) -> double{
        Vec3d p1 = Vec3d(tetmesh.vertex(v1).data());
        Vec3d p2 = Vec3d(tetmesh.vertex(v2).data());
        Vec3d p3 = Vec3d(tetmesh.vertex(v3).data());
        Vec3d p4 = Vec3d(tetmesh.vertex(v4).data());

        Mat3d mat;
        mat << p1-p4, p2-p4, p3-p4;

        double res = 1./6. * mat.determinant();
        // std::cout << "vol " << res << std::endl;
        return res;
    };
    auto volume_prop = tetmesh.request_cell_property<double>("volume", 0);
    auto deg_prop = tetmesh.request_cell_property<bool>("degenerate", false);
    tetmesh.set_persistent(volume_prop);
    tetmesh.set_persistent(deg_prop);

    int counter = 0;

    for (auto c_it = tetmesh.cells_begin(); c_it != tetmesh.cells_end(); ++c_it){
        auto cell = *c_it;
        auto verts = tetmesh.get_cell_vertices(cell);
        double vol = std::abs(volume(verts[0], verts[1], verts[2], verts[3]));
        volume_prop[cell] = vol;
        if (vol < 1e-8){
            deg_prop[cell] = true;
            counter++;
            // tetmesh.delete_cell(cell);
            additional.push_back(insert_from_degenerate(tetmesh.vertex(verts[0]), tetmesh.vertex(verts[1]), tetmesh.vertex(verts[2])));
            std::cout << "vert 1 " << tetmesh.vertex(verts[0]) << std::endl;
            std::cout << "vert 2 " << tetmesh.vertex(verts[1]) << std::endl;
            std::cout << "vert 3 " << tetmesh.vertex(verts[2]) << std::endl;
            std::cout << "vert 4 " << tetmesh.vertex(verts[3]) << std::endl;
        }
    }

    std::cout << "Degenerate cells counter " << counter << std::endl;

    

}

// element wise maximum
OVM::Vec3d vec3_max(OVM::Vec3d a, OVM::Vec3d b){
    return OVM::Vec3d(std::max(a[0], b[0]), std::max(a[1], b[1]), std::max(a[2], b[2]));
}


// element wise minimum
OVM::Vec3d vec3_min(OVM::Vec3d a, OVM::Vec3d b){
    return OVM::Vec3d(std::min(a[0], b[0]), std::min(a[1], b[1]), std::min(a[2], b[2]));
}

template<typename V2, typename V1>
    V2 toVec(const V1 &v){
        return V2(v[0], v[1], v[2]);
    };

std::vector<Point> get_bounds(TM &tetmesh){
    OVM::Vec3d min = OVM::Vec3d(__DBL_MAX__);
    OVM::Vec3d max = OVM::Vec3d(__DBL_MIN__);

    for (auto vh : tetmesh.vertices()){
        auto point = tetmesh.vertex(vh);
        min = vec3_min(point, min);
        max = vec3_max(point, max);
    }

    OVM::Vec3d center = (max+min)/2.0;

    // std::cout << "center " << center << std::endl;

    double max_dim = (max-min).max();
    double bounds_hs = max_dim * 0.7;

    OVM::Vec3d hs_diag = OVM::Vec3d(bounds_hs); // all same value

    //AABB bounds
    OVM::Vec3d minVec = center - hs_diag;
    OVM::Vec3d maxVec = center + hs_diag;

    std::vector<Point> bounding_points;
    
    

    bounding_points.push_back(Point(minVec[0], minVec[1], minVec[2])); // 0 0 0
    bounding_points.push_back(Point(minVec[0], minVec[1], maxVec[2])); // 0 0 1
    bounding_points.push_back(Point(minVec[0], maxVec[1], minVec[2])); // 0 1 0
    bounding_points.push_back(Point(minVec[0], maxVec[1], maxVec[2])); // 0 1 1
    bounding_points.push_back(Point(maxVec[0], minVec[1], minVec[2])); // 1 0 0
    bounding_points.push_back(Point(maxVec[0], minVec[1], maxVec[2])); // 1 0 1
    bounding_points.push_back(Point(maxVec[0], maxVec[1], minVec[2])); // 1 1 0
    bounding_points.push_back(Point(maxVec[0], maxVec[1], maxVec[2])); // 1 1 1 

    std::array<std::array<std::array<Vec3d,2>,2>,2> box;
    box[0][0][0] = toVec<Vec3d>(bounding_points.at(0));
    box[0][0][1] = toVec<Vec3d>(bounding_points.at(1));
    box[0][1][0] = toVec<Vec3d>(bounding_points.at(2));
    box[0][1][1] = toVec<Vec3d>(bounding_points.at(3));
    box[1][0][0] = toVec<Vec3d>(bounding_points.at(4));
    box[1][0][1] = toVec<Vec3d>(bounding_points.at(5));
    box[1][1][0] = toVec<Vec3d>(bounding_points.at(6));
    box[1][1][1] = toVec<Vec3d>(bounding_points.at(7));

    auto b = box;

    Vec3d f1 = (b[0][0][0] + b[0][1][0] + b[0][1][1] + b[0][0][1]) / 4.0;
    Vec3d f2 = (b[1][0][0] + b[1][1][0] + b[1][1][1] + b[1][0][1]) / 4.0;
    Vec3d f3 = (b[0][0][0] + b[0][0][1] + b[1][0][1] + b[1][0][0]) / 4.0;
    Vec3d f4 = (b[0][1][0] + b[0][1][1] + b[1][1][1] + b[1][1][0]) / 4.0;
    Vec3d f5 = (b[0][0][0] + b[0][1][0] + b[1][1][0] + b[1][0][0]) / 4.0;
    Vec3d f6 = (b[0][0][1] + b[0][1][1] + b[1][1][1] + b[1][0][1]) / 4.0;

    bounding_points.push_back(toVec<Point>(f1));
    bounding_points.push_back(toVec<Point>(f2));
    bounding_points.push_back(toVec<Point>(f3));
    bounding_points.push_back(toVec<Point>(f4));
    bounding_points.push_back(toVec<Point>(f5));
    bounding_points.push_back(toVec<Point>(f6));

    return bounding_points;
    
}


// void refine_mesh(Delaunay &delaunay){
//     CGAL::Mesh_criteria_3<Delaunay> criteria(
//         CGAL::parameters::cell_radius_edge_ratio(2.0)
//     );
//     C3T3 c3t3 = CGAL::make_mesh_3<C3T3>(delaunay, criteria);

//     delaunay = c3t3.triangulation();
// }

int main(int argc, char *argv[]) {
    std::cout << CGAL_VERSION_MAJOR << "." << CGAL_VERSION_MINOR << "." << CGAL_VERSION_PATCH << std::endl;


    std::string input_file;
    std::string output_file;
    if (argc != 3){
        std::cerr << "Usage: ./input_handler path/to/input.ovm path/to/output.ovm";
    } else {
        input_file = argv[1];
        output_file = argv[2];
    }


    // Step 1: Read the tetrahedral mesh from OpenVolumeMesh
    TM openMesh;
    // std::string input_file = "../../cube_field/meshes/cube_1.ovm";

    readInTetmesh(openMesh, input_file);
    // populateTestMesh(openMesh);
    saveToFile(openMesh, "mesh_before.ovm");

    std::cout << "Loaded in OpenVolumeMesh with " << openMesh.n_vertices()
              << " vertices and " << openMesh.n_cells() << " tetrahedra." << std::endl;

    // Step 2: Convert OpenVolumeMesh vertices to CGAL points
    std::vector<Point> cgal_points;
    std::map<OpenVolumeMesh::VertexHandle, int> ovm_to_cgal_map;
    int id = 0;

    for (auto vh : openMesh.vertices()) {
        auto point = openMesh.vertex(vh);
        cgal_points.push_back(Point(point[0], point[1], point[2]));
        ovm_to_cgal_map[vh] = id++;
    }

    std::cout << "id " << id << std::endl;

    // Step 3: Insert points into CGAL Delaunay triangulation
    Delaunay delaunay;
    delaunay.insert(cgal_points.begin(), cgal_points.end());

    std::cout << "CGAL before triangulation has " << delaunay.number_of_vertices()
              << " vertices and " << delaunay.number_of_cells() << " cells." << std::endl;

    // Add additional points around the mesh (points of the bounding box)
    auto bounds = get_bounds(openMesh);
    delaunay.insert(bounds.begin(), bounds.end());
    

    std::cout << "CGAL triangulation has " << delaunay.number_of_vertices()
              << " vertices and " << delaunay.number_of_cells() << " cells." << std::endl;
    
    TM newOpenMesh;

    // STEP REFINEMENT OF MESH
    std::vector<Point> additional_points;
    bool added_points = false;
    do {
        added_points = false;

        // Step 4: Convert CGAL triangulation back to OpenVolumeMesh
        newOpenMesh.clear();

        // Map CGAL vertices to new OpenVolumeMesh vertex handles
        std::map<Point, OpenVolumeMesh::VertexHandle> cgal_to_ovm_map;

        for (auto v_it = delaunay.finite_vertices_begin(); v_it != delaunay.finite_vertices_end(); ++v_it) {

            Point p = (*v_it).point();
            auto vh = newOpenMesh.add_vertex(OpenVolumeMesh::Geometry::Vec3d(p.x(), p.y(), p.z()));
            cgal_to_ovm_map[p] = vh;
        }

        // Add tetrahedra to the OpenVolumeMesh
        for (auto cell_it = delaunay.finite_cells_begin(); cell_it != delaunay.finite_cells_end(); ++cell_it) {
            std::vector<OpenVolumeMesh::VertexHandle> vhandles;
            for (int i = 0; i < 4; ++i) {

                Point p = cell_it->vertex(i)->point();
                vhandles.push_back(cgal_to_ovm_map[p]);
            }
            newOpenMesh.add_cell(vhandles);
        }

        

        mark_degenerate(newOpenMesh, additional_points);
        std::cout << "Converted back to OpenVolumeMesh with " << newOpenMesh.n_vertices()
                << " vertices and " << newOpenMesh.n_cells() << " tetrahedra." << std::endl;


        if (additional_points.size() > 0) {
            added_points = true;
            // add new points to resolve degeneracy
            delaunay.insert(additional_points.begin(), additional_points.end());
            additional_points.clear();
        }
    } while (added_points);
    

    // Step 5: Write the new mesh to file
    // std::string outputMeshFile = "output_mesh.ovm";
    saveToFile(newOpenMesh, output_file);

    std::cout << "New mesh written to " << output_file << std::endl;
    

    return EXIT_SUCCESS;
}
