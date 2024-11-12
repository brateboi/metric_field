#ifndef MESHELEMENT_HH
#define MESHELEMENT_HH

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

namespace OVM = OpenVolumeMesh;


enum MeshElementType : char {
    ELEMENT_NULL_TYPE = -1,
    ELEMENT_VERTEX_TYPE = 0,
    ELEMENT_EDGE_TYPE = 1,
    ELEMENT_FACE_TYPE = 2,
    ELEMENT_CELL_TYPE = 3
};

class MeshElement {
public:
    MeshElement(MeshElementType type, int id) {set(type,id);}
    MeshElement() : MeshElement(ELEMENT_NULL_TYPE, -1) {}
    MeshElement(OVM::VertexHandle vh) {set(vh);}
    MeshElement(OVM::EdgeHandle eh) {set(eh);}
    MeshElement(OVM::FaceHandle fh) {set(fh);}
    MeshElement(OVM::CellHandle ch) {set(ch);}
    MeshElement(OVM::HalfFaceHandle hfh) {set(hfh);}
    MeshElement(OVM::HalfEdgeHandle heh) {set(heh);}

    inline void set(MeshElementType type, int idx) {
        this->_type=type;
        this->_idx=idx;
    }

    inline void set(const OVM::VertexHandle& vh) {set(ELEMENT_VERTEX_TYPE,_idx=vh.idx());}
    inline void set(const OVM::EdgeHandle& eh) {set(ELEMENT_EDGE_TYPE,eh.halfedge_handle(0).idx());}
    inline void set(const OVM::FaceHandle& fh) {set(ELEMENT_FACE_TYPE,fh.halfface_handle(0).idx());}
    inline void set(const OVM::CellHandle& ch) {set(ELEMENT_CELL_TYPE,ch.idx());}
    inline void set(const OVM::HalfFaceHandle& hfh) {set(ELEMENT_FACE_TYPE,hfh.idx());}
    inline void set(const OVM::HalfEdgeHandle& heh) {set(ELEMENT_EDGE_TYPE,heh.idx());}
    inline void set(const MeshElement& elem) {set(elem._type,elem._idx);}

    inline bool is_valid() const {return _type!=ELEMENT_NULL_TYPE&&_idx!=-1;}
    inline bool is_vertex() const {return _type==ELEMENT_VERTEX_TYPE;}
    inline bool is_edge() const {return _type==ELEMENT_EDGE_TYPE;}
    inline bool is_face() const {return _type==ELEMENT_FACE_TYPE;}
    inline bool is_cell() const {return _type==ELEMENT_CELL_TYPE;}

    inline OVM::VertexHandle vh() const {return OVM::VertexHandle(is_vertex()? _idx : -1);}
    inline OVM::HalfEdgeHandle heh() const {return OVM::HalfEdgeHandle(is_edge()? _idx : -1);}
    inline OVM::HalfFaceHandle hfh() const {return OVM::HalfFaceHandle(is_face()? _idx : -1);}
    inline OVM::EdgeHandle eh() const {return is_edge()? OVM::HalfEdgeHandle(_idx).edge_handle() : OVM::EdgeHandle(-1);}
    inline OVM::FaceHandle fh() const {return is_face()? OVM::HalfFaceHandle(_idx).face_handle() : OVM::FaceHandle(-1);}
    inline OVM::CellHandle ch() const {return OVM::CellHandle(is_cell()? _idx : -1);}

    inline constexpr char dim() const {return _type;}
    inline constexpr int idx() const {return _idx;}
    inline constexpr MeshElementType type() const {return _type;}

    constexpr bool operator<(const MeshElement& elem) const { return (this->dim() < elem.dim()); }
    constexpr bool operator>(const MeshElement& elem) const { return (this->dim() > elem.dim()); }
    constexpr bool operator==(const MeshElement& elem) const { return elem._idx == this->_idx && this->_type == elem._type; }
    constexpr bool operator!=(const MeshElement& elem) const { return elem._idx != this->_idx || this->_type != elem._type;  }

private:
    MeshElementType _type;
    int _idx;

}; // class Element

#endif // MESHELEMENT_HH
