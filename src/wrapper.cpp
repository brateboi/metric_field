#include <metric_field.hpp>
#include "wrapper.h"

namespace MetricField {

// Define the MetricField struct
struct MetricFieldStruct {
    MetricField* instance;
};

MetricFieldStruct* create_metric_field(const char* mesh){
    MetricFieldStruct* wrapper = new MetricFieldStruct;
    wrapper->instance = new MetricField(mesh);
    return wrapper;
} 



void call_hello_world(){
    hello_world();
}


// C++ wrapper for extern C call
void call_compute_coeff_improved(MetricFieldStruct* _f, const int _cell_handle_start, const double* _q, const double* _p, double* result){
    MetricField* f = _f->instance;

    Vec3d q = Vec3d(_q);
    Vec3d p = Vec3d(_p);
    OVM::CellHandle cell_handle_start = OVM::CellHandle(_cell_handle_start);

    Vec3d ea;
    Mat3d rotation = f->computeCoeffImproved(cell_handle_start, q, p);

    // maybe this should be 2,1,0, not sure
    ea = rotation.eulerAngles(0,1,2);
    // ea = rotation.eulerAngles(2,1,0);

    result[0] = ea(0);
    result[1] = ea(1);
    result[2] = ea(2);

    return;
}

void new_direction(MetricFieldStruct* _f, const int _cell_handle, const double* _dir, const double* _point, double* result){
    MetricField* f = _f->instance;
    Vec3d dir = Vec3d(_dir);
    Vec3d point = Vec3d(_point);

    Mat3d g = f->metricAtPoint(OVM::CellHandle(_cell_handle), point);
    Vec3d temp = g * dir;

    result[0] = temp(0);
    result[1] = temp(1);
    result[2] = temp(2);
    return;
}


void new_normal(MetricFieldStruct* _f, const int _cell_handle, const double* _normal, const double* _point, double* result){
    MetricField* f = _f->instance;
    Vec3d normal = Vec3d(_normal);
    Vec3d point = Vec3d(_point);

    Mat3d g_inverse = (f->metricAtPoint(OVM::CellHandle(_cell_handle), point)).inverse();

    // std::cout << "g_inverse" << g_inverse << std::endl;

    Vec3d temp = g_inverse * normal;
    result[0] = temp(0);
    result[1] = temp(1);
    result[2] = temp(2);

    // std::cout << "temp " << temp << std::endl;
    return;
}

void transform_to_metric_space(double* _direction, double* _metric){ // transforms the direction with the metric in place
    Vec3d direction = Vec3d(_direction);

    Mat3d metric;
    metric << _metric[0], _metric[1], _metric[2], _metric[1], _metric[3], _metric[4], _metric[2], _metric[4], _metric[5];

    Vec3d temp = metric * direction;
    _direction[0] = temp(0);
    _direction[1] = temp(1);
    _direction[2] = temp(2);
    return;
}


void metric_at_point(MetricFieldStruct* _f, const int _cell_handle, const double* midpoint, double* metric_result /*where to write the metric*/){
    MetricField* f = _f->instance;
    Vec3d point = Vec3d(midpoint);

    Mat3d metric = f->metricAtPoint(OVM::CellHandle(_cell_handle) ,point);
    metric_result[0] = metric.coeff(0,0);
    metric_result[1] = metric.coeff(1,0);
    metric_result[2] = metric.coeff(2,0);
    metric_result[3] = metric.coeff(1,1);
    metric_result[4] = metric.coeff(2,1);
    metric_result[5] = metric.coeff(2,2);
    return;
}


int walk_to_point(MetricFieldStruct* _f, const int _cell_handle_start, const double* _start_position, const double* _end_position){
    MetricField* f = _f->instance;
    auto cell_handle_start = OVM::CellHandle(_cell_handle_start);

    Vec3d start_position = Vec3d(_start_position);
    Vec3d end_position = Vec3d(_end_position);

    // std::cout << "start position " << start_position << std::endl;
    // std::cout << "end position " << end_position << std::endl;
    int res = f->startCellFinder(OVM::CellHandle(_cell_handle_start), start_position, end_position).idx();
    return res; 
    
}

int walk_to_point_initial(MetricFieldStruct* _f, const double* _end_position){
    MetricField* f = _f->instance;
    auto cell_handle_start = OVM::CellHandle(0);

    Vec3d start_position = toVec3d(f->get_tetmesh().barycenter(cell_handle_start));
    Vec3d end_position = Vec3d(_end_position);

    // std::cout << "Walking walking walking " << std::endl;
    const int cell_handle = f->startCellFinder(cell_handle_start, start_position, end_position).idx();

    return cell_handle;
    // return 0;
    
}

} // end namespace MetricField
