#include <metric_field.hpp>
#include "wrapper.h"

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


// put this into metric field file, cant be here to compile to C
Mat3d retrieve_metric(MetricFieldStruct* _f, const uint64_t node_code, const double* _point){
    MetricField* f = _f->instance;
    
    if (f->metricHash.contains(node_code)) {
        return f->metricHash[node_code];
    } else {
        Mat3d metric;
        Vec3d point = Vec3d(_point);
        metric = f->metricAtPoint(point);
        f->metricHash[node_code] = metric;

        return metric;
    }
}


// C++ wrapper for extern C call
void call_compute_coeff(MetricFieldStruct* _f, const double* _q, const double* _p, double* result){
    MetricField* f = _f->instance;

    Vec3d q = Vec3d(_q);
    Vec3d p = Vec3d(_p);

    Vec3d ea;
    Mat3d rotation = f->computeCoeff(q, p);

    // maybe this should be 2,1,0, not sure
    ea = rotation.eulerAngles(0,1,2);
    // ea = rotation.eulerAngles(2,1,0);

    result[0] = ea(0);
    result[1] = ea(1);
    result[2] = ea(2);

    return;
}

void new_direction(MetricFieldStruct* _f, const double* _dir, const double* _point, double* result){
    MetricField* f = _f->instance;
    Vec3d dir = Vec3d(_dir);
    Vec3d point = Vec3d(_point);

    // FIXME TODO Metric Field. find correct CellHandle, CACHE IT
    Mat3d g = f->metricAtPoint(point);
    Vec3d temp = g * dir;

    result[0] = temp(0);
    result[1] = temp(1);
    result[2] = temp(2);
    return;
}


void new_normal(MetricFieldStruct* _f, const double* _normal, const double* _point, double* result){
    MetricField* f = _f->instance;
    Vec3d normal = Vec3d(_normal);
    Vec3d point = Vec3d(_point);

    // FIXME TODO Metric Field. find correct CellHandle, CACHE IT
    Mat3d g_inverse = (f->metricAtPoint(point)).inverse();

    Vec3d temp = g_inverse * normal;
    result[0] = temp(0);
    result[1] = temp(1);
    result[2] = temp(2);

    return;
}

// double* transform_to_metric_space(MetricFieldStruct* _f, const double* start_tangent, const uint64_t node_code, const double* node_midpoint){
//     // transforms the given vector to metric space with g^{1/2}*vec
    
//     Vec3d direction = Vec3d(start_tangent);

//     Mat3d metric = retrieve_metric(_f, node_code, node_midpoint);

//     Vec3d res = metric * direction;
//     return res.data();

// }

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


void metric_at_point(MetricFieldStruct* _f, const double* midpoint, double* metric_result /*where to write the metric*/){
    MetricField* f = _f->instance;
    Vec3d point = Vec3d(midpoint);

    Mat3d metric = f->metricAtPoint(point);
    metric_result[0] = metric.coeff(0,0);
    metric_result[1] = metric.coeff(1,0);
    metric_result[2] = metric.coeff(2,0);
    metric_result[3] = metric.coeff(1,1);
    metric_result[4] = metric.coeff(2,1);
    metric_result[5] = metric.coeff(2,2);
    return;
}

// double vec3_abs_dot_metric(MetricFieldStruct* _f, const double* _dir1, const double* _dir2, const uint64_t node_code, const double* node_midpoint){
//     return(abs(vec3_dot_metric(_f, _dir1, _dir2, node_code, node_midpoint)));
// }

// double vec3_dot_metric(MetricFieldStruct* _f, const double* _dir1, const double* _dir2, const uint64_t node_code, const double* node_midpoint){
//     Vec3d dir1 = Vec3d(_dir1);
//     Vec3d dir2 = Vec3d(_dir2);

//     Mat3d metric = retrieve_metric(_f, node_code, node_midpoint);

//     return(dir1.transpose()*metric*metric*dir2);
// }



