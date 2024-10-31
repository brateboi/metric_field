#ifdef __cplusplus
extern "C" {
#endif

typedef struct MetricFieldStruct MetricFieldStruct;


void call_hello_world();
// double* call_compute_coeff(MetricField* _f, const double* _q, const double* _p);
MetricFieldStruct* create_metric_field(const char* mesh);
void call_compute_coeff(MetricFieldStruct* _f, const double* _q, const double* _p, double* result);

void new_direction(MetricFieldStruct* _f, const double* _dir, const double* _point, double* result);
void new_normal(MetricFieldStruct* _f, const double* _normal, const double* _point, double* result);

// double vec3_abs_dot_metric(MetricFieldStruct* _f, const double* dir1, const double* dir2, const uint64_t node_code, const double* node_midpoint);
// double vec3_dot_metric(MetricFieldStruct* _f, const double* _dir1, const double* _dir2, const uint64_t node_code, const double* node_midpoint);

void transform_to_metric_space(double* direction, double* metric);

void metric_at_point(MetricFieldStruct* _f, const double* midpoint, double* metric_result); 

#ifdef __cplusplus
}
#endif
