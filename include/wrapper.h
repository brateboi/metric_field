#ifdef __cplusplus
extern "C" {
#endif

typedef struct MetricFieldStruct MetricFieldStruct;


void call_hello_world();
MetricFieldStruct* create_metric_field(const char* mesh);
void call_compute_coeff_improved(MetricFieldStruct* _f, const int _cell_handle_start, const double* _q, const double* _p, double* result);

void new_direction(MetricFieldStruct* _f, const int _cell_handle, const double* _dir, const double* _point, double* result);
void new_normal(MetricFieldStruct* _f, const int _cell_handle, const double* _normal, const double* _point, double* result);


void transform_to_metric_space(double* direction, double* metric);

void metric_at_point(MetricFieldStruct* _f, const int _cell_handle, const double* midpoint, double* metric_result);

// walks to end_position
int walk_to_point(MetricFieldStruct* _f, const int _cell_handle_start, const double* _start_position, const double* _end_position);

// walks to some point without prior information where to start. Starts at center of CellHandle(0)
int walk_to_point_initial(MetricFieldStruct* _f, const double* _end_position);

#ifdef __cplusplus
}
#endif
