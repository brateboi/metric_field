# Things/Challenges/Experiments
- Query points (i.e. node centers) lie outside the mesh
    - Need to have mesh covering the target and enclose it --> How to do this automatically?
- How to specify the metric at the vertices --> OpenFlipper Plugin
    - Specify metric at select vertices and interpolate
    


- Is it even useful. Metric experiments. What could be interesting?



## Singularity Graph looks weird even without metric?
### Examples
 - s04c_tetrahedron_depth_9_sg.obj
 - n09c_pyramid_depth_9_sg.obj

## Timings
- Extra effort worth it? Slowdown is about factor 5

|Mesh      | Simone  | With metric |factor|
|----------|---------|-------------|--|
|n09c_pyramid_depth_9   |2s 176ms 912µs        |10s 357ms 409µs | x4.756 |
|s04c_tetrahedron_depth_9|1s 393ms 865µs | 6s 994ms 859µs | x5.018|