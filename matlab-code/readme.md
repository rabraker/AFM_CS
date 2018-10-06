## Overview
These classes represent a (hopefully) extendable basic framework to implement what I am calling CS-Entities.

It is my understanding that the common theme in all the CS modalities is some type of scan, followed by an XY jump, followed by another scan. These scans might be mu-paths, single points, spirals etc. 
I term these scan patterns *Measurement Entities*(ME). For every CS modality, the respective MEs share common characteristics:
1. x&y start coordinates
* an index, ie, ME number 235
* some kind of metric of how "long" we should spend at that MI.

The ME only defines the system behaivior while we are actually taking a measurement. The other crucial aspect is how we move *between*  measurements. I call these Move-Entities (MVE). 

Perhaps the easiest way to see what I'm talking about is to dig into `cs-trajectories.m`. 

The Classes
```
CsEntity.m
```
is an abstract base class which implementations of both MVE and ME should subclass. In both cases, we need a vector of x-coods, y-coords, and a (constant) index vector. How exactly the vectors of x-coords and y-coords are defined is left up to the implementation in the sub-class. The easiest example to understand is probly the pair 
```
MeasEntityStatic.m
MoveEntityStatic.m
```
The class `MeasEntityStatic` is for CS at single points only. `MoveEntityStatic` implements step-input defined moves. The essential difference between these two classes is that `MeasEntityStatic` defines an index, (ie, Measurement Entitity number 235) while `MoveEntitystatic` automatically sets the index to 0. Another idea would be to set the MVE index to -j, if j is the index of the associated ME.

Anyway, `CsEntity` itself implements methods: the first, `as_vector(self)` interleaves xref_vec, yref_vec, and index_vec into a single vector, so they can be stuffed into a single FIFO buffer. 

The second method is `as_matrix(self)` which will return a 3-row by N-col matrix of the data.

A bit of oddness is the `factory` methods in the class implementations, e.g., MeasEntityStatic and MoveEntityStatic. This is a static method which just returns a handle to the class itself. This allows us to parameterize the class over different variables, in this case, the number of sample periods the entity is defined for, N. This same idea is used to parameterize, for example, the *scan rate* in `MeasEntityMuPath`. This seems kind of clunky, but I couldn't figure out a better way... 

The class `MasterTrajster` glues all of this together. It is instantiated with the prototype

```
self = MasterTrajster(XR, YR, MoveEntity, MeasEntity)
```
where
----
XR: a vector of x-coords
YR: a vector of y-coords
MoveEntity: a *handle* to an object which sub-classes CsEntity
MeasEntity: a *handle* to an object which sub-classes CsEntity

An example of all of this put together can be found in the script `me_mve_example.m`.

Note: by parameterizing the ME and MVE classes, the MasterTrajster class is completely ignorant of, for example N_mve. This is important because a better implementation will either have N_mve depend on the move size, or even better, not exist at all. Ideally, we will detect on the FPGA when the Move has settled. This basically means we could use a single value, instead of a vector for each sample instant, for MVE. 


