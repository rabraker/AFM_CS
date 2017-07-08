## Overview
These classes represent a (hopefully) extendable basic framework to implement what I am calling CS-Entities.

It is my understanding that the common theme in all the CS modalities is some type of scan, followed by and XY jump, followed by another scan. As I mentioned before, an ME might be a single pixel measurement point, a mu-scan etc. But in general, we will have N MIs. For every CS modality, the respective MEs chare common characteristics:
1. x&y start coordinates
* an index, ie, ME number 235
* some kind of metric of how "long" we should spend at that MI.

The ME only defines the behaivior while we are acutally taking a measurement. The other crucial aspect is how we move *between*  measurements. I call these Move-Entities (MVE). 


The Class
```
CsEntity.m
```
is an abstract base class for both MVE and ME. In both cases, we need a vector of x-coods, y-coords, and a (constant) index vector. How exactly the vectors of x-coords and y-coords are defined is left up to the implementation in the sub-class. As an example, I have written
```
MeasEntityStatic.m
MoveEntityStatic.m
```
The only difference between these two classes is that `MeasEntityStatic` defines an index, while `MoveEntitystatic` automatically sets the index to 0. Another idea would be to set the MVE index to -j, if j is the index of the associated ME.

Anyway, `CsEntity` implements a single method, which interleaves xref_vec, yref_vec, and index_vec into a single vector, so they can be stuffed into a single FIFO buffer. 

A bit of oddness the `factory` method in both Measentitystatic and MoveEntitystatic. This is a static method which just returns a handle to the class itself, which allows us to parameterize the class over different variables, in this case, the number of sample periods the entity is defined for, N. This same idea will be used to parameterize, for example, the *scan rate* in a futre `MeasEntityMuPath`. This seems kind of clunky, but I couldn't figure out a better way... 

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


