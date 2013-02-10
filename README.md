FiberAllocation
=======

Assigns fibers to galaxy targets in the sky. It takes as inputs a set of galaxies and fibers and returns
a match fiber-galaxy. The fibers are arranged in an hexagonal tile.

Languages
=========

Python

In Context
===========

Open Issues to Continue Development
===================================

* The code assumes a flat sky.
* Option for writing to disk the match fiber-galaxy has to be implemented
* Include the priority information of different galaxy populations into the matching routine.
* Write tests for quality assessment.
* Implement Simulated Annealing.

Inputs/Outputs
==============

Required galaxy catalog input:
Is an ASCII file with four columns separated by spaces. Each line
corresponding to a galaxy with: 
* unique galaxy_ID in the sky (scalar, int)
* galaxy type (scalar, int). Identifies three kinds of galaxies: ELG, LRG, QSO.
* ra (scalar, float) position in degrees
* dec (scalar, float) position in degrees

Required inputs to construct the fibers:
* fiber_pitch (scalar, float) default 0.01. Distance between fibers in
degrees. Same number for all fibers. 
* Ndiameter (scalar, integer) default 73. Number of fibers along the
major axis in a hexagonal tile 
* center_x (scalar, float) default 0.0. Position of the hexagon's
center in the ra-direction. In degrees 
* center_y (scalar, float) default 0.0. Position of the hexagon's
center in the dec-direction. In degrees. 


Required input to feed the allocation routine:
* tile_ID (scalar, integer) default 1. A unique ID during the
survey. Identifies the tile that triggered this instance of fiber
allocation. The module will take as input a whole array of tile_ids,
one for each galaxy (or -1).  
* visit_ID (scalar, integer) default 1. A unique ID for each
tile. Identifies the number of times that a particular tile_ID has
been observed. It is the number of times that a tile has been visited  
* patrol_radius (scalar, float) default 0.01. patrol radius of each
fiber. It is the same for all fibers. In degrees. 
* exclusion_radius (scalar, float) default 0.0. Minimum separation
between the centers of two fibers are allowed to have before having a
collision. In degrees. 
* fiber_size (scalar, float) default 0.0. Diameter of the fibers. The
exclusion_radius parameter should be equal or larger than half the
fiber_size. It is the same for all fibers. In degrees. 


Data that is produced inside the routine as initialization:
A list of fibers with the following information:
* A unique fiber ID (scalar, int)
* The unperturbed ra-dec coordinates (scalar, floats). In degrees.

Data that is produced inside the routine after the allocation is completed.
For each fiber:
* The corresponding galaxy_ID that the fiber will observe (scalar, int). 
* The value is -1 if the fiber is free.
* The perturbed positions in ra-dec coordinates (scalar, floats) after
the allocation. They correspond to the galaxy's coordinates allocated
to this fiber. 

Outputs
The output is a list of galaxies with 7 columns:
* unique galaxy_ID in the sky (scalar, int)
* galaxy type (scalar, int). Corresponds to ELG, LRG, QSO.
* ra (scalar float). Same as input.
* dec (scalar, float). Same as input.
* tile_ID (scalar, int). Sames as input if galaxy was allocated. -1 if not.
* visit_ID (scalar, int). Same as input if galaxy was allocated. -1 if not.
* fiber_ID (scalar, int). fiber_ID of the fiber allocated to it. -1 it none.
 

Function Layout
===============

### User-defined functions
* class FiberSet
* class MockGalaxyCatalog
* function close_match_xy
* function reset_fiber_collisions
* function make_fiber_allocation

### Modules required
* numpy
* sys

Code Procedure
==============

Execution
=========
Currently the code can be executed as
```
# Creates a set of fibers
fibers = FiberSet(Ndiameter=73, fiber_pitch=0.1)

# Creates a set of galaxies
gals = MockGalaxyCatalog(Ngalaxies=50000, radius_fov=fibers.radius,priority_levels=3) 

# Runs the fiber allocation
make_fiber_allocation(fibers, gals, tile_visit_ID=340)
```
