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
* Option for reading the galaxy positions from an external file has to be implemented.
* Option for writing to disk the match fiber-galaxy has to be implemented
* Include the priority information of different galaxy populations into the matching routine.

Inputs/Outputs
==============

Required inputs to construct the fibers are: 
* fiber_pitch (scalar, float) default 0.01. fiber pitch in degrees
* Ndiameter (scalar, integer) default 73. Number of fibers along the major axis in the hexagonal tile
* center_x (scalar, float) default 0.0. Position of the central fiber in the x-direction
* center_y (scalar, float) default 0.0. Position of the central fiber in the y-direction

Required inputs to construct a mock catalog of galaxies are:
* radius_fov (scatlar, float) default 3.0. The radius in degrees of the area covered by the galaxies.
* Ngalaxies (scalar, integer) default 1000. The number of galaxies in the catalog.
* priority_levels (scalar, integer) default 3. The number of different priority levels from 0 to (priority_levels-1)

The matching routine takes as inputs 
* The Fibers as constructed in the class FiberSet
* The Galaxies as constructed in the class MockGalaxyCatalog
* epsilon (scalar, float) default 0.1. The radius around which the fibers look for galaxies. This should be set to the fiber_pitch value.
* 
* 
