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

Inputs/Outputs
==============

The required inputs to construct the fibers are: 
* fiber_pitch (scalar, float) default 0.01. fiber pitch in degrees
* Ndiameter (scalar, integer) default 73. Number of fibers along the major axis in the hexagonal tile
* center_x (scalar, float) default 0.0. Position of the central fiber in the x-direction
* center_y (scalar, float) default 0.0. Position of the central fiber in the y-direction

The required inputs to construct a mock catalog of galaxies are:
* The radius of the area covered by the galaxies.
* Ngalaxies (scalar, integer) the number of galaxies in the catalog.
* The 
