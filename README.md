# Aim of *2DTissue*

This physics engine *2DTissue* was developed to simulate the collective behaviour of active (self-propelled) particles on arbitrary curved closed surfaces.

## Init the Project

run `make init` to initialize the project. With this the required packages for C++ and Julia will be installed.

## Compile the C++ code

1. register your C++ script in the CMakeLists.txt file
2. open a terminal in the root of the project
3. run `make build` to compile the C++ code

## Requirements

For C++ please install following libraries:
1. CGAL (version 5.5.1)
2. Boost (version 1.80.0)
3. Eigen (version 3.4.0_1)

## Theoretical Model

The model described is a Vicsek type model (Vicsek et al. 1995, Physical review letters 75(6): 1226) of spherical active particles with a fixed radius confined to the surface of an ellipsoid. Particle interactions are modelled through forces between neighbouring particles that tend to align their velocities (adapted from Szabo et al. 2006, Physical Review E 74(6): 061908).

### Movie: Vicsek model on a Ellipsoid

https://github.com/MorphoPhysics/2DTissue/assets/78916218/21d177b0-b150-459c-92f2-8b50eb44fa2f
