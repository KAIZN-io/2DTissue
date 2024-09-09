# Aim of *2DTissue*

This physics engine *2DTissue* was developed to simulate the collective behaviour of active (self-propelled) particles on arbitrary curved closed surfaces.

## How to compile the C++ code

Simply run `make` in your terminal in the root of this project and everything will get installed and compiled.

## How to run the simulation

```bash
./build/main --step-count 50 --particle-count 100 --step-time 0.02
```

or if you need help:

```bash
./build/main --help
```

## Format the code

CppCheck will help you catch errors and bugs in your C++ code through static analysis.

Navigate into the src folder and run the following command:

```bash
cppcheck --enable=all --inconclusive --force --suppress=missingIncludeSystem ./simulation
```

and

```bash
find . -type f \( -name "*.cpp" -o -name "*.h" \) ! -name "argparse.hpp" -exec clang-format -i {} \;
```

## Theoretical Model

The model described is a Vicsek type model (Vicsek et al. 1995, Physical review letters 75(6): 1226) of spherical active particles with a fixed radius confined to the surface of an ellipsoid. Particle interactions are modelled through forces between neighbouring particles that tend to align their velocities (adapted from Szabo et al. 2006, Physical Review E 74(6): 061908).

### Movie: Vicsek model on a Ellipsoid

https://github.com/MorphoPhysics/2DTissue/assets/78916218/54e727c0-85f1-4cc0-9e52-5fbebdf1470c
