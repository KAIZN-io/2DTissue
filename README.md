# Aim of *2DTissue*

This physics engine *2DTissue* was developed to simulate the collective behaviour of active (self-propelled) particles on arbitrary curved closed surfaces.

## How to compile the C++ code

Make sure you have `nix` installed.

Run `nix-build`. The built binary can be found in `./result/bin/main`.

## Theoretical Model

The model described is a Vicsek type model (Vicsek et al. 1995, Physical review letters 75(6): 1226) of spherical active particles with a fixed radius confined to the surface of an ellipsoid. Particle interactions are modelled through forces between neighbouring particles that tend to align their velocities (adapted from Szabo et al. 2006, Physical Review E 74(6): 061908).

### Movie: Vicsek model on a Ellipsoid

https://github.com/MorphoPhysics/2DTissue/assets/78916218/4867e8fb-8a03-42c8-9770-d94dc20f6894
