# Polar active particles on arbitrary curved surfaces

This model was developed to explore the effect of surface curvature on the collective behaviour of active (self-propelled) particles.

It is part of a wider research initiative in the context of pattern formation in biological systems.

---

## Init the Project

run `make init` to initialize the project. With this the required packages for C++ and Julia will be installed.

---


## Compile the C++ code for accessing it in Julia

1. register your C++ script in the CMakeLists.txt file
2. open a terminal in the root of the project
3. run `make build` to compile the C++ code


---

## Package management in Julia

The `Project.toml` and `Manifest.toml` contain the package definitions for Julia.  
This allwos us to manage packages with Julias built-in package manager.  
This is how you manage packages with it:

- Open a terminal in the root of the project
- Run `julia -t 8` for parallel computing -> use at least 8 cores
- type `]` (closing square bracket)
- run `activate .`
- use `add <PackageName>` to add a new package
- use `rm <PackageName>` to remove a package
- use `up <PackageName>` to update a package to a newer version

All your modifications to the packages will be reflected in the `Project.toml` and `Manifest.toml` respectively.

---

## Run the simulation in Julia

Execute the following command in the project folder inside your shell: `julia -t 8 --project main.jl`


---

## Trouble shooting:  
1. If you can't use CxxWrap
- please delete all the artifacts inside `.julia/artifacts.`
- Print out the Prefix_path using `CxxWrap.prefix_path()` inside the Julia REPL.
- If this works, you can run `make build` inside your shell.
2. If you are on a Mac and get the warning: `dylib (/usr/local/Cellar/julia/{Version}/lib/libjulia.{Version}.dylib) was built for newer macOS version (X.Y) than being linked (X.Z)`
- reinstall xcode-select by running: `sudo rm -rf /Library/Developer/CommandLineTools && xcode-select --install`

---

## Requirements 

For C++ please install following libraries:
1. CGAL (version 5.5.1)
2. Boost (version 1.80.0)
3. Eigen (version 3.4.0_1)

---

## Theoretical Model

The model described is a Vicsek type model (Vicsek et al. 1995, Physical review letters 75(6): 1226) of spherical active particles with a fixed radius confined to the surface of an ellipsoid. Particle interactions are modelled through forces between neighbouring particles that tend to align their velocities (adapted from Szabo et al. 2006, Physical Review E 74(6): 061908).
