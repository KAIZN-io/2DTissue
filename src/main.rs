//! # Main Module for the 2DTissue Project
//!
//! This module is the entry point for the 2DTissue project.
//!
//! ## Metadata
//!
//! - **Author:** Jan-Piotraschke
//! - **Date:** 2023-Nov-14
//! - **License:** [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
//!
//! ## Current Status
//!
//! - **Bugs:** None known at this time.
//! - **Todo:** Further development tasks to be determined.


use std::env;
use std::path::PathBuf;
use autocxx::prelude::*;

include_cpp! {
    #include "input.h"
    safety!(unsafe_ffi)
    generate!("do_math")
}

/// Main function
fn main() {
    // let step_count = 20;
    // let save_data = false;
    // let particle_innenleben = false;
    // let optimized_monotile_boundary = false;
    println!("{}", ffi::do_math(12, 13));

    let mesh_cartography_lib_dir_str = env::var("MeshCartographyLib_DIR").expect("MeshCartographyLib_DIR not set");
    let mesh_cartography_lib_dir = PathBuf::from(mesh_cartography_lib_dir_str);
    let new_path = mesh_cartography_lib_dir.join("meshes/ellipsoid_x4.off");
    println!("Parent dir: {:?}", new_path);

    for particle_count in (1000..=1000).step_by(100) {
        println!("Particle count: {}", particle_count);
        let start = std::time::Instant::now();
        let duration = start.elapsed();
        println!("Time taken: {:?} seconds", duration.as_secs_f64());
    }
}
