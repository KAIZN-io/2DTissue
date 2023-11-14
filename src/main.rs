use std::env;
use std::path::PathBuf;

fn main() {
    // let step_count = 20;
    // let save_data = false;
    // let particle_innenleben = false;
    // let optimized_monotile_boundary = false;

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
