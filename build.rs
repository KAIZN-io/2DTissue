fn main() -> miette::Result<()>  {
    let include_path = std::path::PathBuf::from("src");
    // let pmp_path = std::path::PathBuf::from("MeshCartographyLib/pmp-library/src");

    let mut builder = autocxx_build::Builder::new("src/main.rs", &[&include_path])
        .build()?;
    builder.flag_if_supported("-std=c++17").compile("autocxx-demo");

    println!("cargo:rerun-if-changed=src/main.rs");

    Ok(())
}
