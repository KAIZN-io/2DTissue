with import <nixpkgs> {};

stdenv.mkDerivation rec {
  pname = "2DTissue";
  version = "0.0";

  src = ./.;

  buildInputs = [ cmake llvm assimp ninja boost eigen gtest cgal_5 gmp mpfr ];

  configurePhase = ''
    cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$out \
      -GNinja \
      .
  '';

  buildPhase = ''
    ninja
  '';
}
