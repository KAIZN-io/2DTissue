.PHONY: all clean build

all: build

build:
	cmake -S /Users/jan-piotraschke/git_repos/Confined_active_particles -B /Users/jan-piotraschke/git_repos/Confined_active_particles/build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O2" -DCMAKE_PREFIX_PATH="/Users/jan-piotraschke/.julia/artifacts/434257acf70ce7407e3946f2c8452e814d7cca64"
	cmake --build /Users/jan-piotraschke/git_repos/Confined_active_particles/build --config Release -j 1
	@echo "Build finished. The binaries are in build/lib"

clean:
	rm -rf $(REPOSITORY)
	rm -rf build