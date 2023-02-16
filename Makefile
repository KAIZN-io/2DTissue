.PHONY: all clean build

all: build

build:
	cmake -S ${HOME}/git_repos/Confined_active_particles -B ${HOME}/git_repos/Confined_active_particles/build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O2" -DCMAKE_PREFIX_PATH=`julia --project=@. -e 'using CxxWrap; CxxWrap.prefix_path() |> print'`
	cmake --build ${HOME}/git_repos/Confined_active_particles/build --config Release -j 1
	@echo "Build finished. The binaries are in build/lib"

clean:
	rm -rf $(REPOSITORY)
	rm -rf build