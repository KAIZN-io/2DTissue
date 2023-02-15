build:
	cmake -S /Users/jan-piotraschke/git_repos/Confined_active_particles/ -B /Users/jan-piotraschke/git_repos/Confined_active_particles/build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O2" -DCMAKE_PREFIX_PATH=`julia --project=@. -e 'using CxxWrap; CxxWrap.prefix_path() |> print'`
	cmake --build /Users/jan-piotraschke/git_repos/Confined_active_particles/build --config Release -j 1
	# g++ -std=c++17 -shared -fPIC -o $@ $^  \
	# -I /Users/jan-piotraschke/git_repos/Confined_active_particles/libcxxwrap-julia/include \
	# -I /opt/homebrew/Cellar/julia/1.8.5/include/julia \
	# -DJulia_EXECUTABLE=/opt/homebrew/bin/julia ../libcxxwrap-julia\
	# -DBUILD_JULIA=ON


clean:
	rm -f *.o *.dylib
	rm -rf build