libhello.dylib: libhello.cpp
	g++ -std=c++14 -shared -fPIC -o $@ $^  \
	-I /Users/jan-piotraschke/git_repos/Confined_active_particles/libcxxwrap-julia/include \
	-I /opt/homebrew/Cellar/julia/1.8.5/include/julia \
	-DJulia_EXECUTABLE=/opt/homebrew/bin/julia \
	-DBUILD_JULIA=ON


clean:
	rm -f *.o *.dylib