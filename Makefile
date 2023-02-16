# Important for the development process
# ============================================================
# SHELL := /bin/bash


REPOSITORY=libcxxwrap-julia

.PHONY: all
all: build


########################################################################################################################
# DEPENDENCIES                                                                                                         #
########################################################################################################################

# Init Cxxwrap
.PHONY: init
init: clone
	cd $(REPOSITORY)-build && \
	cmake -DJulia_EXECUTABLE=/opt/homebrew/bin/julia  ../$(REPOSITORY)
	cmake --build $(REPOSITORY)-build --config Release

# Clone the repository
.PHONY: clone
clone:
ifeq ($(wildcard ./$(REPOSITORY)/.),)
	git clone https://github.com/JuliaInterop/$(REPOSITORY).git
	mkdir $(REPOSITORY)-build
else
	@echo "already exists"
endif

# Build your own C++ library
.PHONY: build
build:
	# -S: path to the source repository which contains CMakeLists.txt; -B: path to the build directory
	cmake -S ${HOME}/git_repos/Confined_active_particles \
		  -B ${HOME}/git_repos/Confined_active_particles/build \
		  -DCMAKE_BUILD_TYPE=Release \
		  -DCMAKE_CXX_FLAGS="-O2" \
		  -DCMAKE_PREFIX_PATH=`julia --project=@. -e 'using CxxWrap; CxxWrap.prefix_path() |> print'`
	cmake --build ${HOME}/git_repos/Confined_active_particles/build --config Release -j 1
	@echo "Build finished. The binaries are in build/lib"


########################################################################################################################
# SETUP                                                                                                                #
########################################################################################################################

# Remove everything
.PHONY: clean
clean:
	rm -rf $(REPOSITORY)
	rm -rf build
	rm -rf $(REPOSITORY)-build