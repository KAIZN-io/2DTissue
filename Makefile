# Makefile for building C++ libraries

SHELL := /bin/bash

# Path to project directory
PROJECT_DIR := $(shell pwd)

# CMake flags
CMAKE_FLAGS := -DCMAKE_BUILD_TYPE=Release \
               -DCMAKE_CXX_FLAGS="-O2" \
				-DCMAKE_C_COMPILER=$(shell brew --prefix llvm)/bin/clang \
			   -DCMAKE_CXX_COMPILER=$(shell brew --prefix llvm)/bin/clang++;

.PHONY: all
all: check_dependencies build_cgal build

# Check if LLVM and Emscripten are installed, if not, install using apt-get
.PHONY: check_dependencies
check_dependencies:
	@echo "Checking dependencies..."
	@OS=$$(uname -s); \
	if [ "$$OS" == "Darwin" ]; then \
		LLVM_PATH=$$(brew --prefix llvm); \
		if [ -z "$$LLVM_PATH" ]; then \
			echo "Installing LLVM via Homebrew..."; \
			brew install llvm; \
			LLVM_PATH=$$(brew --prefix llvm); \
		fi; \
		export PATH="$$LLVM_PATH/bin:$$PATH"; \
		export LDFLAGS="-L$$LLVM_PATH/lib $$LDFLAGS"; \
		export CPPFLAGS="-I$$LLVM_PATH/include $$CPPFLAGS"; \
		CMAKE_FLAGS="$$CMAKE_FLAGS -DCMAKE_C_COMPILER=$$LLVM_PATH/bin/clang -DCMAKE_CXX_COMPILER=$$LLVM_PATH/bin/clang++"; \
		which emcc >/dev/null || (echo "Installing Emscripten via Homebrew..."; brew install emscripten); \
		which assimp >/dev/null || (echo "Installing Assimp via Homebrew..."; brew install assimp); \
		which yarn >/dev/null || (echo "Installing Yarn via Homebrew..."; brew install yarn); \
	elif [ "$$OS" == "Linux" ]; then \
		sudo apt-get update; \
		MAKEFILE_DEPS="g++ llvm clang emscripten yarn cmake libboost-all-dev libeigen3-dev libgmp-dev libmpfr-dev googletest libgtest-dev libomp-dev libassimp-dev"; \
		for DEP in $$MAKEFILE_DEPS; do \
			dpkg -s $$DEP >/dev/null 2>&1 || (echo "Installing $$DEP via package manager..."; sudo apt-get install -y $$DEP); \
		done; \
	elif [ "$$OS" == "MINGW64_NT-10.0" ]; then \
		@echo "Please ensure you have installed LLVM, Emscripten, Assimp, and Yarn manually, and they are available in the PATH."; \
	else \
		@echo "Unsupported OS. Please install LLVM, Emscripten, Assimp, and Yarn manually."; \
	fi
	@echo "Dependencies check complete."

# Build and install CGAL
.PHONY: build_cgal
build_cgal:
	@OS=$$(uname -s); \
	if [ "$$OS" == "Linux" ]; then \
		if [ ! -d "cgal" ]; then \
			git clone -b 'v5.5.2' --single-branch --depth 1 https://github.com/CGAL/cgal.git; \
			mkdir -p build/cgal; \
			cd build/cgal && cmake ../../cgal -DCMAKE_BUILD_TYPE=Release -DCGAL_HEADER_ONLY=OFF && make && sudo make install; \
		fi; \
	fi

.PHONY: build
build:
	cmake -S $(PROJECT_DIR) \
		  -B $(PROJECT_DIR)/build \
		  $(CMAKE_FLAGS)
	$(MAKE) -C $(PROJECT_DIR)/build -j $(shell nproc)
	@echo "Build finished. The binaries are in $(PROJECT_DIR)/build/lib"

########################################################################################################################
# Setup                                                                                                                #
########################################################################################################################

.PHONY: clean
clean:
	rm -rf $(PROJECT_DIR)/build

.PHONY: distclean
distclean: clean
	rm -rf $(PROJECT_DIR)/build/*
