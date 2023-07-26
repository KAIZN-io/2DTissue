# Makefile for building C++ libraries

SHELL := /bin/bash

# Path to project directory
PROJECT_DIR := $(shell pwd)
DATA_DIR := $(PROJECT_DIR)/data
ASSETS_DIR := $(PROJECT_DIR)/assets
ARCHITECTURE := arm64

.PHONY: all
all: check_dependencies build_cgal build_libsbml build

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
		which emcc >/dev/null || (echo "Installing Emscripten via Homebrew..."; brew install emscripten); \
		which assimp >/dev/null || (echo "Installing Assimp via Homebrew..."; brew install assimp); \
		which yarn >/dev/null || (echo "Installing Yarn via Homebrew..."; brew install yarn); \
		which ninja >/dev/null || (echo "Installing Ninja via Homebrew..."; brew install ninja); \
		brew --prefix sundials >/dev/null || (echo "Installing SUNDIALS via Homebrew..."; brew install sundials); \
	elif [ "$$OS" == "Linux" ]; then \
		sudo apt-get update; \
		MAKEFILE_DEPS="g++ llvm clang emscripten yarn cmake libboost-all-dev libeigen3-dev libgmp-dev libmpfr-dev googletest libgtest-dev libomp-dev libassimp-dev libsundials-dev ninja-build"; \
		for DEP in $$MAKEFILE_DEPS; do \
			dpkg -s $$DEP >/dev/null 2>&1 || (echo "Installing $$DEP via package manager..."; sudo apt-get install -y $$DEP); \
		done; \
	elif [ "$$OS" == "MINGW64_NT-10.0" ]; then \
		@echo "Please ensure you have installed LLVM, Emscripten, Assimp, SUNDIALS and Yarn manually, and they are available in the PATH."; \
	else \
		@echo "Unsupported OS. Please install LLVM, Emscripten, Assimp, SUNDIALS and Yarn manually."; \
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

# Build and install libSBML
.PHONY: build_libsbml
build_libsbml:
	@echo "Installing libSBML from source..."
	@if [ ! -d "libsbml" ]; then \
		git clone -b 'v5.20.0' --single-branch --depth 1 https://github.com/sbmlteam/libsbml.git; \
		mkdir -p libsbml/build; \
		cd libsbml/build && cmake -G Ninja ../../libsbml -DENABLE_COMP=ON && ninja && sudo ninja install; \
	fi

# Build and install libRoadRunner
.PHONY: build_libroadrunner
build_libroadrunner:
	@echo "Installing libRoadRunner from source..."
	if [ ! -d "roadrunner" ]; then \
		git clone https://github.com/sys-bio/roadrunner.git; \
	fi; \
	cd roadrunner; \
	mkdir -p build-release; \
	cd build-release; \
	cmake -GNinja -DCMAKE_INSTALL_PREFIX="../install-release" \
	    -DLLVM_INSTALL_PREFIX="../../llvm-13.x/install-release" \
	    -DRR_DEPENDENCIES_INSTALL_PREFIX="../../libroadrunner-deps/install-release" \
	    -DCMAKE_BUILD_TYPE="Release" \
	    -DCMAKE_CXX_STANDARD=17 \
		-DCMAKE_OSX_ARCHITECTURES=$(ARCHITECTURE) ..; \
	ninja; \
	ninja install


.PHONY: build_libroadrunner_deps
build_libroadrunner_deps:
	@echo "Installing libRoadRunner dependencies from source..."
	if [ ! -d "libroadrunner-deps" ]; then \
		git clone https://github.com/sys-bio/libroadrunner-deps.git --recurse-submodules; \
	fi; \
	cd libroadrunner-deps; \
	mkdir -p build; \
	cd build; \
	cmake -GNinja -DCMAKE_INSTALL_PREFIX="../install-release" \
		-DCMAKE_BUILD_TYPE="Release" \
		-DCMAKE_CXX_STANDARD=17 \
		-DCMAKE_OSX_ARCHITECTURES=$(ARCHITECTURE) ..; \
	ninja; \
	ninja install


.PHONY: build_llvm_13
build_llvm_13:
	@echo "Installing LLVM 13 from source..."; \
	if [ ! -d "llvm-13.x" ]; then \
		git clone https://github.com/sys-bio/llvm-13.x.git; \
	fi; \
	cd llvm-13.x; \
	mkdir -p build; \
	cd build; \
	cmake -GNinja -DCMAKE_INSTALL_PREFIX="../install-release" \
		-DCMAKE_BUILD_TYPE="Release" \
		-DCMAKE_OSX_ARCHITECTURES=$(ARCHITECTURE) ../llvm; \
	ninja; \
	ninja install

.PHONY: build
build: $(DATA_DIR)
	@OS=$$(uname -s); \
	if [ "$$OS" == "Darwin" ]; then \
		cmake -S $(PROJECT_DIR) \
			-B $(PROJECT_DIR)/build \
			-DCMAKE_BUILD_TYPE=Release \
			-DCMAKE_C_COMPILER=$(shell brew --prefix llvm)/bin/clang \
			-DCMAKE_CXX_COMPILER=$(shell brew --prefix llvm)/bin/clang++ \
			-DCMAKE_CXX_STANDARD=20 \
			-DCMAKE_OSX_ARCHITECTURES=$(ARCHITECTURE) \
			-GNinja; \
		ninja -C $(PROJECT_DIR)/build -j $$(sysctl -n hw.logicalcpu); \
	elif [ "$$OS" == "Linux" ]; then \
		echo "Building for Linux..."; \
		cmake -S $(PROJECT_DIR) \
			-B $(PROJECT_DIR)/build \
			-DCMAKE_BUILD_TYPE=Release \
			-DCMAKE_C_COMPILER=/usr/bin/gcc \
			-DCMAKE_CXX_COMPILER=/usr/bin/g++ \
			-DCMAKE_CXX_STANDARD=20 \
			-GNinja; \
		ninja -C $(PROJECT_DIR)/build -j $$(nproc); \
	fi

$(DATA_DIR):
	mkdir -p $(DATA_DIR)
	mkdir -p $(ASSETS_DIR)

########################################################################################################################
# Setup                                                                                                                #
########################################################################################################################

.PHONY: clean
clean:
	rm -rf $(PROJECT_DIR)/build $(DATA_DIR) $(ASSETS_DIR)

.PHONY: clean_data
clean_data:
	rm -rf $(DATA_DIR) $(ASSETS_DIR)

.PHONY: distclean
distclean: clean
	rm -rf $(PROJECT_DIR)/build/*
