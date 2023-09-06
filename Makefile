# Makefile for building 2DTissue

SHELL := /bin/bash

# Paths
PROJECT_DIR := $(shell pwd)
EXTERNAL_DIR := $(PROJECT_DIR)/external
DATA_DIR := $(PROJECT_DIR)/data
ASSETS_DIR := $(PROJECT_DIR)/assets
ARCHITECTURE := $(shell uname -m)

# Platform selection
PLATFORM ?= executive
ifeq ($(PLATFORM), wasm)
	CMAKE_CMD = emcmake cmake
else
	CMAKE_CMD = cmake
endif

# Determine OS
OS := $(shell uname -s)

# Compiler paths
ifeq ($(OS), Darwin)
	C_COMPILER=$(shell brew --prefix llvm)/bin/clang
	CXX_COMPILER=$(shell brew --prefix llvm)/bin/clang++
else ifeq ($(OS), Linux)
	C_COMPILER=/usr/bin/gcc
	CXX_COMPILER=/usr/bin/g++
endif

.PHONY: all
all: check_dependencies build_cgal build_libroadrunner_deps build_llvm_13 build_libroadrunner build

# Check if LLVM and Emscripten are installed, if not, install using apt-get
.PHONY: check_dependencies
check_dependencies:
	@echo "Checking dependencies..."
ifeq ($(OS), Darwin)
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
	which ninja >/dev/null || (echo "Installing Ninja via Homebrew..."; brew install ninja)
else ifeq ($(OS), Linux)
	sudo apt-get update; \
	MAKEFILE_DEPS="g++ llvm clang emscripten yarn cmake libboost-all-dev libeigen3-dev libgmp-dev libmpfr-dev googletest libgtest-dev libomp-dev libassimp-devninja-build"; \
	for DEP in $$MAKEFILE_DEPS; do \
		dpkg -s $$DEP >/dev/null 2>&1 || (echo "Installing $$DEP via package manager..."; sudo apt-get install -y $$DEP); \
	done
else ifeq ($(OS), MINGW64_NT-10.0)
	@echo "Please ensure you have installed LLVM, Emscripten, Assimp and Yarn manually, and they are available in the PATH."
else
	@echo "Unsupported OS. Please install the packages manually."
endif
	@echo "Dependencies check complete."

# Build and install CGAL
.PHONY: build_cgal
build_cgal:
	@echo "Building CGAL..."
ifeq ($(OS), Linux)
	if [ ! -d "$(EXTERNAL_DIR)/cgal" ]; then \
		git clone -b 'v5.5.2' --single-branch --depth 1 https://github.com/CGAL/cgal.git $(EXTERNAL_DIR)/cgal; \
		mkdir -p build/cgal; \
		cd build/cgal && $(CMAKE_CMD) $(EXTERNAL_DIR)/cgal -DCMAKE_BUILD_TYPE=Release -DCGAL_HEADER_ONLY=OFF && make && sudo make install; \
	fi
endif

.PHONY: build_libroadrunner_deps
build_libroadrunner_deps:
	@echo "Installing 2DTissue-deps from source..."
	if [ ! -d "$(EXTERNAL_DIR)" ]; then \
		git clone --recursive https://github.com/MorphoPhysics/2DTissue-deps.git external; \
	fi;
	cd $(EXTERNAL_DIR)/libroadrunner-deps; \
	mkdir -p build; \
	cd build; \
	$(CMAKE_CMD) -GNinja -DCMAKE_INSTALL_PREFIX="../install-release" \
		-DCMAKE_BUILD_TYPE="Release" \
		-DCMAKE_CXX_STANDARD=17 \
		-DCMAKE_OSX_ARCHITECTURES=$(ARCHITECTURE) ..; \
	ninja; \
	ninja install

.PHONY: build_llvm_13
build_llvm_13:
	@echo "Installing LLVM 13 from source..."; \
	cd $(EXTERNAL_DIR)/libroadrunner-deps/third_party/llvm-13.x; \
	mkdir -p build; \
	cd build; \
	$(CMAKE_CMD) -GNinja -DCMAKE_INSTALL_PREFIX="../install-release" \
		-DCMAKE_BUILD_TYPE="Release" \
		-DCMAKE_CXX_STANDARD=17 \
		-DCMAKE_OSX_ARCHITECTURES=$(ARCHITECTURE) ../llvm; \
	ninja; \
	ninja install

# Build and install libRoadRunner
.PHONY: build_libroadrunner
build_libroadrunner:
	cd $(EXTERNAL_DIR)/roadrunner; \
    mkdir -p build-release; \
    cd build-release; \
    $(CMAKE_CMD) -GNinja -DCMAKE_INSTALL_PREFIX="$(EXTERNAL_DIR)/roadrunner/install-release" \
        -DLLVM_INSTALL_PREFIX="$(EXTERNAL_DIR)/libroadrunner-deps/third_party/llvm-13.x/install-release" \
        -DRR_DEPENDENCIES_INSTALL_PREFIX="$(EXTERNAL_DIR)/libroadrunner-deps/install-release" \
        -DCMAKE_BUILD_TYPE="Release" \
        -DCMAKE_CXX_STANDARD=17 \
        -DCMAKE_OSX_ARCHITECTURES=$(ARCHITECTURE) ..; \
    ninja; \
    ninja install

.PHONY: build
build: $(DATA_DIR)
	echo "Building for platform: $(PLATFORM)"; \
	$(CMAKE_CMD) -S $(PROJECT_DIR) \
			-B $(PROJECT_DIR)/build \
			-DCMAKE_BUILD_TYPE=Release \
			-DCMAKE_C_COMPILER=$(C_COMPILER) \
			-DCMAKE_CXX_COMPILER=$(CXX_COMPILER) \
			-DCMAKE_CXX_STANDARD=20 \
			-DCMAKE_OSX_ARCHITECTURES=$(ARCHITECTURE) \
			-GNinja
ifeq ($(OS), Darwin)
	ninja -C $(PROJECT_DIR)/build -j $(shell sysctl -n hw.logicalcpu)
else ifeq ($(OS), Linux)
	ninja -C $(PROJECT_DIR)/build -j $(shell nproc)
endif

$(DATA_DIR):
	mkdir -p $(DATA_DIR)
	mkdir -p $(ASSETS_DIR)

# Cleaning
.PHONY: clean
clean:
	rm -rf $(PROJECT_DIR)/build $(DATA_DIR) $(ASSETS_DIR)

.PHONY: clean_data
clean_data:
	rm -rf $(DATA_DIR) $(ASSETS_DIR)

.PHONY: distclean
distclean: clean
	rm -rf $(PROJECT_DIR)/build/*
