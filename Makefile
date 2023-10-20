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
	BUILD_CMD = emmake ninja
else
	CMAKE_CMD = cmake
	BUILD_CMD = ninja
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
all: check_dependencies check_submodule build

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
	which yarn >/dev/null || (echo "Installing Yarn via Homebrew..."; brew install yarn); \
	which ninja >/dev/null || (echo "Installing Ninja via Homebrew..."; brew install ninja)
else ifeq ($(OS), Linux)
	sudo apt-get update; \
	MAKEFILE_DEPS="g++ llvm clang emscripten yarn cmake libeigen3-dev libgmp-dev libmpfr-dev googletest libgtest-dev ninja-build"; \
	for DEP in $$MAKEFILE_DEPS; do \
		dpkg -s $$DEP >/dev/null 2>&1 || (echo "Installing $$DEP via package manager..."; sudo apt-get install -y $$DEP); \
	done
else ifeq ($(OS), MINGW64_NT-10.0)
	@echo "Please ensure you have installed LLVM, Emscripten and Yarn manually, and they are available in the PATH."
else
	@echo "Unsupported OS. Please install the packages manually."
endif
	@echo "Dependencies check complete."

.PHONY: check_submodule
check_submodule:
	@if [ ! "$(shell git submodule status | grep MeshCartographyLib | cut -c 1)" = "-" ]; then \
		echo "MeshCartographyLib submodule already initialized and updated."; \
	else \
		echo "MeshCartographyLib submodule is empty. Initializing and updating..."; \
		git submodule update --init -- MeshCartographyLib; \
	fi

.PHONY: update_submodule
update_submodule:
	@echo "Updating MeshCartographyLib submodule..."; \
	git submodule update --remote MeshCartographyLib; \

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
	$(BUILD_CMD) -C $(PROJECT_DIR)/build -j $(shell sysctl -n hw.logicalcpu)
else ifeq ($(OS), Linux)
	$(BUILD_CMD) -C $(PROJECT_DIR)/build -j $(shell nproc)
endif

$(DATA_DIR):
	mkdir -p $(DATA_DIR)
	mkdir -p $(ASSETS_DIR)


# Optional
.PHONY: install_analysis
install_analysis:
	@echo "Installing analysis dependencies..."
	git clone https://github.com/Jan-Piotraschke/2DTissue-Analysis.git $(EXTERNAL_DIR)/2DTissue-Analysis

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
