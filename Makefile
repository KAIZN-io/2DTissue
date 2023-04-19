# Makefile for building C++ libraries with CxxWrap and Julia

SHELL := /bin/bash

# Path to project directory
PROJECT_DIR := $(shell pwd)

# Path to external directory
EXTERNAL_DIR := $(PROJECT_DIR)/external

# Get JlCxx path
JLCXX_PATH := $(shell julia --project=@. -e 'using Pkg; Pkg.instantiate(); using CxxWrap; println(CxxWrap.prefix_path())')

# CMake flags
CMAKE_FLAGS := -DCMAKE_BUILD_TYPE=Release \
               -DCMAKE_CXX_FLAGS="-O2" \
			   -DCMAKE_C_COMPILER=$(shell brew --prefix llvm)/bin/clang \
			   -DCMAKE_CXX_COMPILER=$(shell brew --prefix llvm)/bin/clang++ \
               -DCMAKE_PREFIX_PATH=$(JLCXX_PATH);

# Extern repository name
REPOSITORY := libcxxwrap-julia

.PHONY: all
all: check_dependencies build

# Check if LLVM and Emscripten are installed, if not, install using Homebrew
.PHONY: check_dependencies
check_dependencies:
	@echo "Checking dependencies..."
	# ... (rest of the check_dependencies target)

########################################################################################################################
# Dependencies                                                                                                         #
########################################################################################################################

.PHONY: init
init: clone
	mkdir -p $(EXTERNAL_DIR) && \
	cd $(EXTERNAL_DIR)/$(REPOSITORY)-build && \
	cmake -DJulia_EXECUTABLE=$(shell which julia) ../$(REPOSITORY)
	$(MAKE) -C $(EXTERNAL_DIR)/$(REPOSITORY)-build -j $(shell sysctl -n hw.ncpu)

.PHONY: clone
clone:
ifeq ($(wildcard $(EXTERNAL_DIR)/$(REPOSITORY)/.),)
	mkdir -p $(EXTERNAL_DIR) && \
	git clone https://github.com/JuliaInterop/$(REPOSITORY).git $(EXTERNAL_DIR)/$(REPOSITORY)
	mkdir $(EXTERNAL_DIR)/$(REPOSITORY)-build
else
	@echo "Repository already exists"
endif

.PHONY: build
build:
	cmake -S $(PROJECT_DIR) \
		  -B $(PROJECT_DIR)/build \
		  $(CMAKE_FLAGS)
	$(MAKE) -C $(PROJECT_DIR)/build -j $(shell sysctl -n hw.ncpu)
	@echo "Build finished. The binaries are in $(PROJECT_DIR)/build/lib"


########################################################################################################################
# Setup                                                                                                                #
########################################################################################################################

.PHONY: clean
clean:
	rm -rf $(REPOSITORY)
	rm -rf $(PROJECT_DIR)/build
	rm -rf $(REPOSITORY)-Build

.PHONY: clean_env
clean_env:
	julia --project=@. -e 'using Pkg; for pkg in keys(Pkg.project().dependencies) Pkg.rm(pkg) end'

.PHONY: distclean
distclean: clean
	rm -rf $(PROJECT_DIR)/build/*
