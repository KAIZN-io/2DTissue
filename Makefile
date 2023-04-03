# Makefile for building C++ libraries with CxxWrap and Julia

SHELL := /bin/bash

# Path to project directory
PROJECT_DIR := $(shell pwd)

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
all: build


########################################################################################################################
# Dependencies                                                                                                         #
########################################################################################################################

.PHONY: init
init: clone
	cd $(REPOSITORY)-build && \
	cmake -DJulia_EXECUTABLE=$(shell which julia) ../$(REPOSITORY)
	$(MAKE) -C $(REPOSITORY)-build -j $(shell nproc)

.PHONY: clone
clone:
ifeq ($(wildcard ./$(REPOSITORY)/.),)
	git clone https://github.com/JuliaInterop/$(REPOSITORY).git
	mkdir $(REPOSITORY)-build
else
	@echo "Repository already exists"
endif

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
	rm -rf $(REPOSITORY)
	rm -rf $(PROJECT_DIR)/build
	rm -rf $(REPOSITORY)-Build

.PHONY: clean_env
clean_env:
	julia --project=@. -e 'using Pkg; for pkg in keys(Pkg.project().dependencies) Pkg.rm(pkg) end'

.PHONY: distclean
distclean: clean
	rm -rf $(PROJECT_DIR)/build/*
