# Makefile for building C++ libraries with CxxWrap and Julia

SHELL := /bin/bash

# Path to project directory
PROJECT_DIR := $(HOME)/git_repos/Confined_active_particles

# CMake flags
CMAKE_FLAGS := -DCMAKE_BUILD_TYPE=Release \
               -DCMAKE_CXX_FLAGS="-O2" \
               -DCMAKE_PREFIX_PATH=`julia --project=@. -e 'using CxxWrap; CxxWrap.prefix_path() |> print'`

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
		  -B $(PROJECT_DIR)/src/build \
		  $(CMAKE_FLAGS)
	$(MAKE) -C $(PROJECT_DIR)/src/build -j $(shell nproc)
	@echo "Build finished. The binaries are in $(PROJECT_DIR)/src/build/lib"


########################################################################################################################
# Setup                                                                                                                #
########################################################################################################################

.PHONY: clean
clean:
	rm -rf $(REPOSITORY)
	rm -rf $(PROJECT_DIR)/src/build
	rm -rf $(REPOSITORY)-build

.PHONY: distclean
distclean: clean
	rm -rf $(PROJECT_DIR)/src/build/*
