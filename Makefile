# Makefile for building 2DTissue

SHELL := /bin/bash

# Paths
PROJECT_DIR := $(shell pwd)
EXTERNAL_DIR := $(PROJECT_DIR)/external
DATA_DIR := $(PROJECT_DIR)/data
ASSETS_DIR := $(PROJECT_DIR)/assets
ARCHITECTURE := $(shell uname -m)

# Platform selection
BUILD_DIR = build
CMAKE_CMD = cmake
BUILD_CMD = ninja

# Determine OS
OS := $(shell uname -s)

# Compiler paths
ifeq ($(OS), Darwin)
	C_COMPILER=$(shell xcrun --find clang)
	CXX_COMPILER=$(shell xcrun --find clang++)
	SQLITE3_LIB=$(shell brew --prefix sqlite3)/lib
	SQLITE3_INCLUDE=$(shell brew --prefix sqlite3)/include
	VCPKG_ROOT := $(PROJECT_DIR)/vcpkg
	VCPKG_TOOLCHAIN := $(VCPKG_ROOT)/scripts/buildsystems/vcpkg.cmake
else ifeq ($(OS), Linux)
	C_COMPILER=/usr/bin/gcc
	CXX_COMPILER=/usr/bin/g++
	SQLITE3_LIB=/usr/local/lib
	SQLITE3_INCLUDE=/usr/local/include
	VCPKG_ROOT := $(PROJECT_DIR)/vcpkg
	VCPKG_TOOLCHAIN := $(VCPKG_ROOT)/scripts/buildsystems/vcpkg.cmake
endif

.PHONY: all
all: check_dependencies check_submodule init_vcpkg build

# Check if LLVM and other dependencies are installed
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
	which yarn >/dev/null || (echo "Installing Yarn via Homebrew..."; brew install yarn); \
	which ninja >/dev/null || (echo "Installing Ninja via Homebrew..."; brew install ninja); \
	which sqlite3 >/dev/null || (echo "Installing SQLite3 via Homebrew..."; brew install sqlite3)
else ifeq ($(OS), Linux)
	sudo apt-get update; \
	MAKEFILE_DEPS="g++ llvm clang yarn cmake libeigen3-dev libgmp-dev libmpfr-dev googletest libgtest-dev ninja-build"; \
	for DEP in $$MAKEFILE_DEPS; do \
		dpkg -s $$DEP >/dev/null 2>&1 || (echo "Installing $$DEP via package manager..."; sudo apt-get install -y $$DEP); \
	done
else ifeq ($(OS), MINGW64_NT-10.0)
	@echo "Please ensure you have installed LLVM, Yarn, and SQLite manually, and they are available in the PATH."
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
	git submodule update --remote MeshCartographyLib;

# Initialize and set up vcpkg
.PHONY: init_vcpkg
init_vcpkg:
	@if [ ! "$(shell git submodule status | grep vcpkg | cut -c 1)" = "-" ]; then \
		echo "vcpkg submodule already initialized and updated."; \
	else \
		echo "vcpkg submodule is empty. Initializing and updating..."; \
		git submodule update --init -- vcpkg; \
	fi
	@echo "Initializing vcpkg..."
	cd $(VCPKG_ROOT) && ./bootstrap-vcpkg.sh
	@echo "Integrating vcpkg with system..."
	$(VCPKG_ROOT)/vcpkg integrate install
	@echo "Installing librdkafka via vcpkg..."
	$(VCPKG_ROOT)/vcpkg install librdkafka
	@echo "vcpkg initialization and library installation complete."

.PHONY: build
build: $(DATA_DIR)
	echo "Building for platform: $(PLATFORM)"; \
	$(CMAKE_CMD) -S $(PROJECT_DIR) \
			-B $(PROJECT_DIR)/$(BUILD_DIR) \
			-DCMAKE_BUILD_TYPE=Release \
			-DCMAKE_C_COMPILER=$(C_COMPILER) \
			-DCMAKE_CXX_COMPILER=$(CXX_COMPILER) \
			-DCMAKE_CXX_STANDARD=20 \
			-DCMAKE_OSX_ARCHITECTURES=$(ARCHITECTURE) \
			-DSQLITE3_INCLUDE_DIR=$(SQLITE3_INCLUDE) \
			-DSQLITE3_LIBRARY=$(SQLITE3_LIB) \
			-DCMAKE_TOOLCHAIN_FILE=$(VCPKG_TOOLCHAIN) \
			-GNinja
ifeq ($(OS), Darwin)
	$(BUILD_CMD) -C $(PROJECT_DIR)/$(BUILD_DIR) -j $(shell sysctl -n hw.logicalcpu)
else ifeq ($(OS), Linux)
	$(BUILD_CMD) -C $(PROJECT_DIR)/$(BUILD_DIR) -j $(shell nproc)
endif

$(DATA_DIR):
	mkdir -p $(DATA_DIR)
	mkdir -p $(ASSETS_DIR)

.PHONY: start_test_kafka
start_test_kafka:
	@echo "Starting Zookeeper and Kafka services using Homebrew..."; \
	brew services start zookeeper; \
	brew services start kafka; \
	echo "Waiting for Kafka to start..."; \
	sleep 10; \
	echo "Running Kafka consumer..."; \
	/opt/homebrew/opt/kafka/bin/kafka-console-consumer --bootstrap-server localhost:9092 --topic simulation_topic

.PHONY: stop_test_kafka
stop_test_kafka:
	@echo "Stopping Kafka and Zookeeper services using Homebrew..."; \
	brew services stop kafka; \
	brew services stop zookeeper; \
	echo "Verifying services have stopped..."; \
	STATUS_KAFKA=$$(brew services list | grep kafka | awk '{print $$2}'); \
	STATUS_ZOOKEEPER=$$(brew services list | grep zookeeper | awk '{print $$2}'); \
	if [ "$$STATUS_KAFKA" = "none" ] && [ "$$STATUS_ZOOKEEPER" = "none" ]; then \
		echo "Kafka and Zookeeper have successfully stopped."; \
	else \
		echo "Warning: One or both services are still running."; \
	fi

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
