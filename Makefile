# Directories
INCL_DIR := $(CURDIR)/include
SRC_DIR := $(CURDIR)/src
BUILD_DIR := $(CURDIR)/build
LIB_DIR := $(CURDIR)/lib

# Version header
VERSION := $(shell git describe --tags --always)
VERSION_HEADER := $(INCL_DIR)/version.h
.PHONY: $(VERSION_HEADER) clean all debug]
$(VERSION_HEADER):
	@echo "Updating version header with: $(VERSION)"
	@mkdir -p $(INCL_DIR)
	@echo "#pragma once" > $@
	@echo "#define VERSION \"$(VERSION)\"" >> $@

# Conda environment directories
CONDA_PREFIX := $(shell echo $$CONDA_PREFIX)
CONDA_INCL_DIR := $(CONDA_PREFIX)/include
CONDA_LIB_DIR := $(CONDA_PREFIX)/lib

# Compiler and Flags
CXX := g++
CXXFLAGS := -std=c++17 -g -I$(INCL_DIR) -I$(CONDA_INCL_DIR) -Wall -Wextra -pedantic

# Linker Flags
# Ensure that the library paths are set correctly for linking
LDFLAGS := -L$(LIB_DIR) -L$(CONDA_LIB_DIR) -Wl,-rpath=$(CONDA_LIB_DIR)  # Add rpath for shared libraries
LDLIBS := -lhts  # Link with libhts.a or libhts.so

# Sources and Output
SOURCES := $(filter-out $(SRC_DIR)/swig_wrapper.cpp, $(wildcard $(SRC_DIR)/*.cpp))  # Filter out the SWIG wrapper from the sources
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
TARGET := $(BUILD_DIR)/contextsv

# Default target
all: $(VERSION_HEADER) $(TARGET)

# Debug target
debug: CXXFLAGS += -DDEBUG
debug: all

# Link the executable
$(TARGET): $(OBJECTS) $(VERSION_HEADER)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(VERSION_HEADER)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean the build directory
clean:
	rm -rf $(BUILD_DIR)
	