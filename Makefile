INCL_DIR := $(CURDIR)/include
SRC_DIR := $(CURDIR)/src
LIB_DIR := $(CURDIR)/lib


all:
	# Generate the SWIG wrapper (C++ -> Python)
	swig -c++ -python -o $(SRC_DIR)/swig_wrapper.cpp -outdir $(SRC_DIR) $(SRC_DIR)/swig_wrapper.i

	# Compile the SWIG wrapper using setuptools
	python setup.py build_ext --build-lib $(LIB_DIR)
