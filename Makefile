INCL_DIR := $(CURDIR)/include
SRC_DIR := $(CURDIR)/src
LIB_DIR := $(CURDIR)/lib


all:
	# Generate the SWIG wrapper (C++ -> Python)
	swig -c++ -python -I$(INCL_DIR) -o $(SRC_DIR)/swig_wrapper.cpp -outdir $(LIB_DIR) $(SRC_DIR)/swig_wrapper.i

	# Compile the SWIG wrapper using setuptools
	python3 setup.py build_ext --build-lib $(LIB_DIR)
