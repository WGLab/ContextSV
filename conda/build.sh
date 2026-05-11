#!/bin/bash

set -e

echo "Building ContextSV..."
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PREFIX}/lib
export CONDA_PREFIX=$PREFIX
export CXXFLAGS="-I$PREFIX/include $CXXFLAGS"
export LDFLAGS="-L$PREFIX/lib $LDFLAGS"

echo "Checking for HTSLib..."
ls -la $PREFIX/include/htslib/ || echo "HTSLib headers not found"
pkg-config --exists htslib && echo "✓ HTSLib found" || echo "⚠ HTSLib not via pkg-config"

echo "Compiling ContextSV..."
make

echo "Installing ContextSV..."
mkdir -p ${PREFIX}/bin
cp build/contextsv ${PREFIX}/bin/
chmod +x ${PREFIX}/bin/contextsv
cp python/cnv_plots_json.py ${PREFIX}/bin/contextsv-cnv-plot
chmod +x ${PREFIX}/bin/contextsv-cnv-plot

echo "Verifying ContextSV installation..."
$PREFIX/bin/contextsv --help
$PREFIX/bin/contextsv --version
