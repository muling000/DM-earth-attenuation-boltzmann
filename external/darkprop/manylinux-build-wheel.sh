#!/bin/bash
#
# This script runs in docker container manylinux2014_x86_64.
# At the root directory of the source code, the container can be started as
#
# $ dokcer run -it -v `pwd`:/io quay.io/pypa/manylinux_2_28_x86_64
#
# Then in the container, cd to /io and run this script
#
# # cd /io
# # ./manylinux-build-wheel.sh
#
# The wheels will be generated in the wheelhouse folder.

dnf update -y
dnf install -y epel-release
dnf install -y gsl-devel hdf5-devel spdlog-devel

OUTDIR=dist

mkdir ${OUTDIR}

for PYVERSION in "cp38-cp38" "cp39-cp39" "cp310-cp310" "cp311-cp311" "cp312-cp312"; do
    PYTHON="/opt/python/${PYVERSION}/bin/python3"
    $PYTHON -m build --wheel --outdir ${OUTDIR}
done

auditwheel repair ${OUTDIR}/*.whl

exit
