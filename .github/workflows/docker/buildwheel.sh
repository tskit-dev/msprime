#!/bin/bash
DOCKER_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$DOCKER_DIR/shared.env"

set -e -x

ARCH=`uname -p`
echo "arch=$ARCH"
# The GSL version available from yum install is too old so we manually install.
curl -o gsl-${GSL_VERSION}.tar.gz "ftp://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz"
tar -zxf gsl-${GSL_VERSION}.tar.gz
cd gsl-${GSL_VERSION}
./configure --prefix=/usr
make -j 2
make install
cd ..

for V in "${PYTHON_VERSIONS[@]}"; do
    PYBIN=/opt/python/$V/bin
    rm -rf build/       # Avoid lib build by narrow Python is used by wide python
    # Instead of letting setup.py install a newer numpy we install it here
    # using the oldest supported version for ABI compatibility
    $PYBIN/pip install oldest-supported-numpy
    $PYBIN/python setup.py build_ext --inplace
    $PYBIN/python setup.py bdist_wheel
done

cd dist
for whl in *.whl; do
    auditwheel -v repair "$whl"
    rm "$whl"
done
