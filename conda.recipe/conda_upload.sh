PKG_NAME=aeolus
USER=dennissergeev

PYCODE="import ${PKG_NAME}; print(${PKG_NAME}.__version__)"
export VERSION=`python -c "${PYCODE}"`

# OS=$TRAVIS_OS_NAME-64
export CONDA_BLD_PATH=$HOME/conda-bld
mkdir -p $CONDA_BLD_PATH
conda config --set anaconda_upload no
conda config --add channels conda-forge
conda build --no-test .
PKG_FULL_NAME=`conda build --output .`
if [[ $VERSION == *"+"* ]]; then
    LABEL="nightly";
else
    LABEL="main";
fi;
echo ""
echo "label: $LABEL"
echo "version: $VERSION"
echo ""
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l $LABEL --force $PKG_FULL_NAME
