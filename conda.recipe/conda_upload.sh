PKG_NAME=aeolus
USER=dennissergeev

OS=$TRAVIS_OS_NAME-64
mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
export VERSION=`python -c 'import aeolus; print(aeolus.__version__)'`
conda build --no-test .
PKG_FULL_NAME=`conda build --output .`
# echo 'PKG_FULL_NAME='$PKG_FULL_NAME
if [[ $VERSION == *"+"* ]]; then
    LABEL="nightly";
else
    LABEL="main";
fi;
echo ""
echo "label: $LABEL"
echo "version: $VERSION"
echo ""
anaconda --token $CONDA_UPLOAD_TOKEN upload -u $USER -l $LABEL --force $PKG_FULL_NAME
