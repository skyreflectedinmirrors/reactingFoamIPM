# user defined input paths

# path to cantera installation
export CT_PATH="$CONDA_PREFIX/lib/"
export CT_INC_PATH="$CONDA_PREFIX/include/"
# opencl header dir
export CL_INC_DIR=/opt/opencl-headers/
# path to sundials installation, leave blank if Cantera was compiled with a private
# sundials copy
export SUN_PATH="$CONDA_PREFIX/lib/"
export SUN_LIBS=""
if [ ! -z $SUN_PATH ]; then
    # sundials libraries for version >= 3.0
    export SUN_LIBS="-lsundials_sunlinsollapackband -lsundials_sunmatrixdense\
 -lsundials_sunmatrixband"
    # sundials libraries for version < 3.0
    # export SUN_LIBS="-lsundials_nvecserial -lsundials_cvodes"
fi


CF_PATH=`pwd`/chemFoamIPM/ \
WRAPPER_PATH=`pwd`/BatchedODE/AccelerInt/ \
PYJAC_PATH=`pwd`/include/ \
PYJAC_LIB_PATH=`pwd`/lib/ \
ACCELERINT_PATH=`pwd`/include/ \
ACCELERINT_LIB_PATH=`pwd`/lib/ \
./Allwmake -j12
