# reactingFoamIPM
An OpenFOAM solver based on reactingFoam, with profiling enabled by the IPM library and vectorized ODE integration.


Generating the pyJac files
==========================
```
python -m pyjac -l opencl -w 8 -se -p intel -jf full -o C -b reactingFoamIPM/include/opencl -i foam/bluffbody_LES/non-reacting/LES/chemkin/grimech30.dat -t foam/bluffbody_LES/non-reacting/LES/chemkin/thermo30.dat
```

Generating the pyJac library
============================
python -m pyjac.libgen -l opencl -so /store/ncurtis/reactingFoamIPM/include/opencl -out /store/ncurtis/reactingFoamIPM/lib/


Generating the accelerInt library
=================================
```
cp -r mechanism/ /store/ncurtis/reactingFoamIPM/include/
scons mechanism_dir=/store/ncurtis/reactingFoamIPM/include/ prefix=/store/ncurtis/reactingFoamIPM/ opencl-wrapper
scons install-opencl
```

Building the reactingFoamIPM exe for the model
==============================================
gri:
```
BASE_PATH=`pwd` \
CL_INC_DIR=/opt/opencl-headers/ \
PYJAC_PATH=`pwd`/include/opencl/ \
PYJAC_LIB_PATH=`pwd`/lib/ \
ACCELERINT_PATH=`pwd`/include \
ACCELERINT_LIB_PATH=`pwd`/lib/ \
MODEL_NAME=gri ./Allwmake -j12
```

skeletal:
```
BASE_PATH=`pwd` \
CL_INC_DIR=/opt/opencl-headers/ \
PYJAC_PATH=`pwd`/skeletal/include/opencl/ \
PYJAC_LIB_PATH=`pwd`/skeletal/lib/ \
ACCELERINT_PATH=`pwd`/skeletal/include \
ACCELERINT_LIB_PATH=`pwd`/skeletal/lib/ \
MODEL_NAME=skeletal ./Allwmake -j12
```
