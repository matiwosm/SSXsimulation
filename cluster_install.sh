#
# Dedalus Installation Script for clusters
#
# This script is designed to create a fully isolated Python installation
# with the dependencies you need to run dedalus.
#
# This script is designed to install on PSC/Bridges and is based on the vanilla Dedalus install script.
#
# There are a few options, but you only need to set *one* of them.  And
# that's the next one, DEST_DIR.  But, if you want to use an existing HDF5
# installation you can set HDF5_DIR, or if you want to use some other
# subversion checkout of Dedalus, you can set DEDALUS_DIR, too.  (It'll already
# check the current directory and one up.
#
date

DEST_SUFFIX="dedalus"
DEST_DIR="`pwd`/${DEST_SUFFIX/ /}"   # Installation location
BRANCH="tip" # This is the branch to which we will forcibly update.

MAKE_PROCS="-j" # change this to set the number of cores you use to build packages

# If you need to supply arguments to the NumPy or SciPy build, supply them here
# This one turns on gfortran manually:
#NUMPY_ARGS="--fcompiler=gnu95"
# If you absolutely can't get the fortran to work, try this:
#NUMPY_ARGS="--fcompiler=fake"
#NUMPY_ARGS="--fcompiler=fake"

# Packages

INST_OPENMPI=0 # by default, don't build OpenMPI. If you're on linux, use your package manager.
INST_HDF5=1 # by default, do build HDF5.
INST_ATLAS=0 # by default, we will not build our own ATLAS. If you're on OSX, you'll want to use accelerate anyway.
INST_SCIPY=1
INST_IPYTHON=0 # by default, don't build ipython
INST_FTYPE=0 # by default, don't install freetype
INST_PNG=0 # by default, don't install libpng
INST_PKGCFG=0 # by default, don't install pkg-config
INST_OPENSSL=0 #by default, don't install openssl
INST_ZLIB=0 # by default, don't install zlib
INST_MERCURIAL=0 # by default, don't install mercurial

## Packages to install from source and version number
PYTHON_VERSION='3.6.4'
PYTHON_VERSION_SHORT='3.6'
FFTW_VERSION='3.3.7'
NUMPY_VERSION='1.14.0'
SCIPY_VERSION='1.0.0'
HDF5_VERSION_MAJOR='1.10'
HDF5_VERSION=$HDF5_VERSION_MAJOR".1" 
MERCURIAL_VERSION='3.9.1'

if [ ${REINST_DEDALUS} ] && [ ${REINST_DEDALUS} -eq 1 ] && [ -n ${DEDALUS_DEST} ]
then
    DEST_DIR=${DEDALUS_DEST}
fi
# Make sure we are NOT being run as root
if [[ $EUID -eq 0 ]]
then
   echo "******************************************************"
   echo "*                                                    *"
   echo "*                                                    *"
   echo "*  IT IS A BAD IDEA TO RUN THIS SCRIPT AS ROOT!!!!   *"
   echo "*                                                    *"
   echo "*                                                    *"
   echo "******************************************************"
   echo
   echo "If you really want to do this, you must manually edit"
   echo "the script to re-enable root-level installation.  Sorry!"
   exit 1
fi
if [[ ${DEST_DIR%/} == /usr/local ]]
then
   echo "******************************************************"
   echo "*                                                    *"
   echo "*                                                    *"
   echo "*  THIS SCRIPT WILL NOT INSTALL TO /usr/local !!!!   *"
   echo "*                                                    *"
   echo "*                                                    *"
   echo "******************************************************"
   exit 1
fi

if type -P wget &>/dev/null
then
    echo "Using wget"
    export GETFILE="wget -nv"
else
    echo "Using curl"
    export GETFILE="curl -sSO"
fi

if type -P sha512sum &> /dev/null
then
    echo "Using sha512sum"
    export SHASUM="sha512sum"
elif type -P shasum &> /dev/null
then
    echo "Using shasum -a 512"
    export SHASUM="shasum -a 512"
else
    echo
    echo "Unable to locate any shasum-like utility."
    echo "ALL FILE INTEGRITY IS NOT VERIFIABLE."
    echo "THIS IS PROBABLY A BIG DEAL."
    echo
    echo "(I'll hang out for a minute for you to consider this.)"
    sleep 60
fi

function help_needed
{
    [ -e $1 ] && return
    echo
    echo "WE NEED YOUR HELP!" 
    echo
    echo "We do not have testing facilities on $1 yet. While we're confident "
    echo "that Dedalus will work on $1, to get it working, "
    echo "we need to know the list of packages from $1 Dedalus requires."
    echo
    echo "If you are familiar with $1, please look over the package list for Ubuntu"
    echo "in this install script and help us translate it into $1 packages."
    echo
    echo "If you'd like to help, please don't hesitate to email the dev list,"
    echo "dedalus-dev@googlegroups.com"
    echo
    echo "   --the Dedalus team"
    echo
    echo
}

function get_dedalusproject
{
    [ -e $1 ] && return
    echo "Downloading $1 from dedalus-project.org"
    ${GETFILE} "http://dedalus-project.org/dependencies/$1" || do_exit
    ( ${SHASUM} -c $1.sha512 2>&1 ) 1>> ${LOG_FILE} || do_exit
}
function get_file
{
    [ -e $2 ] && return
    echo "Downloading $2 from $1"
    ${GETFILE} "$1/$2" || do_exit
}
function git_repo
{
    [ -e $2 ] && return
    echo "Gitting $2 from $1"
    git clone "$1" "$2" || do_exit
}
function do_setup_py
{
    [ -e $1/done ] && return
    LIB=$1
    shift
    if [ -z "$@" ]
    then
        echo "Installing $LIB"
    else
        echo "Installing $LIB (arguments: '$@')"
    fi
    [ ! -e $LIB/extracted ] && tar xfz $LIB.tar.gz
    touch $LIB/extracted
    cd $LIB
    if [ ! -z `echo $LIB | grep h5py` ]
    then
	( ${DEST_DIR}/bin/python3 setup.py build --hdf5=${HDF5_DIR} $* 2>&1 ) 1>> ${LOG_FILE} || do_exit
    else
        ( ${DEST_DIR}/bin/python3 setup.py build   $* 2>&1 ) 1>> ${LOG_FILE} || do_exit
    fi
    ( ${DEST_DIR}/bin/python3 setup.py install    2>&1 ) 1>> ${LOG_FILE} || do_exit
    touch done
    cd ..
}
function do_exit
{
    echo "********************************************"
    echo "        FAILURE REPORT:"
    echo "********************************************"
    echo
    tail -n 10 ${LOG_FILE}
    echo
    echo "********************************************"
    echo "********************************************"
    echo "Failure.  Check ${LOG_FILE}.  The last 10 lines are above."
    exit 1
}

function host_specific
{
    MYHOST=`hostname -s`  # just give the short one, not FQDN
    MYHOSTLONG=`hostname --long` # FQDN
    MYOS=`uname -s`       # A guess at the OS
    #BLAS="/usr/lib/"
    #LAPACK="/usr/lib/"
    HDF5_USE_GCC=0
    MATPLOTLIB_USE_CPP=0
    INSECURE_HG=0
    PY_IPO=""

    if [[ $MYHOSTLONG == *"bridges.psc.edu" ]]
    then
	echo "Looks like you are running on PSC/Bridges"
	echo
	echo "We assume here that you have the intel compiler suite loaded,"
	echo "and a module list like the following (versions omitted):"
	echo "Currently Loaded Modulefiles:"
        echo "  1) psc_path    2) slurm   3) zlib      4) icc"
	echo
	echo "Here are your module files:"
	module list
	echo
	MPICC="mpiicc"
	MPICXX="mpiicpc"
	MPIF90="mpiifort"
	ARCH_CONF="-axCORE-AVX2 -xSSE4.2"
	INST_MERCURIAL=1
	
	#critical for mpi4py install; otherwise we hit GCC
	export I_MPI_CC=icc

	IS_PSC_BRIDGES=1

	MKL_LIB="/opt/packages/intel/compilers_and_libraries/linux/mkl/lib/intel64"
	MKL_INCLUDE="/opt/packages/intel/compilers_and_libraries/linux/mkl/include"

	# points to IntelMPI location
	MPI_PATH=$I_MPI_ROOT

	BUILD_TIME="an hour"
	
    elif [[ $MYHOSTLONG == *"rc."*"colorado.edu" ]]
    then
	echo "Looks like you are running on CU/Summit"
	echo
	echo "We assume here that you have the intel compiler suite loaded,"
	echo "and a module list like the following (versions omitted):"
	echo "Currently Loaded Modulefiles:"
        echo "  1) intel/intel-version      2) impi/impi-version   3) mkl/mkl-version"
	echo
	echo "Here are your module files:"
	module list
	echo
	echo
	MPICC="mpiicc"
	MPICXX="mpiicpc"
	MPIF90="mpiifort"
	CC="mpiicc"

	ARCH_CONF="-axCORE-AVX2 -xSSE4.2"
	INST_MERCURIAL=1
	
	IS_CU_SUMMIT=1
	PY_IPO=""

	MKL_LIB=$MKLROOT"/lib/intel64/"
	MKL_INCLUDE=$MKLROOT"/include"

	# points to IntelMPI location
	MPI_PATH=$I_MPI_ROOT
	
	# use insecure wget
#	export GETFILE="wget -nv --no-check-certificate"

	BUILD_TIME="an unkown time"
	
        elif [[ $MYHOSTLONG == "hyades.ucsc.edu" ]]
        then
	echo "Looks like you are running on UCSC/Hyades"
	echo
	echo "We assume here that you have the intel compiler suite loaded,"
	echo "and a module list like the following (versions omitted):"
	echo "Currently Loaded Modulefiles:"
        echo "  1) impi/impi-version      2) intel/intel-version"
	echo
	echo "Here are your module files:"
	module list
	echo
	echo
	MPICC="mpiicc"
	MPICXX="mpiicpc"
	MPIF90="mpiifort"
	CC="mpiicc"
        export MPICC_CC=icc
	export MPICXX_CXX=icpc

	ARCH_CONF="-xAVX"
	
	IS_CU_SUMMIT=1
	PY_IPO=""
        INSECURE_HG=1
        
	MKL_LIB=$MKLROOT"/lib/intel64/"
	MKL_INCLUDE=$MKLROOT"/include"

	# points to IntelMPI location
	MPI_PATH=$I_MPI_ROOT
	
	BUILD_TIME="an unkown time"
	
    elif [[ $MYHOSTLONG == *"pfe"* && $MYHOSTLONG == *"nas.nasa.gov" ]]
    then
	echo "Looks like you are running on NAS/Pleaides (NASA)"
	echo
	echo "We assume here that you have the intel compiler suite loaded,"
	echo "and a module list like the following (versions omitted):"
	echo "Currently Loaded Modulefiles:"
        echo "  	1) mpi-sgi/mpt.version     2) mpi-sgi/mpt             3) comp-intel/version   4) git/version               5) openssl/version"
	echo
	echo "Here are your module files:"
	module list
	echo

	# proper wrappers for using Intel rather than GNU compilers,
	# Thanks to Daniel Kokron at NASA.
	export MPICC_CC=icc
	export MPICXX_CXX=icpc

	export CC=mpicc

	MPICC="mpicc"
	MPICXX="mpicxx"
	MPIF90="mpif90"
	ARCH_CONF="-axCORE-AVX512 -xSSE4.2"
	HDF5_USE_GCC=1
	MATPLOTLIB_USE_CPP=1
	INST_MERCURIAL=0
	
	IS_NAS_PLEIADES=1
	
	MKL_LIB=$MKLROOT"/lib/intel64/"
	MKL_INCLUDE=$MKLROOT"/include"

	# points to MPI-SGI location
	MPI_PATH=$MPI_ROOT

	# avoid contention for forking on login nodes
	MAKE_PROCS="" 

	# newer versions of HDF5 are core dumping
	#HDF5_VERSION='1.8.15-patch1'

	BUILD_TIME="45 minutes"

	if [ ! -e ~/.pip/pip.conf ]
	then
	    [ ! -e ~/.pip ] && mkdir ~/.pip
	    cat >> ~/.pip/pip.conf << EOF
[global]
cert = /etc/ssl/ca-bundle.pem
EOF
	fi
    elif [[ $HOSTNAME == "mira"* ]]
    then
	echo "Looks like you are running on ANL/Mira (DOE)"
	echo

	MPICC="mpicc"
	MPICXX="mpicxx"
	MPIF90="mpif90"
	ARCH_CONF=""
	
	IS_ANL_MIRA=1
	
	BUILD_TIME="unknown"
	
    fi
    PYCONF_ARGS="CC="$MPICC"  CXX="$MPICXX" F90="$MPIF90
    #some parts of next line critical for compiling on Janus.  Deprecate?
    PYCONF_ARGS=$PYCONF_ARGS" --with-cxx-main="$MPICXX" --with-gcc="$MPICC" --enable-shared  --with-system-ffi"
    PY_LDFLAGS="-lpthread"
    PY_OPT="-w -vec-report0 -qopt-report0"
    PY_FLAGS="-mkl -O3 "$PY_IPO" "$ARCH_CONF" -fPIC"
    FFTWCONF_ARGS="CC="$MPICC" \
                   CXX="$MPICXX" \
                   F77="$MPIF90" \
                   MPICC="$MPICC" MPICXX="$MPICXX" \
                   --enable-shared \
                   --enable-mpi --enable-openmp --enable-threads"
    
    FFTW_FLAGS="-O3 "$ARCH_CONF
    HDF5_FLAGS="-w" # suppress gcc warnings
    HDF5CONF_ARGS="CC="$MPICC" CXX="$MPICXX" F77="$MPIF90" MPICC="$MPICC" MPICXX="$MPICXX" CFLAGS="$HDF5_FLAGS" CPPFLAGS="$HDF5_FLAGS" --enable-shared --enable-parallel"

    echo
    echo "If this is correct, then we are ready to continue."
    echo
    echo "Building the stack will take about $BUILD_TIME to complete."
    echo
    
    if [ ! -z "${CFLAGS}" ]
    then
        echo "******************************************"
        echo "******************************************"
        echo "**                                      **"
        echo "**    Your CFLAGS is not empty.         **"
        echo "**    This can break h5py compilation.  **"
        echo "**                                      **"
        echo "******************************************"
        echo "******************************************"
    fi
}
export MPI_PATH

ORIG_PWD=`pwd`
echo "+++++++++"
echo "Greetings, human. Welcome to the Dedalus Install Script for Clusters"
echo "            (warning: hand tuned.  Your mileage may vary)"
echo
host_specific
echo
echo
read -p "[hit enter] "
echo
echo

LOG_FILE="${DEST_DIR}/dedalus_install.log"

mkdir -p ${DEST_DIR}/src
cd ${DEST_DIR}/src

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${DEST_DIR}/lib/"

PYTHON='Python-'$PYTHON_VERSION
FFTW='fftw-'$FFTW_VERSION
NUMPY='numpy-'$NUMPY_VERSION
SCIPY='scipy' # git repo now; version handled by checking out specific tag
HDF5='hdf5-'$HDF5_VERSION
HDF5_MAJOR='hdf5-'$HDF5_VERSION_MAJOR
MERCURIAL='mercurial-'$MERCURIAL_VERSION

# get the files
get_file https://www.python.org/ftp/python/$PYTHON_VERSION/ $PYTHON.tgz
get_file http://www.fftw.org/ $FFTW.tar.gz

#git_repo git://github.com/numpy/numpy.git numpy 
#git_repo git://github.com/scipy/scipy.git scipy

#[ $INST_OPENMPI -eq 1 ] && get_dedalusproject $OPENMPI.tar.gz

[ $INST_HDF5 -eq 1 ] && get_file http://www.hdfgroup.org/ftp/HDF5/releases/$HDF5_MAJOR/$HDF5/src/ $HDF5.tar.gz
[ $INST_MERCURIAL -eq 1 ] && get_file https://www.mercurial-scm.org/release/ $MERCURIAL.tar.gz

# first, OpenMPI, if we're doing that
# this is a bit tricky on clusters, and needs to be done on the compute nodes to ensure that
# the right memory space is mapped (e.g., Pleiades issues early on); this is also where the
# largest inter-cluster variability occurs in compile flags/options/etc.
#
# for now, plan to handle this outside of the script.
if [ $INST_OPENMPI -eq 1 ]
then
    if [ ! -e $OPENMPI/done ]
    then
        [ ! -e $OPENMPI ] && tar xfz $OPENMPI.tar.gz
        echo "Installing OPENMPI"
        cd $OPENMPI
        ( ./configure CPPFLAGS=-I${DEST_DIR}/include CFLAGS=-I${DEST_DIR}/include --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    OPENMPI_DIR=${DEST_DIR}
    export LDFLAGS="${LDFLAGS} -L${OPENMPI_DIR}/lib/ -L${OPENMPI_DIR}/lib64/"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${OPENMPI_DIR}/lib/"
    PATH="${OPENMPI_DIR}/bin:${PATH}"
    export MPI_PATH=${OPENMPI_DIR}
fi


# python3 

if [ ! -e $PYTHON/done ]
then
    echo "Installing Python."
    [ ! -e $PYTHON ] && tar xfz $PYTHON.tgz
    cd $PYTHON
    if [ $INST_OPENSSL -eq 1 ]
    then
       export PY_CFLAGS="-I${DEST_DIR}/include"
    fi
    echo "PYCONF_ARGS = $PYCONF_ARGS"
    echo "PY_FLAGS = $PY_FLAGS"
    echo "PY_LDFLAGS = $PY_LDFLAGS"
    
    ( ./configure --prefix=${DEST_DIR}/ ${PYCONF_ARGS} OPT="$PY_OPT" FOPT="$PY_OPT" \
		  CFLAGS="$PY_FLAGS" CPPFLAGS="$PY_FLAGS" F90FLAGS="$PY_FLAGS" LDFLAGS="$PY_LDFLAGS" 2>&1 ) 1>> ${LOG_FILE} || do_exit

    ( make ${MAKE_PROCS} 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
    ( ${DEST_DIR}/bin/pip3 install --upgrade pip 2>&1) 1>> ${LOG_FILE} || do_exit
    touch done
    cd ..

fi



export PYTHONPATH=${DEST_DIR}/lib/python$PYTHON_VERSION_SHORT/site-packages/


# FFTW3
if [ ! -e $FFTW/done ]
then
    echo "Installing FFTW."
    [ ! -e $FFTW ] && tar xzf $FFTW.tar.gz 
    cd $FFTW

    ( ./configure --prefix=${DEST_DIR}/ ${FFTWCONF_ARGS} CFLAGS="$FFTW_FLAGS" CPPFLAGS="$FFTW_FLAGS" F90FLAGS="$FFTW_FLAGS" 2>&1 ) 1>> ${LOG_FILE} || do_exit

    ( make ${MAKE_PROCS} 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
    # what is this line?
    ( ln -sf ${DEST_DIR}/bin/python2.7 ${DEST_DIR}/bin/pyyt 2>&1 ) 1>> ${LOG_FILE}
    ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
    touch done
    cd ..

fi
export FFTW_PATH=${DEST_DIR}/

# HDF5, if we're doing that.
if [ $INST_HDF5 -eq 1 ]
then
    if [ ! -e $HDF5/done ]
    then
	echo "Installing HDF5"
        [ ! -e $HDF5 ] && tar xzf $HDF5.tar.gz
        cd $HDF5
	if [ $HDF5_USE_GCC -eq 1 ]
	then
	    # fall back to GCC; otherwise mpi-sgi stack compile fails on hdf5
	    OLD_MPICC_CC=$MPICC_CC
	    OLD_MPICXX_CXX=$MPICXX_CXX
	    export MPICC_CC=
	    export MPICXX_CXX=
	fi
	echo "HDF5CONF_ARGS=$HDF5CONF_ARGS"
        ( ./configure --prefix=${DEST_DIR}/ ${HDF5CONF_ARGS} 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
	if [ $HDF5_USE_GCC -eq 1 ]
	then
	    export MPICC_CC=$OLD_MPICC_CC
	    export MPICXX_CXX=$OLD_MPICXX_CXX
	fi
    fi
    export HDF5_DIR=${DEST_DIR}

fi


# cython
echo "pip installing cython."
( ${DEST_DIR}/bin/pip3 install cython 2>&1 ) 1>> ${LOG_FILE} || do_exit

# numpy
if [ ! -e ~/.numpy-site.cfg ]
then
    cat >> ~/.numpy-site.cfg << EOF
[mkl]
library_dirs = $MKL_LIB
include_dirs = $MKL_INCLUDE
mkl_libs = mkl_rt
lapack_libs =
EOF
fi
echo "pip installing numpy"
( ${DEST_DIR}/bin/pip3 install numpy --no-binary :all: 2>&1 ) 1>> ${LOG_FILE} || do_exit

echo "Executing numpy speed test"
cat >> numpy_speed_test.py << EOF
import numpy
import sys
import timeit

print("numpy version:", numpy.__version__)
print("     location:", numpy.__file__)
print()
print(numpy.__config__.show())
print()
print("numpy speed test")
setup = "import numpy; x = numpy.random.random((1000,1000))"
count = 100

t = timeit.Timer("numpy.dot(x, x.T)", setup=setup)
print("dot:", t.timeit(count)/count, "sec")
print()
EOF
${DEST_DIR}/bin/python3 numpy_speed_test.py

#scipy
echo "pip installing scipy"
( ${DEST_DIR}/bin/pip3 install scipy 2>&1 ) 1>> ${LOG_FILE} || do_exit

# nose
echo "pip installing nose."
( ${DEST_DIR}/bin/pip3 install nose 2>&1 ) 1>> ${LOG_FILE} || do_exit

# mpi4py
echo "pip installing mpi4py."
( ${DEST_DIR}/bin/pip3 install mpi4py 2>&1 ) 1>> ${LOG_FILE} || do_exit

export HDF5_DIR=${DEST_DIR}

echo "pip installing h5py."
echo "HDF5_DIR=$HDF5_DIR"
( ${DEST_DIR}/bin/pip3 install h5py 2>&1 ) 1>> ${LOG_FILE} || do_exit

# matplotlib
PATH=$DEST_DIR/bin/:$PATH
if [ $INST_FTYPE -eq 1 ]
then
    echo "manually installing matplotlib."
# Now we set up the basedir for matplotlib:
    mkdir -p ${DEST_DIR}/src/$MATPLOTLIB
    echo "[directories]" >> ${DEST_DIR}/src/$MATPLOTLIB/setup.cfg
    echo "basedirlist = ${DEST_DIR}" >> ${DEST_DIR}/src/$MATPLOTLIB/setup.cfg
    if [ `uname` = "Darwin" ]
    then
	echo "[gui_support]" >> ${DEST_DIR}/src/$MATPLOTLIB/setup.cfg
	echo "macosx = False" >> ${DEST_DIR}/src/$MATPLOTLIB/setup.cfg
    fi
    
    if [ ! -e $MATPLOTLIB/extracted ]
    then
	tar xfz $MATPLOTLIB.tar.gz
	touch $MATPLOTLIB/extracted
    fi
    cd $MATPLOTLIB
    patch -b setupext.py <<EOF
960c960
<             'freetype2', 'ft2build.h',
---
>             'freetype2', 'freetype2/ft2build.h',
EOF
    cd ..
    do_setup_py $MATPLOTLIB
    
else
    echo "pip installing matplotlib."
    if [ $MATPLOTLIB_USE_CPP -eq 1 ]
    then
	# workaround for qhull compilation failure; force intel CC to be C++ compiler
	OLD_CC=$CC
	export CC=icpc
    fi
    ( ${DEST_DIR}/bin/pip3 install matplotlib 2>&1 ) 1>> ${LOG_FILE} || do_exit
    [ $MATPLOTLIB_USE_CPP -eq 1 ] && export CC=$OLD_CC
fi

# ipython
if [ $INST_IPYTHON -eq 1 ]
then
    echo "pip installing ipython."
    ( ${DEST_DIR}/bin/pip3 install ipython 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( ${DEST_DIR}/bin/pip3 install pyzmq 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( ${DEST_DIR}/bin/pip3 install jinja2 2>&1 ) 1>> ${LOG_FILE} || do_exit
fi

# We assume that hg can be found in the path.
if type -P hg &>/dev/null
then
    export HG_EXEC=hg
else
    if [ $INST_MERCURIAL -eq 1 ]
    then
	if [ ! -e $MERCURIAL/done ]
	then
	    if [ ! -e $MERCURIAL/extracted ]
	    then
               [ ! -e $MERCURIAL ] && tar xzf $MERCURIAL.tar.gz
	       touch $MERCURIAL/extracted
	    fi
            echo "Installing mercurial"
            cd $MERCURIAL
	    module swap intel gcc >> ${LOG_FILE} || do_exit
	    module load python >> ${LOG_FILE} || do_exit
	    module list
	    OLD_CC=$CC
	    export CC=gcc
	    ( make install PREFIX=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit
	    module unload python >> ${LOG_FILE} || do_exit
	    module swap gcc intel >> ${LOG_FILE} || do_exit
	    export CC=$OLD_CC
	    touch done
	    cd ..
	else
	    echo "Cannot find mercurial.  Please make sure it is installed, or request that the script build mercurial."
	    do_exit
	fi
    fi
    export HG_EXEC=hg
fi

if [ -z "$DEDALUS_DIR" ]
then
    DEDALUS_DIR="$PWD/dedalus/"
    if [ ! -e dedalus ]
    then
	if [ $INSECURE_HG -eq 1 ]
	then
	    echo "WARNING: INSECURE HG (check SSL)"
	    echo "WARNING: INSECURE HG (check SSL)" >> ${LOG_FILE}
	    ( ${HG_EXEC} --debug clone --insecure https://bitbucket.org/dedalus-project/dedalus/ dedalus 2>&1 ) 1>> ${LOG_FILE}
	else
            ( ${HG_EXEC} --debug clone https://bitbucket.org/dedalus-project/dedalus/ dedalus 2>&1 ) 1>> ${LOG_FILE}
	fi
        # Now we update to the branch we're interested in.
        ( ${HG_EXEC} -R ${DEDALUS_DIR} up -C ${BRANCH} 2>&1 ) 1>> ${LOG_FILE}
    fi
    echo Setting DEDALUS_DIR=${DEDALUS_DIR}
fi


## afterwards
# Add the environment scripts
( cp ${DEDALUS_DIR}/docs/activate ${DEST_DIR}/bin/activate 2>&1 ) 1>> ${LOG_FILE}
sed -i.bak -e "s,__DEDALUS_DIR__,${DEST_DIR}," ${DEST_DIR}/bin/activate
( cp ${DEDALUS_DIR}/docs/activate.csh ${DEST_DIR}/bin/activate.csh 2>&1 ) 1>> ${LOG_FILE}
sed -i.bak -e "s,__DEDALUS_DIR__,${DEST_DIR}," ${DEST_DIR}/bin/activate.csh

echo "Doing Dedalus update, wiping local changes and updating to branch ${BRANCH}"
MY_PWD=`pwd`
cd $DEDALUS_DIR
( ${HG_EXEC} pull 2>1 && ${HG_EXEC} up -C 2>1 ${BRANCH} 2>&1 ) 1>> ${LOG_FILE}

echo "Installing Dedalus"

( export PATH=$DEST_DIR/bin:$PATH ; ${DEST_DIR}/bin/python3 setup.py build_ext --inplace 2>&1 ) 1>> ${LOG_FILE} || do_exit
touch done
cd $MY_PWD

if !( ( ${DEST_DIR}/bin/python3 -c "import readline" 2>&1 )>> ${LOG_FILE})
then
    echo "Installing pure-python readline"
    ( ${DEST_DIR}/bin/pip3 install readline 2>&1 ) 1>> ${LOG_FILE}
fi


function print_afterword
{
    echo
    echo
    echo "========================================================================"
    echo
    echo "dedalus is now installed in $DEST_DIR ."
    echo
    echo "To run from this new installation, use the activate script for this "
    echo "environment."
    echo
    echo "    $ source $DEST_DIR/bin/activate"
    echo
    echo "This modifies the environment variables DEDALUS_DEST, PATH, PYTHONPATH, and"
    echo "LD_LIBRARY_PATH to activate dedalus.  If you use csh, just"
    echo "append .csh to the above."
    echo
    echo "The source for dedalus is located at:"
    echo "    $DEDALUS_DIR"
    echo
    echo "For support, see the website and join the mailing list:"
    echo
    echo "    http://dedalus-project.org/"
    echo "    http://dedalus-project.readthedocs.org/       (Docs)"
    echo
    echo "    https://groups.google.com/forum/#!forum/dedalus-users"
    echo
    echo "========================================================================"
    echo
    echo "Good luck, and email the user list if you run into any problems."
}

print_afterword
print_afterword >> ${LOG_FILE}

echo "dedalus dependencies were last updated on" > ${DEST_DIR}/.dedalus_update
date >> ${DEST_DIR}/.dedalus_update

date
