#!/bin/bash
################
# WHERE WE ARE #
################
COPPSHOME=$(pwd)
#############
## VERSION ##
#############
echo "Please, select what to do:"
OPTIONS="build-serial build-parallel clean exit"
select opt in $OPTIONS; do
	if [ "$opt" = "build-serial" ]; then
		CC="g++"
		FF="gfortran"
		COMPOPT="-O2 -funroll-loops -frerun-loop-opt"
		break
	elif [ "$opt" = "build-parallel" ]; then
		CC="mpicxx"
		FF="mpif77"
		COMPOPT="-O2 -funroll-loops -frerun-loop-opt -D_MPI_"
		break
	elif [ "$opt" = "clean" ]; then
		echo "This operation cleans all libraries and executables leaving only the source code. Continue? [yes | no] "
		read answer
		if [ "$answer" == "yes" ]; then

			cd ${COPPSHOME}/src
			make clean

			cd ${COPPSHOME}/lib/
			rm ./*.a

			cd src
			rm -rf cminpack-1.0.2 cminpack_jac-1.0.2 cqp cubature-1.0 jemdi-v1.2 levmar-2.4 sparselib_1_6

			cd ../../

			rm ${COPPSHOME}/Makefile.in

			echo ""
			echo "C++OPPS build dir completely clean."
			echo ""
			exit
		else
			echo ""
			echo "Nothing to be done. Build script stopped."
			echo ""
			exit
		fi
		break
	elif [ "$opt" = "exit" ]; then
		echo "Build script interrupted"
		exit
		break
	else
		echo "Bad option - please answer 1 to build, 2 to clean the package, or 3 to exit"
	fi
done
#####################
# WHERE IS LAPACK? ##
#####################
read -p "Enter the location of liblapack.so and libblas.so libraries [/usr/lib64/]: " dir
dir=${dir:-/usr/lib64/}
export LAPACKDIR=$dir
###########################################
# PREPARATION OF THE make.COMPILERS FILE ##
###########################################
cat > ./Makefile.in << EOF
CC      = ${CC}
FF      = ${FF}
FLAG    = ${COMPOPT}
FLAGF   = ${COMPOPT}
LL=-L$(pwd)/lib/ -L$(pwd)/lib/src/jemdi-v1.2/lib/ -L$dir
II=-I$(pwd)/include -I$(pwd)/lib/src/jemdi-v1.2/include/
LIBS=-std=c++11 -ldite2 -lzmat -larmadillo -D_ARMA_USE_LAPACK -llapack -lblas
FORTLIB=-lgfortran #-lifcore
LAPACKLIB = -llapack -lblas
EOF
########################
# BUILD THE CQP IBRARY #
########################
cd ${COPPSHOME}/lib/src/
tar zxf cqp.tar.gz
cd cqp
make -j4
CQPPATH=$(pwd)
# Check if CQP has been built
log=`find -name "libcqp++.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the CQP library (lib/src/cqp/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
cp libcqp++.a ../../
##################################
# BUILD THE CMINPACK_JAC LIBRARY #
##################################
cd ${COPPSHOME}/lib/src/
tar zxf cminpack_jac-1.0.2.tgz
cd cminpack_jac-1.0.2/
make -j4
# Check if CMINPACK_JAC has been built
log=`find -name "libminpack_jac.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the CMINPACK_JAC library (lib/src/cminpack_jac-1.0.2/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
cp libminpack_jac.a ../../
##################################
# BUILD THE SPARSE LIBRARY #
##################################
cd ${COPPSHOME}/lib/src/
tar zxf sparselib_1_6-mod.tar.gz
cd sparselib_1_6/
make all -j4
# Check if SPARSE has been built
cd lib
log=`find -name "libsparse.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the SPARSE library (lib/src/sparselib_1_6/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
log=`find -name "libspblas.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the SPARSE library (lib/src/sparselib_1_6/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
log=`find -name "libmv.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the SPARSE library (lib/src/sparselib_1_6/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
cd ..
cp lib/*.a ../../
#############################
# BUILD THE LEVAMAR LIBRARY #
#############################
cd ${COPPSHOME}/lib/src/
tar zxf levmar-2.4.tgz
cd levmar-2.4/
make -j4
# Check if LEVMAR has been built
log=`find -name "liblevmar.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the LEVMAR library (lib/src/levmar-2.4/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
cp liblevmar.a ../../
##############################
# BUILD THE CUBATURE LIBRARY #
##############################
cd ${COPPSHOME}/lib/src/
tar zxf cubature-1.0.tgz
cd cubature-1.0/
make -j4
# Check if CUBATURE has been built
log=`find -name "libcubature.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the CUBATURE library (lib/src/cubature-1.0/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
cp libcubature.a ../../
###########################
# BUILD THE JEMDI LIBRARY #
###########################
cd ${COPPSHOME}/lib/src/
tar zxf jemdi-v1.2.tgz
cd jemdi-v1.2/
cat > ./build_input << EOF
2
1
EOF
sh ./build.sh < build_input
# Check if JEMDI has been built
cd lib
log=`find -name "libcmrg.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the JEMDI/SPRNG library (lib/src/jemdi-1.2/sprng/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
log=`find -name "libjemdi.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the JEMDI library (lib/src/jemdi-1.2/. Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
cd ..
#################
# BUILD C++OPPS #
#################
cd ${COPPSHOME}/src/
make -j4
# Check if C++OPPS has been built
log=`find -name "copps"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the C++OPPS-v2.2 stand alone (src/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
cd ../
#######
# END #
#######
echo ""
echo "====================================================================================================================="
echo " ** C++OPPS v2.2 succesfully compiled ** "
echo ""
echo " The 'copps' executable is in the $(pwd)/src/ directory"
echo ""
echo "Please, find in the $(pwd)/doc/input_examples/ subdirectory some"
echo "sample scripts useful to configure your first C++OPPS run."
echo ""
echo "====================================================================================================================="
echo ""

