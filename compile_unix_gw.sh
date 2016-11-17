#!/bin/bash
# compile_unix_gw.sh
#
# purpose:
# - compile libgwlib.so , a shared library needed by pest groundwater utilies
# - compile pest groundwater utilities
#
#
# usage:
# simply type "bash compile_unix_gw" in the folder where it come out
#
# needs:
# - src_utils src_lib and src_dep to be present in the same folder
# - ifort as F90 compiler , if you want tu use another do a find and replace ;)
#
# workflow:
# 1) creates ./bin ./obj and ./lib
# 2) compiles the dependencies defn.f90 and inter.f90
# 3) create the makefile for libgwlib.so
# 4) compiles it
# 5) create the makefile for the utilities
# 6) compiles them
# 7) wipe out the mess
#
#
# IMPORTANT NOTE:
# you need to put "./lib/libgwlib.so" somewhere in your LD_LIBRARY_PATH to use the utilities
#
#
# known issues:
# - it works with the intel ifort compiler, there where many 
#  compilations issues due to gfortran, if you want to use gfortran, 
#  expect spending many hours unsuccessfully like I did before ifort ;)
#
# v0.1 : 
# 	- first trial, let's see if it works for someone else...
#
# LEGAL DISCLAIMER:
# This software is an experimental research prototype. The software is
# provided by the copyright holders and contributors "as is" without
# warranty of any kind, express or implied, including, but not limited to,
# the implied warranties of merchantability and fitness for a particular
# purpose and non infringement. In no event shall the authors or copyright
# owner or contributors be liable for any claim, any direct, indirect,
# incidental, special, exemplary, or consequential damages (including, but
# not limited to, procurement of substitute goods or services; loss of use,
# data, or profits; or business interruption) however caused and on any
# theory of liability, whether in contract, strict liability, or tort
# (including negligence or otherwise) arising in any way out of the use of
# or in connection with this software, or other dealing in the software
# even if advised of the possibility of such damage.#
#
# 
# date: September 2011
# Author: Andrea Borghi, University of Neuchatel, Switzerland
#
# contact : andrea.borghi@unine.ch


###################################
###      DOEXECMAKEFILE        ####
###################################
function doExecMakeFile() {

printf "BINDIR=./bin\n"
printf "OBJDIR=./obj\n\n"
echo "F90  = ifort"
echo "F90 = ifort"

printf "\n"
printf "\n"

echo "LIB= libgwlib.so"

echo "FFLAGS = -O0 \$(OBJDIR)/defn.mod  \$(OBJDIR)/inter.mod"
echo "LDFLAGS= \$(LIB)"

printf "SRC="
for infile in $(ls ${SRCDIR}/*.f90)
do
	printf "\t%s" $infile
done
for infile in $(ls ${SRCDIR}/*.f)
do
	printf "\t%s" $infile
done

printf "\n"
printf "OBJ= "
for infile in $(ls ${SRCDIR}/*.f90)
do
	printf "\t\${OBJDIR}/%s.o" $(basename $infile .f90)
done

for infile in $(ls ${SRCDIR}/*.f)
do
	printf "\t\${OBJDIR}/%s.o" $(basename $infile .f)
done

printf "\n"
printf "EXEC= "
for infile in $(ls ${SRCDIR}/*.f90)
do
	printf "\t\${BINDIR}/%s" $(basename $infile .f90)
done

for infile in $(ls ${SRCDIR}/*.f)
do
	printf "\t\${BINDIR}/%s" $(basename $infile .f)
done

printf "\n"
printf "\n"
printf "all: \$(EXEC)\n\n"

# only .f file
printf "\${BINDIR}/linpos :\n\t\$(F90) \$(LDLAGS) -o \${OBJDIR}/linpos.o \$(OBJDIR)/defn.o \$(OBJDIR)/inter.o -lm \$(LIB)\n\n"
printf "\${OBJDIR}/linpos.o :\n\t\$(F90) -c \$(FFLAGS) -o \${OBJDIR}/linpos.o ./src_utils/linpos.f\n\n"
	
for infile in $(ls ${SRCDIR}/*.f90)
do
	objName=$(printf "\${OBJDIR}/%s.o" $(basename $infile .f90))
	execName=$(printf "\${BINDIR}/%s" $(basename $infile .f90))
	
	if [ $(basename $infile .f90) = ppsamp ]
	then
		printf "${objName}: ${infile}\n"
		printf "\t\$(F90) -c \$(FFLAGS) -o ${objName} ${infile}\n"
		printf "\n" 
		# this is a particular case
		printf "${execName} :\t\$(OBJDIR)/ppsamp.o \$(OBJDIR)/linpos.o\n\t\$(F90) \$(LDLAGS) -o ${execName} ${objName} \$(OBJDIR)/defn.o \$(OBJDIR)/inter.o \$(OBJDIR)/linpos.o -lm \$(LIB)\n\n"
		
	else
	
		printf "${objName}: ${infile}\n"
		printf "\t\$(F90) -c \$(FFLAGS) -o ${objName} ${infile}\n"
		printf "\n"
	
		printf "${execName}: ${objName}\n"
		printf "\t\$(F90) \$(LDLAGS) -o ${execName} ${objName} \$(OBJDIR)/defn.o \$(OBJDIR)/inter.o -lm \$(LIB)\n\n"
		printf "\n"
	fi
done

}

###################################
###     DOEXECMAKEFILEMAIN     ####
###################################
function doExecMakeFileMain(){
SRCDIR=./src_utils

if [ ! -d ${SRCDIR} ]
then
    printf "${SRCDIR} no such file or directory, abort"
    exit 1
fi

if [ 1 -eq $# ]
then 
	makefileName=$1
else
	makefileName="makefile"
fi
if [ ! -d ./obj ]
then
	mkdir ./obj
fi

doExecMakeFile > ${makefileName}  
}

###################################
###      DOLIBMAKEFILE         ####
###################################
function doLibMakeFile() {
echo "FC_O  = ifort"
echo "FC_SO = ifort"

echo "FFLAGS_O  = -O0 -m64 -fPIC"
echo "FFLAGS_SO = -shared"

printf "\n"
printf "\n"

echo "LIB= libgwlib.so"

printf "OBJDIR=./obj\n\n"
printf "SRC="
for infile in $(ls ${SRCDIR}/*.f90)
do
	printf "\t%s" $infile
done

for infile in $(ls ${SRCDIR}/*.f)
do
	printf "\t%s" $infile
done

printf "\n"
printf "OBJ= "
for infile in $(ls ${SRCDIR}/*.f90)
do
	printf "\t\${OBJDIR}/%s.o" $(basename $infile .f90)
done

for infile in $(ls ${SRCDIR}/*.f)
do
	printf "\t\${OBJDIR}/%s.o" $(basename $infile .f)
done

printf "\n"
printf "\n"
echo "all: \${LIB}"
printf "\n"

printf "\${LIB}: \${OBJ}\n\t\${FC_SO} \${FFLAGS_SO} -o \${LIB} \${OBJ}\n"
printf "\n"


for infile in $(ls ${SRCDIR}/*.f90)
do
	objName=$(printf "\${OBJDIR}/%s.o" $(basename $infile .f90))

	printf "${objName}: ${infile}\n"
	printf "\t\${FC_O} -c \${FFLAGS_O} -o ${objName} ${infile}\n"
	printf "\n"
done

for infile in $(ls ${SRCDIR}/*.f)
do
	objName=$(printf "\${OBJDIR}/%s.o" $(basename $infile .f))

	printf "${objName}: ${infile}\n"
	printf "\t\${FC_O} -c \${FFLAGS_O} -o ${objName} ${infile}\n"
	printf "\n"
done

}

###################################
###     DOLIBMAKEFILEMAIN      ####
###################################
function doLibMakeFileMain(){
SRCDIR=./src_lib

if [ ! -d ${SRCDIR} ]
then
    printf "${SRCDIR} no such file or directory, abort"
    exit 1
fi

if [ 1 -eq $# ]
then 
	makefileName=$1
else
	makefileName="makefile"
fi
if [ ! -d ./obj ]
then
	mkdir ./obj
fi

doLibMakeFile > ${makefileName}
}


###################################
###            MAIN            ####
###################################

F90=ifort
currDir=$(pwd)

if [ ! -d ./lib ]
then
	mkdir lib
fi
if [ ! -d ./bin ]
then
	mkdir ./bin
fi
if [ ! -d ./obj ]
then
	mkdir ./obj
fi

# first compile dependencies: defn.f90 and inter.f90 
$F90 -c ./src_dep/defn.f90
$F90 -c  ./src_dep/inter.f90

cp *.o ./obj
cp *.mod ./obj
cp *.mod ./bin

# creating the shared library
doLibMakeFileMain library.mak
make -f library.mak

# creating the executables
doExecMakeFileMain exec.mak
make -f exec.mak all

# wipe out the mess
mv libgwlib.so ./lib
rm *.o *.mod *.mak
