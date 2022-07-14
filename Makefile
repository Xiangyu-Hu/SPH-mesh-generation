################################################################################
# Build script for MRSPH program
# Writen by Zhe Ji
# Date: 23/03/2018
################################################################################
include Makefile.in
#include Makefile.eu.in

ifeq ($(MACHINE), LRZ)
TET_DIR=/home/hpc/pr53vu/ga86wux3/tetgen1.5.1
#TET_DIR=/home/hpc/t7841/ga86wux2/tetgen1.5.1
endif

DIRS := $(shell find . -maxdepth 3 -type d)
rmtrashfile=rm -rf $(foreach d,$(DIRS), $d/*~)

# name of the executable
EXECUTABLE	:= MRSPH 

# compiler
ifeq ($(MACHINE), LRZ)
  CXX := mpiCC
else
  #CXX := /scratch/prefix-mpich-icpc/bin/mpicxx
  CXX := mpicxx
endif

# includes
INCLUDES	:= -I$(SRC) -I$(SRC)/$(UTILS) -I$(SRC)/$(DIR_SOLVER) -I$(SRC)/$(DIR_BASE)

# source files
CCFILES	:= $(wildcard $(SRC)/*.cpp)
CCFILES += $(wildcard $(SRC)/$(UTILS)/*.cpp)
CCFILES += $(wildcard $(SRC)/$(DIR_SOLVER)/*.cpp)
CCFILES += $(wildcard $(SRC)/$(DIR_BASE)/*.cpp)
ifeq ($(PARALLEL), MULTICPU)
  CCFILES	+= $(wildcard $(SRC)/$(DIR_DD)/*.cpp)
  CCFILES	+= $(wildcard $(SRC)/$(DIR_GRAPH)/*.cpp)
  INCLUDES += -I$(SRC)/$(DIR_DD) -I$(SRC)/$(DIR_GRAPH)
endif
ifeq ($(SOLVER), mesh_generation)
  CCFILES += $(wildcard $(SRC)/$(DIR_LEVEL_SET)/*.cpp)
  CCFILES += $(wildcard $(SRC)/$(DIR_BOUNDARY)/*.cpp)
  CCFILES += $(wildcard $(SRC)/$(DIR_BASE_MESH)/*.cpp)
  CCFILES += $(wildcard $(SRC)/$(DIR_GEOMETRY)/*.cpp)
  INCLUDES += -I$(SRC)/$(DIR_LEVEL_SET)
  INCLUDES += -I$(SRC)/$(DIR_BOUNDARY)
  INCLUDES += -I$(SRC)/$(DIR_BASE_MESH)
  INCLUDES += -I$(SRC)/$(DIR_GEOMETRY)
  INCLUDES += -I${TET_DIR}
endif
ifeq ($(SOLVER), eularian_cat)
  CCFILES += $(wildcard $(SRC)/$(DIR_BOUNDARY)/*.cpp)
  INCLUDES += -I$(SRC)/$(DIR_BOUNDARY)
endif

# compiler flags
CCFLAGS   := -std=c++0x -m64

# libraries
ifeq ($(MACHINE), LRZ)
  BOOST_PREFIX56=${HOME}/prefix-boost-156
  BOOST_PREFIX =${HOME}/prefix-boost-154
  VORO_PREFIX = ${HOME}/voro++
  BOOST_MPI_INCLUDES := -I${BOOST_PREFIX}/include -I${BOOST_PREFIX56}/include -I${VORO_PREFIX}/include
  BOOST_MPI_FLAGS   := -L${BOOST_PREFIX}/lib -lboost_mpi -lboost_serialization -Wl,-rpath,${BOOST_PREFIX}/lib
  BOOST_MPI_FLAGS   += -L${VORO_PREFIX}/lib -lvoro++ -Wl,-rpath,${VORO_PREFIX}/lib
  TBB_PREFIX=${TBB_BASE}
  TBB_INCLUDE:=-I${TBB_PREFIX}/include
  TBB_MPI_FLAGS:=${TBB_SHLIB}
else
  # older version
  # TBB_INCLUDE:=-I${TBB_PREFIX}/include
  # TBB_MPI_FLAGS:=${TBB_SHLIB}
  # BOOST_PREFIX56=/scratch/prefix-boost-1.56
  # BOOST_PREFIX = /scratch/prefix-boost-1.54
  # VORO_PREFIX = /scratch/voro++
  # BOOST_MPI_INCLUDES := -I${BOOST_PREFIX}/include -I${BOOST_PREFIX56}/include -I${VORO_PREFIX}/include
  # BOOST_MPI_FLAGS   := -L${BOOST_PREFIX}/lib -lboost_mpi -lboost_serialization -Wl,-rpath,${BOOST_PREFIX}/lib
  # BOOST_MPI_FLAGS   += -L${VORO_PREFIX}/lib -lvoro++ -Wl,-rpath,${VORO_PREFIX}/lib

  TBB_INCLUDE:=-I${TBB_DIR}/include/tbb
  TBB_MPI_FLAGS:=-L${TBB_DIR}/lib/intel64/gcc4.7 -ltbb -ltbbmalloc -Wl,-rpath,${TBB_DIR}/lib/intel64/gcc4.7
  BOOST_PREFIX:=${BOOST_DIR}
  VORO_PREFIX=${VORO_DIR}
  BOOST_MPI_INCLUDES := -I${BOOST_PREFIX}/ -I${VORO_PREFIX}/include/voro++
  BOOST_MPI_FLAGS   := -L${BOOST_PREFIX}/stage/lib -lboost_mpi -lboost_serialization -Wl,-rpath,${BOOST_PREFIX}/stage/lib
  BOOST_MPI_FLAGS   += -L${VORO_PREFIX}/lib -lvoro++ -Wl,-rpath,${VORO_PREFIX}/lib
endif

# debug build flags
ifeq ($(dbg),1)
	CCFLAGS   += -g -O0
	LIBPARALLEL = -ltbb_debug -ltbbmalloc_debug
else
	CCFLAGS   += -O3
	LIBPARALLEL = -ltbb -ltbbmalloc
endif

# warning flags
# GCCWARN_FLAGS := -W -Wall -Wformat -Wparentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wreturn-type -Wno-unused-function 
GCCWARN_FLAGS :=

# PREPARE AND COMPILE:
CXXFLAGS = $(CCFLAGS) $(DFLAGS) $(GCCWARN_FLAGS)

# common includes and LD_flags
INCLUDES      += $(BOOST_MPI_INCLUDES) ${TBB_INCLUDE}
LDFLAGS       := $(LIBPARALLEL) -lrt $(LIBS) ${BOOST_MPI_FLAGS} ${TBB_MPI_FLAGS}
ifeq ($(SOLVER), mesh_generation)
LDFLAGS       += -L${TET_DIR} -ltet
endif
# Location of output and obj files
OUTDATADIR	?= ./outdata
TPOUTDATADIR	?= ./tpout
OBJECTDIR 	?= ./obj
RESTARTDIR	?= ./rstfile
TESTDIR         ?= ./test_log
ifeq ($(SOLVER), mesh_generation)
RESTARTPTLDIR   ?= ./rstfile_particles
endif

# Set up object files
OBJECTS = $(patsubst %.cpp,$(OBJECTDIR)/%.o, $(notdir $(CCFILES)))

# Target rules
all: directories $(EXECUTABLE)

# Link commands:
$(OBJECTDIR)/%.o : $(SRC)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

$(OBJECTDIR)/%.o : $(SRC)/$(UTILS)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

$(OBJECTDIR)/%.o : $(SRC)/$(DIR_SOLVER)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

$(OBJECTDIR)/%.o : $(SRC)/$(DIR_BASE)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

ifeq ($(SOLVER), mesh_generation)
  $(OBJECTDIR)/%.o : $(SRC)/$(DIR_LEVEL_SET)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<
  $(OBJECTDIR)/%.o : $(SRC)/$(DIR_BOUNDARY)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<
  $(OBJECTDIR)/%.o : $(SRC)/$(DIR_BASE_MESH)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<
  $(OBJECTDIR)/%.o : $(SRC)/$(DIR_GEOMETRY)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<
endif
ifeq ($(SOLVER), eularian_cat)
  $(OBJECTDIR)/%.o : $(SRC)/$(DIR_BOUNDARY)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<
endif
ifeq ($(PARALLEL), MULTICPU)
  $(OBJECTDIR)/%.o : $(SRC)/$(DIR_DD)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<
  $(OBJECTDIR)/%.o : $(SRC)/$(DIR_GRAPH)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<
endif

$(EXECUTABLE): directories $(OBJECTS) Makefile
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LDFLAGS)

directories:
	mkdir -p $(OBJECTDIR)
	mkdir -p $(TPOUTDATADIR)
	mkdir -p $(OUTDATADIR)
	mkdir -p $(RESTARTDIR)
	mkdir -p $(TESTDIR)
ifeq ($(SOLVER), mesh_generation)
	mkdir -p $(RESTARTPTLDIR)
endif
run: all
	mpiexec -np $(NP) ./MRSPH

clean:
	rm -rf $(EXECUTABLE) $(OBJECTDIR)/* *.bin
	$(call rmtrashfile,$(SRC)/)

cleandata:
	rm -rf $(OUTDATADIR)/* $(TPOUTDATADIR)/* $(TESTDIR)/*
	$(call rmtrashfile,$(SRC)/)

reset:
	rm -rf $(EXECUTABLE) $(OBJECTDIR)/* *.bin $(OBJECTDIR) $(OUTDATADIR) $(TPOUTDATADIR) $(TESTDIR) 
	$(call rmtrashfile,$(SRC)/)
