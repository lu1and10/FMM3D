# Makefile for FMM3D
#
# This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 

# compiler, and linking from C, fortran

CC=gcc
FC=gfortran

FFLAGS=-fPIC -O3 -funroll-loops -march=native 
CFLAGS= -std=c99 
CFLAGS+= $(FFLAGS) 

CLINK = -L/usr/local/Cellar/gcc/8.3.0/lib/gcc/8 -lgfortran -lm

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS = -fopenmp
MOMPFLAGS = -lgomp -D_OPENMP

# flags for MATLAB MEX compilation..
MFLAGS=-largeArrayDims -DMWF77_UNDERSCORE1 -lgfortran -lgomp -D_OPENMP -lm
MWFLAGS=-c99complex 

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../mwrap-0.33/mwrap

# For your OS, override the above by placing make variables in make.inc
-include make.inc

# multi-threaded libs & flags needed
ifeq ($(OMP),on)
CFLAGS += $(OMPFLAGS)
FFLAGS += $(OMPFLAGS)
endif

LIBNAME=libfmm3d
DYNAMICLIB = lib/$(LIBNAME).so
STATICLIB = lib-static/$(LIBNAME).a

# objects to compile
#
# Common objects
COM = src/Common
COMOBJS = $(COM)/besseljs3d.o $(COM)/cdjseval3d.o $(COM)/dfft.o \
	$(COM)/fmmcommon.o $(COM)/legeexps.o $(COM)/prini.o \
	$(COM)/rotgen.o $(COM)/rotproj.o $(COM)/rotviarecur.o \
	$(COM)/tree_lr_3d.o $(COM)/yrecursion.o 
 
# Helmholtz objects
HELM = src/Helmholtz
HOBJS = $(HELM)/h3dcommon.o $(HELM)/h3dterms.o $(HELM)/h3dtrans.o \
	$(HELM)/helmrouts3d.o $(HELM)/hfmm3d.o $(HELM)/hfmm3dwrap.o \
	$(HELM)/hfmm3dwrap_vec.o $(HELM)/hpwrouts.o \
	$(HELM)/hwts3.o $(HELM)/numphysfour.o $(HELM)/projections.o \
	$(HELM)/quadread.o

# Laplace objects
LAP = src/Laplace
LOBJS = $(LAP)/lwtsexp_sep1.o $(LAP)/l3dterms.o $(LAP)/l3dtrans.o \
	$(LAP)/laprouts3d.o $(LAP)/lfmm3d.o $(LAP)/lfmm3dwrap.o \
	$(LAP)/lfmm3dwrap_vec.o $(LAP)/lwtsexp_sep2.o \
	$(LAP)/lpwrouts.o

# Test objects
TOBJS = $(COM)/hkrand.o $(COM)/dlaran.o

# C Headers and objects
COBJS = c/cprini.o c/utils.o
CHEADERS = c/cprini.h c/utils.h c/hfmm3d_c.h

OBJS = $(COMOBJS) $(HOBJS) $(LOBJS)

.PHONY: usage lib examples test perftest python all c c-examples

default: usage

all: lib examples test perftest python c c-examples

usage:
	@echo "Makefile for FMM3D. Specify what to make:"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make examples - compile and run fortran in examples/"
	@echo "  make c-examples - compile and run fortran in c/"
	@echo "  make test - compile and run quick math validation tests"
	@echo "  make matlab - compile matlab interfaces"
	@echo "  make python - compile and test python interfaces"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=OFF' for multi-threaded (otherwise single-threaded)"


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

# build the library...
lib: $(STATICLIB) $(DYNAMICLIB)
ifeq ($(OMP),ON)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif
$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)
$(DYNAMICLIB): $(OBJS) 
	$(FC) -shared $(OMPFLAGS) $(OBJS) -o $(DYNAMICLIB)

# matlab..
MWRAPFILE = fmm3d
GATEWAY = $(MWRAPFILE)
matlab:	$(STATICLIB) matlab/$(GATEWAY).c
	$(MEX) matlab/$(GATEWAY).c $(STATICLIB) $(MFLAGS) -output matlab/fmm3d


mex:  $(STATICLIB)
	cd matlab; $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw;\
	$(MEX) $(GATEWAY).c ../$(STATICLIB) $(MFLAGS) -output fmm3d

#python
python:
	cd python && pip install -e . && cd test && pytest -s

# testing routines
#
test: $(STATICLIB) $(TOBJS) test/helmrouts test/hfmm3d test/hfmm3d_vec test/laprouts test/lfmm3d test/lfmm3d_vec
	(cd test/Helmholtz; ./run_helmtest.sh)
	(cd test/Laplace; ./run_laptest.sh)
	cat print_testreshelm.txt
	cat print_testreslap.txt
	rm print_testreshelm.txt
	rm print_testreslap.txt

test/helmrouts: 
	$(FC) $(FFLAGS) test/Helmholtz/test_helmrouts3d.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/test_helmrouts3d 

test/hfmm3d:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/test_hfmm3d

test/hfmm3d_vec:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d_vec.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/test_hfmm3d_vec 

test/laprouts:
	$(FC) $(FFLAGS) test/Laplace/test_laprouts3d.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/Laplace/test_laprouts3d 

test/lfmm3d:
	$(FC) $(FFLAGS) test/Laplace/test_lfmm3d.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/Laplace/test_lfmm3d

test/lfmm3d_vec:
	$(FC) $(FFLAGS) test/Laplace/test_lfmm3d_vec.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/Laplace/test_lfmm3d_vec 



#
##  examples
#

examples: $(STATICLIB) examples/ex1_helm examples/ex2_helm examples/ex1_lap examples/ex2_lap
	time -p ./examples/example1_lap
	time -p ./examples/example2_lap
	time -p ./examples/example1_helm
	time -p ./examples/example2_helm

examples/ex1_lap:
	$(FC) $(FFLAGS) examples/lfmm3d_example.f $(TOBJS) $(COMOBJS) $(LOBJS) -o examples/example1_lap

examples/ex2_lap:
	$(FC) $(FFLAGS) examples/lfmm3d_vec_example.f $(TOBJS) $(COMOBJS) $(LOBJS) -o examples/example2_lap

examples/ex1_helm:
	$(FC) $(FFLAGS) examples/hfmm3d_example.f $(TOBJS) $(COMOBJS) $(HOBJS) -o examples/example1_helm

examples/ex2_helm:
	$(FC) $(FFLAGS) examples/hfmm3d_vec_example.f $(TOBJS) $(COMOBJS) $(HOBJS) -o examples/example2_helm



# C interface
c: $(COBJS) $(OBJS) $(CHEADERS) c/lfmm3d c/hfmm3d

c/lfmm3d:
	$(CC) $(CFLAGS) c/test_lfmm3d.c $(COBJS) $(OBJS) $(CLINK) -o c/test_lfmm3d
	time -p c/test_lfmm3d

c/hfmm3d:
	$(CC) $(CFLAGS) c/test_hfmm3d.c $(COBJS) $(OBJS) $(CLINK) -o c/test_hfmm3d
	time -p c/test_hfmm3d



# C examples
c-examples: $(COBJS) $(OBJS) $(CHEADERS) c/ex1_lap c/ex2_lap c/ex1_helm c/ex2_helm 
	time -p c/example1_lap
	time -p c/example2_lap
	time -p c/example1_helm
	time -p c/example2_helm

c/ex1_lap:
	$(CC) $(CFLAGS) c/lfmm3d_example.c $(COBJS) $(OBJS) $(CLINK) -o c/example1_lap

c/ex2_lap:
	$(CC) $(CFLAGS) c/lfmm3d_vec_example.c $(COBJS) $(OBJS) $(CLINK) -o c/example2_lap

c/ex1_helm:
	$(CC) $(CFLAGS) c/hfmm3d_example.c $(COBJS) $(OBJS) $(CLINK) -o c/example1_helm

c/ex2_helm:
	$(CC) $(CFLAGS) c/hfmm3d_vec_example.c $(COBJS) $(OBJS) $(CLINK) -o c/example2_helm

clean: objclean
	rm -f lib-static/*.a lib/*.so
	rm -f python/*.so
	rm -rf python/build
	rm -rf fmm3dpy.egg-info
	rm -f examples/example1_helm
	rm -f examples/example2_helm
	rm -f examples/example1_lap
	rm -f examples/example2_lap
	rm -f c/example1_helm
	rm -f c/example2_helm
	rm -f c/example1_lap
	rm -f c/example2_lap
	rm -f c/test_hfmm3d
	rm -f c/test_lfmm3d
	
	

objclean: 
	rm -f $(OBJS) $(COBJS) $(TOBJS)
	rm -f test/*.o examples/*.o
