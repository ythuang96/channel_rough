# Makefile to compile the DNS channel code for various architectures
# All architecture-specific variables are defined in <hostname>.mk and included below as needed

# architecture is identified based on the first 8 characters of the host name
hostname := $(shell hostname | cut -b 1-8 -n)
ifeq ($(hostname),eu-login)
    include euler.mk
else ifeq ($(hostname),limoncel)
    include limoncello.mk
else ifeq ($(hostname),richards)
    include richardson.mk
endif
# INCLUDEPATH is defined in architecture-specific makefiles above
FFLAGS = $(INCLUDEPATH) -c -w -O1
LDFLAGS = $(INCLUDEPATH) -w -O1
objects = modules.v7.o \
	  types.o \
	  boundary_planes.o \
	  h5save.o \
	  save_flowfield.o \
	  wall_roughness.o \
	  cfdiff8.v7.o \
	  cftsingle.o \
	  change.allfac.o \
	  cross.o \
	  fou3D.o \
	  laps.v7.o \
	  main.o \
	  rftsingle.o
VPATH = ../src

# delete target if build fails
.DELETE_ON_ERROR :
# executable
channel : $(objects)
	$(FC) $(LDFLAGS) $(objects) -o $@ $(LIBS)
# fixed-form source files
%.o : %.f ctes3D
	$(FC) $(FFLAGS) $< -o $@
# free-form source files
%.o : %.f90 ctes3D
	$(FC) $(FFLAGS) -ffree-form $< -o $@
%.o : %.f03 ctes3D
	$(FC) $(FFLAGS) -ffree-form $< -o $@
# clean compiled files. -rm tells make to continue despite of errors
.PHONY : clean
clean :
	-rm channel *.o *.mod
# print architecture-specific settings
print :
	$(info hostname is $(hostname))
	$(info FC is $(FC))
	$(info INCLUDEPATH is $(INCLUDEPATH))
	$(info LIBS is $(LIBS))
