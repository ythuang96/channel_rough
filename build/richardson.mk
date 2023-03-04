# Makefile variables specific to the Euler cluster
# the following modules are required to compile this code on Euler:
# open_mpi (tested version: 1.6.5)
# szip (tested version: 2.1)
FC = mpif90
INCLUDEPATH = -I ../src -I/home/yh/libraries/hdf5-1.8.12-serial/include
LIBS = -L/home/yh/libraries/hdf5-1.8.12-serial/lib \
/home/yh/libraries/hdf5-1.8.12-serial/lib/libhdf5hl_fortran.a \
/home/yh/libraries/hdf5-1.8.12-serial/lib/libhdf5_hl.a \
/home/yh/libraries/hdf5-1.8.12-serial/lib/libhdf5_fortran.a \
/home/yh/libraries/hdf5-1.8.12-serial/lib/libhdf5.a \
-lz -ldl -lm -Wl,-rpath -Wl,/home/yh/libraries/hdf5-1.8.12-serial/lib


