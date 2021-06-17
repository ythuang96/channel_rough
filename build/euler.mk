# Makefile variables specific to the Euler cluster
# the following modules are required to compile this code on Euler:
# open_mpi (tested version: 1.6.5)
# szip (tested version: 2.1)
FC = mpif90
INCLUDEPATH = -I/cluster/apps/hdf5/1.8.12/x86_64/gcc_4.8.2/openmpi_1.6.5/include
LIBS = -L/cluster/apps/hdf5/1.8.12/x86_64/gcc_4.8.2/openmpi_1.6.5/lib /cluster/apps/hdf5/1.8.12/x86_64/gcc_4.8.2/openmpi_1.6.5/lib/libhdf5hl_fortran.a /cluster/apps/hdf5/1.8.12/x86_64/gcc_4.8.2/openmpi_1.6.5/lib/libhdf5_hl.a /cluster/apps/hdf5/1.8.12/x86_64/gcc_4.8.2/openmpi_1.6.5/lib/libhdf5_fortran.a /cluster/apps/hdf5/1.8.12/x86_64/gcc_4.8.2/openmpi_1.6.5/lib/libhdf5.a -lsz -lz -ldl -lm -Wl,-rpath -Wl,/cluster/apps/hdf5/1.8.12/x86_64/gcc_4.8.2/openmpi_1.6.5/lib

