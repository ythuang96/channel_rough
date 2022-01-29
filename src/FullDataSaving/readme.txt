These two files saves intermediate data from the DNS computation.

Data include: u,v,w,o1,o2,o3,hv,hg in Fourier space from the current time step;

dvdy, do2dy in Fourier space at rk1;
u,v,w,o1,o2,o3 in physical space at rk1;
H1,H2,H3 in Fourier space at rk1;

The data saving at rk1 can be adjusted to save data at any rk step in cross.f