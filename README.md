# Parallel eigenvalue decomposition.
* Serial code is in folder `eigen_serial`.
* MPI code is in folder `eigen_mpi`. Matrix size must be chosen as a multiple of the number of processors. Random matrices may be used, or a matrix may be loaded using HDF5.
* Auxiliary Julia files for scaling analysis are included in folder `analysis`.
* Some sample data is included in `data` and `data_weak_scaling` (more would occupy too much disk space).