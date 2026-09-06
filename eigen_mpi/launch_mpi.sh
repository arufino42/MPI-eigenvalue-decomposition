#! /bin/bash
ml purge; ml intel ; ml intel-oneapi-mpi; ml intel-oneapi-vtune; ml hdf5/1.12.2-mpi

for p in 1 2 4 8 12 18 24 36 72
do
    n=$(echo "a=e(l(5*10^8*$p)/3);scale=0;b=(a/$p);b*$p" | bc -l)
    echo "n=$n, p=$p"
    # echo "test_var = $test_var"
    # echo "cut = $cut"
    
    echo "Running code"
    srun -A phys-743 -q phys-743 -N 2 -n $p -t 06:00:00  ./eigen_mpi_main $n
done