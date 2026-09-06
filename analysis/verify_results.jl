using LinearAlgebra
using HDF5
using Glob
vn=[144  216  288  576  648  1152  1944  2304  4608  5832]
for n in vn
    @show n
    mat_files=glob("initial_matrix_mpi_n=$(n)*","project_eigen/data")
    eigenvector_files=glob("eigenvectors_mpi_n=$(n)*","project_eigen/data")
    eigenvalue_files=glob("eigenvalues_mpi_n=$(n)*","project_eigen/data")
    f=h5open(mat_files[1],"r")
    A=reshape(read(f,"vector"),n,n)
    f=h5open(eigenvector_files[1],"r")
    Z=reshape(read(f,"vector"),n,n)
    f=h5open(eigenvalue_files[1],"r")
    Λ=read(f,"vector")
    @show maximum(abs.(sort(Λ) .-eigvals(A)) )
    @show maximum(abs.(sort(diag(Z*A*Z')) .-eigvals(A)) )
    @show maximum(abs.(Z*A*Z'-diagm(diag(Z*A*Z'))))
end
# begin
#     open("output/initial_matrix_mpi_n=1600_0.txt", "r") do io
#         global A=readdlm(io,',',skipblanks=true)
#     end 
#     A=[Float64(A[i,j]) for i in 1:200, j in 1:200]
#     open("output/eigenvalues_mpi_n=1600_0.txt", "r") do io
#         global eig=readdlm(io,',',skipblanks=true)
#     end 
#     @show maximum(sort(eig[:]).-eigvals(A)) 
# end
