rm test_grid_*
mpirun -np 4 BEC_ground_par
g++ data_merge_gs.c -o merge.out
./merge.out