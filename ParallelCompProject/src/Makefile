build:
	mpicxx -fopenmp -c ParallelCompProject.c -o ParallelCompProject.o
	mpicxx -fopenmp -o parallelMpiOpenMP ParallelCompProject.o  /usr/local/cuda-9.1/lib64/libcudart_static.a -ldl -lrt

clean:
	rm -f *.o ./parallelMpiOpenMP

runIn1:
	mpiexec -np 4 ./parallelMpiOpenMP

runIn2:
	mpiexec -np 4 -machinefile  mf  -map-by  node  ./parallelMpiOpenMP
