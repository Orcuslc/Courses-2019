#include "mpi.h"
#include <iostream>
using namespace std;

int main() {
    // init
    MPI_Init(NULL, NULL);
    
    // get number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // get rank of process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    // make sure all processes run in order via Barrier
    for(int i = 0; i < world_size; i++) {
        if(i == world_rank)
            cout << "Hello World! I am processor " << world_rank << " of " << world_size << " processors." << endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // MPI_Finalize
    MPI_Finalize();
    return 0; 
}
