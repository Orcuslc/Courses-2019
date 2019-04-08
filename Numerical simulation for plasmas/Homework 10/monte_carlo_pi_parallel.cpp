#include "mpi.h"
#include <cstdlib>
#include <iostream>
#include <time.h>
using namespace std;

// int N = 1e8;


double random(double MIN, double MAX) {
    double r = (double)rand() / RAND_MAX;
    return MIN + (MAX-MIN)*r;
}

inline bool in_circle(double x, double y) {
    return x*x + y*y <= 1.0;
}
/*
void get_input(int argc, char* argv[], int rank, int* per_processor) {
    // for root node
    if(rank == 0) {
        int N_processor = atoi(argv[1]);
        *per_processor = N/(N_processor);
    }
    
    // broadcast per_processor from root to other nodes
    MPI_Bcast(per_processor, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
*/

int main(int argc, char** argv) {
    int N = 1e8;
    int count = 0, per_processor, total_count;
    double start, cost, elapsed;
    
    // init
    MPI_Init(&argc, &argv);

    // get number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    per_processor = N/world_size;
    N = per_processor*world_size;

    // get rank of process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // get input information
 //   get_input(argc, argv, world_rank, &per_processor);

    srand(time(NULL)+world_rank);

    // timed
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    for(int i = 0; i < per_processor; i++) { 
        double x = random(-1., 1.);
        double y = random(-1., 1.);
        if(in_circle(x, y))
            count++;
    }
    cost = MPI_Wtime() - start;

    // reduce time 
    MPI_Reduce(&cost, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    // reduce the computed result to root node
    MPI_Reduce(&count, &total_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // operation at root node
    if(world_rank == 0) {
        double res = (double)total_count/(double)N*4.0;
        cout.precision(17);
        cout << res << endl;
        cout << elapsed << endl;
    }

    MPI_Finalize();
    return 0;
}
