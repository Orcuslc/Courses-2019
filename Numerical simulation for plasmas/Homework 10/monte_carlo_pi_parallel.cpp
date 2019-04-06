#include "mpi.h"
#include <iostream>
using namespace std;

int N = 1e8;

double random(double MIN, double MAX) {
    double r = (double)rand() / RAND_MAX;
    return MIN + (MAX-MIN)*r;
}

inline bool in_circle(double x, double y) {
    return x*x + y*y <= 1.0;
}

void get_input(int argc, char* argv[], int rank, int* per_processor) {
    // for root node
    if(rank == 0) {
        int N_processor = atoi(argv[1]);
        *per_processor = N/(N_processor);
        cout << *per_processor << endl;
    }

    
    // broadcast per_processor from root to other nodes
    MPI_Bcast(per_processor, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

int main(int argc, char** argv) {
    int count = 0, per_processor;
    
    // init
    MPI_Init(&argc, &argv);

    // get number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    // get rank of process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // get input information
    get_input(argc, argv, world_rank, &per_processor);

    srand(time(NULL)+world_rank);

    // the node with world rank == 0 is regarded as the master node
    if(world_rank != 0) {
        for(int i = 0; i < per_processor; i++) { 

            if(i % 1000000 == 0) {
                cout << world_rank << endl;
                cout << i << endl;
            }
            double x = random(-1., 1.);
            double y = random(-1., 1.);
            if(in_circle(x, y))
                count++;
        }
    }
    
    // for master node
    int total_count, total_N;

    // reduce the computed result to root node
    MPI_Reduce(&count, &total_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&per_processor, &total_N, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    cout << total_count << " " << total_N << " " << (double)total_count/(double)total_N << endl;
    // operation at root node
    if(world_rank == 0) {
        double res = (double)total_count/(double)total_N*4.0;
        cout << res << endl;
    }

    MPI_Finalize();
    return 0;
}
