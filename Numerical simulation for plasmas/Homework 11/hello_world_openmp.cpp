#include <omp.h>
#include <iostream>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {
    int total_threads, current_thread;

    #pragma omp parallel default(shared) private(total_threads, current_thread)
    {
        total_threads = omp_get_num_threads();
        current_thread = omp_get_thread_num();
        stringstream ss;
        ss << "Hello World! I am thread " << current_thread << " of " << total_threads << " threads.";
        cout << ss.str() << endl;
    }
    return 0;
}
