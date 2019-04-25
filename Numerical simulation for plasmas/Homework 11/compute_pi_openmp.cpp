#include <iostream>
#include <cmath>
using namespace std;

int main(int argc, char* argv[]) {
    int total_threads;
    double pi, sum;
    int N = 1e9;

    #pragma omp parallel for reduction (+:sum)
    for(int i = 0; i < N; i++) {
        sum += pow(-1, i)/(2.*i+1.);
    }

    #pragma omp master
    {
        pi = 4*sum;
        cout.precision(17);
        cout << pi << endl;
    }
    return 0;
}
