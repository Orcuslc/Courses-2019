#include <iostream>
#include <cmath>
using namespace std;

int main(int argc, char* argv[]) {
    int n = 5000;
    double pi = 0.0, c = 1.0;

    for(int i = 0; i < n; i++) {
        pi += c/(2*(double)i + 1.0);
        c *= -1.0;
    }
    pi = pi*4.0;

    double error = M_PI - pi;
    cout.precision(17);
    cout << pi << " " << error << endl;
    return 0;
}
