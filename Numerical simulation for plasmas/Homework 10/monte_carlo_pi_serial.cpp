
#include <iostream>
using namespace std;

double random(double MIN, double MAX) {
    double r = (double)rand() / RAND_MAX;
    return MIN + (MAX-MIN)*r;
}

bool in_circle(double x, double y) {
    return x*x + y*y <= 1.0;
}

int main() {
    int N = 1e8;
    
    int sum = 0.;
    for(int i = 0; i < N; i++) {
        double x = random(-1., 1.);
        double y = random(-1., 1.);
        if(in_circle(x, y)) 
            sum++;
    }
    double res = 4.*(double)sum/(double)N;
    cout.precision(17);
    cout << res << endl;
    return 0;
}
