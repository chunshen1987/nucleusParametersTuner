#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <ctime>
#include <sys/time.h>

#include "Nucleus.h"

using namespace std;

int main()
{
    time_t start, end;
    double cpu_time_used;
    start = clock();

    // init random seed from system time
    timeval a;
    gettimeofday(&a, 0);
    int randomSeed = -1;
    if (randomSeed<0) randomSeed=a.tv_usec; // randomSeed<0 means to use CPU clock
    srand48(randomSeed);

    Nucleus test(208, 6.0, 0.5, false, 0.0, 0.0);
    double x, y, z;
    test.get_nucleon_corrdinate(x, y, z);
    cout << x << "   " << y << "   " << z << endl;

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    cout << "Time elapsed (in seconds): " << cpu_time_used << endl;
    return(0);
}
