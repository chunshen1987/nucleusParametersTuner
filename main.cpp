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
#include "ParameterReader.h"

using namespace std;

int main(int argc, char *argv[])
{
    time_t start, end;
    double cpu_time_used;
    start = clock();

    // Read-in parameters
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.readFromArguments(argc, argv);
    paraRdr.echo();

    // init random seed from system time
    timeval a;
    gettimeofday(&a, 0);
    int randomSeed = paraRdr.getVal("randomSeed");
    if (randomSeed<0) randomSeed=a.tv_usec; // randomSeed<0 means to use CPU clock
    srand48(randomSeed);

    Nucleus test(&paraRdr);
    double x, y, z;
    test.get_nucleon_corrdinate(x, y, z);
    cout << x << "   " << y << "   " << z << endl;

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    cout << "Time elapsed (in seconds): " << cpu_time_used << endl;
    return(0);
}
