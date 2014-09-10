#ifndef NUCLEUS_h
#define NUCLEUS_h

#include <iostream>
#include <cstdlib>
#include <cmath>

class Nucleus
{
    private:
        int atomic_num;    // atomic number for the nucleus
        double atomic_mass;    // atomic mass for the nucleus

        // Woods-Saxon (WS) distribution parameters
        double rho_0;      // WS parameter: central number density [1/fm^3]
        double r_0_std;    // WS parameter: standard nucleus radius [fm]
        double xsi_std;    // WS parameter: standard surface thickness [fm]
        double r_0;        // WS parameter: fit nucleus radius [fm]
        double xsi;        // WS parameter: fit surface thickness [fm]

        bool deformed;     // deformation flag
        double beta2_std, beta4_std; //standard deformation parameters
        double beta2, beta4; //deformation parameters to fit

        double ws_max;     // the maximum density in WS distribution

        double r_max;      // the maximum radius for MC-sampling

    public:
        Nucleus(int atomic_num_in, double r_0_in, double xsi_in, 
                bool deformed_in, double beta2_in, double beta4_in);
        ~Nucleus();

        // main function: returns a random coordinate of nucleon
        double woods_saxon_distribution(double r, double phi, double cos_theta);
        void get_nucleon_corrdinate(double& x, double& y, double& z);
        double spherical_harmonics(int l, double cos_theta);

};

#endif
