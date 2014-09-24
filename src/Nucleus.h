#ifndef NUCLEUS_h
#define NUCLEUS_h

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "ParameterReader.h"
#include "Particle.h"

class Nucleus
{
    private:
        ParameterReader* paraRdr;
        int atomic_num;    // atomic number for the nucleus
        double atomic_mass;    // atomic mass for the nucleus

        double r_min;      // minimum separation between nucleons

        // Woods-Saxon (WS) distribution parameters
        double rho_0;      // WS parameter: central number density [1/fm^3]
        double r_0_std;    // WS parameter: standard nucleus radius [fm]
        double xsi_std;    // WS parameter: standard surface thickness [fm]
        double r_0;        // WS parameter: fit nucleus radius [fm]
        double xsi;        // WS parameter: fit surface thickness [fm]

        int deformed;     // deformation flag
        double beta2_std, beta4_std; //standard deformation parameters
        double beta2, beta4; //deformation parameters to fit

        double ws_max;     // the maximum density in WS distribution

        double r_max;      // the maximum radius for MC-sampling

        vector<Particle*> nucleus;

    public:
        Nucleus(ParameterReader* paraRdr_in);
        ~Nucleus();

        void set_woods_saxon_parameters(double r0_in, double xsi_in);
        void set_woods_saxon_parameters(double r0_in, double xsi_in, double beta2_in, double beta4_in);
        double woods_saxon_distribution(double r, double phi, double cos_theta);
        double standard_woods_saxon_distribution(double r, double phi, double cos_theta);
        void sample_nucleon_corrdinate(double& x, double& y, double& z);
        void generate_nucleus();
        double spherical_harmonics(int l, double cos_theta);

        int get_atomic_num(){return(atomic_num);};
        void get_nucleon_position(int index, double& x, double& y, double& z);
        double get_ws_r0_std() {return(r_0_std);};
        double get_ws_xsi_std() {return(xsi_std);};
        double get_ws_beta2_std() {return(beta2_std);};
        double get_ws_beta4_std() {return(beta4_std);};

};

#endif
