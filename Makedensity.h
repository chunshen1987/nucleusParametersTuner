#ifndef MAKEDENSITY_h
#define MAKEDENSITY_h

#include "Nucleus.h"
#include "ParameterReader.h"

class MakeDensity
{
    private:
        ParameterReader *paraRdr;
        double gauss_nucl_width;
        double siginNN;

        int n_r, n_phi, n_theta;
        double r_max;
        double*** rho;
        double *rho_r, *rho_r_std;
        double *r, *r_weight;
        double *phi, *phi_weight;
        double *sin_phi, *cos_phi;
        double *cos_theta, *cos_theta_weight;
        double *sin_theta;
        
        int n_event;

        Nucleus* test_nucleus;

    public:
        MakeDensity(ParameterReader*);
        ~MakeDensity();
        
        double calculate_density(double ws_r0, double ws_xsi);
        double calculate_rho_r(double ws_r0, double ws_xsi);
        double compare_rho_r_with_standard_ws();

};

#endif
