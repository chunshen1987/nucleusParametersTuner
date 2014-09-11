#ifndef MAKEDENSITY_h
#define MAKEDENSITY_h

#include "Nucleus.h"
#include "ParameterReader.h"
#include <gsl/gsl_vector.h> 

class MakeDensity;

struct CCallbackHolder   // callback structure for gsl
{
   MakeDensity* clsPtr;
   void *params;
};
      
class MakeDensity
{
    private:
        ParameterReader *paraRdr;
        double gauss_nucl_width;          // gaussian width of nucleon size
        double siginNN;                   // inelastic nucleon-nucleon cross section

        int n_r, n_phi, n_theta;
        double r_max;                     // maximum radius in nucleon density distribution
        double*** rho;                    // 3-d nucleon density distribution
        double *rho_r, *rho_r_std;        // 1-d nucleon density distribution along radial direction
        double *r, *r_weight;
        double *phi, *phi_weight;
        double *sin_phi, *cos_phi;
        double *cos_theta, *cos_theta_weight;
        double *sin_theta;
        
        int n_event;                      // number of events used to calculate event averaged nucleon density

        Nucleus* test_nucleus;

    public:
        MakeDensity(ParameterReader*);
        ~MakeDensity();
        
        double calculate_density(double ws_r0, double ws_xsi);
        double calculate_chisq_rho_r(const gsl_vector *v, void *params);
        double minimize_chisq();
        static double CCallback_chisq_rho_r(const gsl_vector *x, void* params)
        {
            CCallbackHolder* h = static_cast<CCallbackHolder*>(params);
            return h->clsPtr->calculate_chisq_rho_r(x, h->params);
        }

};

#endif
