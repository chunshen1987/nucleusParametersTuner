#include <iomanip>
#include <fstream>

#include "ParameterReader.h"
#include "Nucleus.h"
#include "MakeDensity.h"
#include "gauss_quadrature.h"

#include <gsl/gsl_vector.h> 
#include <gsl/gsl_multimin.h>

using namespace std;

MakeDensity::MakeDensity(ParameterReader *paraRdr_in)
{
    paraRdr = paraRdr_in;
    n_event = paraRdr->getVal("n_event");
    
    // ref: http://arxiv.org/pdf/1108.5379.pdf
    gauss_nucl_width = 0.473;  // calculated from sigma_nn_inel at sqrt{s} = 23.5 GeV

    r_max = paraRdr->getVal("r_max");
    n_r = paraRdr->getVal("n_r");
    n_phi = paraRdr->getVal("n_phi");
    n_theta = paraRdr->getVal("n_theta");

    r = new double [n_r];
    r_weight = new double [n_r];
    gauss_quadrature(n_r, 1, 0.0, 0.0, 0.0, r_max, r, r_weight);
    phi = new double [n_phi];
    sin_phi = new double [n_phi];
    cos_phi = new double [n_phi];
    phi_weight = new double [n_phi];
    gauss_quadrature(n_phi, 1, 0.0, 0.0, 0.0, 2*M_PI, phi, phi_weight);
    for(int i = 0; i < n_phi; i++)
    {
        sin_phi[i] = sin(phi[i]);
        cos_phi[i] = cos(phi[i]);
    }
    cos_theta = new double [n_theta];
    sin_theta = new double [n_theta];
    cos_theta_weight = new double [n_theta];
    gauss_quadrature(n_theta, 1, 0.0, 0.0, -1.0, 1.0, cos_theta, cos_theta_weight);
    for(int i = 0; i < n_theta; i++)
        sin_theta[i] = sqrt(1. - cos_theta[i]*cos_theta[i]);

    rho = new double** [n_r];
    rho_r = new double [n_r];
    rho_r_std = new double [n_r];
    rho_r_theta = new double* [n_r];
    rho_r_theta_std = new double* [n_r];
    for(int i = 0; i < n_r; i++)
    {
        rho[i] = new double* [n_phi];
        rho_r[i] = 0.0;
        rho_r_std[i] = 0.0;
        rho_r_theta[i] = new double [n_theta];
        rho_r_theta_std[i] = new double [n_theta];
        for(int k = 0; k < n_theta; k++)
        {
            rho_r_theta[i][k] = 0.0;
            rho_r_theta_std[i][k] = 0.0;
        }
        for(int j = 0; j < n_phi; j++)
        {
            rho[i][j] = new double [n_theta];
            for(int k = 0; k < n_theta; k++)
                rho[i][j][k] = 0.0;
        }
    }
    
    test_nucleus = new Nucleus (paraRdr);
}

MakeDensity::~MakeDensity()
{
    for(int i = 0; i < n_r; i++)
    {
        for(int j = 0; j < n_phi; j++)
            delete[] rho[i][j];
        delete[] rho[i];
        delete[] rho_r_theta[i];
        delete[] rho_r_theta_std[i];
    }
    delete[] rho;
    delete[] rho_r;
    delete[] rho_r_std;
    delete[] rho_r_theta;
    delete[] rho_r_theta_std;

    delete[] r;
    delete[] r_weight;
    delete[] phi;
    delete[] phi_weight;
    delete[] sin_phi;
    delete[] cos_phi;
    delete[] sin_theta;
    delete[] cos_theta;
    delete[] cos_theta_weight;
    delete test_nucleus;
}

double MakeDensity::calculate_density(double ws_r0, double ws_xsi)
{
    // mc sample nucleus configuration and compute event averaged nucleon density distribution
    double prefactor = 1./(pow(gauss_nucl_width, 3)*pow(2.*M_PI, 1.5));
    double sigma_sq = 2.*gauss_nucl_width*gauss_nucl_width;
    test_nucleus->set_woods_saxon_parameters(ws_r0, ws_xsi);
    for(int ir = 0; ir < n_r; ir++)
         for(int iphi = 0; iphi < n_phi; iphi++)
             for(int itheta = 0; itheta < n_theta; itheta++)
                 rho[ir][iphi][itheta] = 0.0;
    for(int i = 0; i < n_event; i ++)
    {
        //cout << "processing event " << i << endl;
        test_nucleus->generate_nucleus();
        for(int inucl = 0; inucl < test_nucleus->get_atomic_num(); inucl++)
        {
            double x0, y0, z0;
            test_nucleus->get_nucleon_position(inucl, x0, y0, z0);
            for(int iphi = 0; iphi < n_phi; iphi++)
            {
                double sin_phi_lcoal = sin_phi[iphi];
                double cos_phi_lcoal = cos_phi[iphi];
                for(int itheta = 0; itheta < n_theta; itheta++)
                {
                    double cos_theta_local = cos_theta[itheta];
                    double sin_theta_local = sin_theta[itheta];
                    for(int ir = 0; ir < n_r; ir++)
                    {
                        double r_local = r[ir];
                        double x = r_local*cos_phi_lcoal*sin_theta_local;
                        double y = r_local*sin_phi_lcoal*sin_theta_local;
                        double z = r_local*cos_theta_local;
                        double r_diff_sq = (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0);
                        if(r_diff_sq < 25.*sigma_sq)
                            rho[ir][iphi][itheta] += prefactor*exp(-r_diff_sq/sigma_sq);
                    }
                }
            }
        }
    }
    for(int ir = 0; ir < n_r; ir++)
         for(int iphi = 0; iphi < n_phi; iphi++)
             for(int itheta = 0; itheta < n_theta; itheta++)
                 rho[ir][iphi][itheta] = rho[ir][iphi][itheta]/n_event;
    return(0.0);
}

double MakeDensity::calculate_density(double ws_r0, double ws_xsi, double ws_beta2, double ws_beta4)
{
    // mc sample nucleus configuration and compute event averaged nucleon density distribution
    double prefactor = 1./(pow(gauss_nucl_width, 3)*pow(2.*M_PI, 1.5));
    double sigma_sq = 2.*gauss_nucl_width*gauss_nucl_width;
    test_nucleus->set_woods_saxon_parameters(ws_r0, ws_xsi, ws_beta2, ws_beta4);
    for(int ir = 0; ir < n_r; ir++)
         for(int iphi = 0; iphi < n_phi; iphi++)
             for(int itheta = 0; itheta < n_theta; itheta++)
                 rho[ir][iphi][itheta] = 0.0;
    for(int i = 0; i < n_event; i ++)
    {
        test_nucleus->generate_nucleus();
        for(int inucl = 0; inucl < test_nucleus->get_atomic_num(); inucl++)
        {
            double x0, y0, z0;
            test_nucleus->get_nucleon_position(inucl, x0, y0, z0);
            for(int iphi = 0; iphi < n_phi; iphi++)
            {
                double sin_phi_lcoal = sin_phi[iphi];
                double cos_phi_lcoal = cos_phi[iphi];
                for(int itheta = 0; itheta < n_theta; itheta++)
                {
                    double cos_theta_local = cos_theta[itheta];
                    double sin_theta_local = sin_theta[itheta];
                    for(int ir = 0; ir < n_r; ir++)
                    {
                        double r_local = r[ir];
                        double x = r_local*cos_phi_lcoal*sin_theta_local;
                        double y = r_local*sin_phi_lcoal*sin_theta_local;
                        double z = r_local*cos_theta_local;
                        double r_diff_sq = (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0);
                        if(r_diff_sq < 25.*sigma_sq)
                            rho[ir][iphi][itheta] += prefactor*exp(-r_diff_sq/sigma_sq);
                    }
                }
            }
        }
    }
    for(int ir = 0; ir < n_r; ir++)
         for(int iphi = 0; iphi < n_phi; iphi++)
             for(int itheta = 0; itheta < n_theta; itheta++)
                 rho[ir][iphi][itheta] = rho[ir][iphi][itheta]/n_event;
    return(0.0);
}

double MakeDensity::calculate_chisq_rho_r(const gsl_vector *v, void* params)
{
    // compute the chi square of the nucleon density distribution along radial direction
    double ws_r0, ws_xsi;
    ws_r0 = gsl_vector_get(v, 0);
    ws_xsi = gsl_vector_get(v, 1);

    cout << ws_r0 << "   " << ws_xsi << endl;
    calculate_density(ws_r0, ws_xsi);
    for(int ir = 0; ir < n_r; ir++)
    {
        rho_r[ir] = 0.0;
        rho_r_std[ir] = 0.0;
    }
    for(int ir = 0; ir < n_r; ir++)
        for(int iphi = 0; iphi < n_phi; iphi++)
            for(int itheta = 0; itheta < n_theta; itheta++)
            {
                 rho_r[ir] += rho[ir][iphi][itheta]*phi_weight[iphi]*cos_theta_weight[itheta];
                 rho_r_std[ir] += test_nucleus->standard_woods_saxon_distribution(r[ir], phi[iphi], cos_theta[itheta])*phi_weight[iphi]*cos_theta_weight[itheta];
            }
    
    double chi_sq = 0.0;
    for(int ir = 0; ir < n_r; ir++)
    {
        chi_sq += pow((rho_r[ir] - rho_r_std[ir]), 2);
    }
    return(chi_sq);
}

double MakeDensity::calculate_chisq_rho_r_theta(const gsl_vector *v, void* params)
{
    // compute the chi square of the nucleon density distribution along radial direction
    double ws_r0, ws_xsi, ws_beta2, ws_beta4;
    ws_r0 = gsl_vector_get(v, 0);
    ws_xsi = gsl_vector_get(v, 1);
    ws_beta2 = gsl_vector_get(v, 2);
    ws_beta4 = gsl_vector_get(v, 3);

    cout << ws_r0 << "   " << ws_xsi << "   " 
         << ws_beta2 << "   " << ws_beta4 << endl;
    calculate_density(ws_r0, ws_xsi, ws_beta2, ws_beta4);
    for(int ir = 0; ir < n_r; ir++)
    {
        for(int itheta = 0; itheta < n_theta; itheta++)
        {
            rho_r_theta[ir][itheta] = 0.0;
            rho_r_theta_std[ir][itheta] = 0.0;
        }
    }
    for(int ir = 0; ir < n_r; ir++)
        for(int itheta = 0; itheta < n_theta; itheta++)
            for(int iphi = 0; iphi < n_phi; iphi++)
            {
                rho_r_theta[ir][itheta] += rho[ir][iphi][itheta]*phi_weight[iphi];
                rho_r_theta_std[ir][itheta] += test_nucleus->standard_woods_saxon_distribution(r[ir], phi[iphi], cos_theta[itheta])*phi_weight[iphi];
            }
    
    double chi_sq = 0.0;
    for(int ir = 0; ir < n_r; ir++)
        for(int itheta = 0; itheta < n_theta; itheta++)
        {
            chi_sq += pow((rho_r_theta[ir][itheta] - rho_r_theta_std[ir][itheta]), 2);
        }

    return(chi_sq);
}

void MakeDensity::fit_nucleus_parameters()
{
    int deformed = paraRdr->getVal("deformed");
    if(deformed == 0)
        minimize_chisq();
    else
        minimize_chisq_deformed();
    return;
}

void MakeDensity::minimize_chisq()
{
    // minimization using gsl routine
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    /* Starting point */
    x = gsl_vector_alloc (2);
    gsl_vector_set (x, 0, test_nucleus->get_ws_r0_std());
    gsl_vector_set (x, 1, test_nucleus->get_ws_xsi_std());

    /* Set initial step sizes */
    ss = gsl_vector_alloc (2);
    gsl_vector_set (ss, 0, 1.0);
    gsl_vector_set (ss, 1, 0.3);

    /* Initialize method and iterate */
    double *paramsPtr = NULL;
    CCallbackHolder *Callback_params = new CCallbackHolder;
    Callback_params->clsPtr = this;
    Callback_params->params = paramsPtr;
    minex_func.n = 2;
    minex_func.f = this->CCallback_chisq_rho_r;
    minex_func.params = Callback_params;

    s = gsl_multimin_fminimizer_alloc (T, 2);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status) 
            break;
        
        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-2);

        if (status == GSL_SUCCESS)
        {
            printf ("converged to minimum at\n");
            ws_r0_best = gsl_vector_get(s->x, 0);
            ws_xsi_best = gsl_vector_get(s->x, 1);
        }
        printf ("%5lu %10.3e %10.3e chi_sq = %7.3f size = %.3f\n", iter, 
                gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), 
                s->fval, size);
    } while (status == GSL_CONTINUE && iter < 100);

    output_final_rho_r_vs_rho_r_std(ws_r0_best, ws_xsi_best);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
    return;
}

void MakeDensity::minimize_chisq_deformed()
{
    // minimization using gsl routine
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    /* Starting point */
    x = gsl_vector_alloc (4);
    gsl_vector_set (x, 0, test_nucleus->get_ws_r0_std());
    gsl_vector_set (x, 1, test_nucleus->get_ws_xsi_std());
    gsl_vector_set (x, 2, test_nucleus->get_ws_beta2_std());
    gsl_vector_set (x, 3, test_nucleus->get_ws_beta4_std());

    /* Set initial step sizes */
    ss = gsl_vector_alloc (4);
    gsl_vector_set (ss, 0, 1.0);
    gsl_vector_set (ss, 1, 0.3);
    gsl_vector_set (ss, 2, 0.3);
    gsl_vector_set (ss, 3, 0.3);

    /* Initialize method and iterate */
    double *paramsPtr = NULL;
    CCallbackHolder *Callback_params = new CCallbackHolder;
    Callback_params->clsPtr = this;
    Callback_params->params = paramsPtr;
    minex_func.n = 4;
    minex_func.f = this->CCallback_chisq_rho_r_theta;
    minex_func.params = Callback_params;

    s = gsl_multimin_fminimizer_alloc (T, 4);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status) 
            break;
        
        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-2);

        if (status == GSL_SUCCESS)
        {
            printf ("converged to minimum at\n");
            ws_r0_best = gsl_vector_get(s->x, 0);
            ws_xsi_best = gsl_vector_get(s->x, 1);
            ws_beta2_best = gsl_vector_get(s->x, 2);
            ws_beta4_best = gsl_vector_get(s->x, 3);
        }
        printf ("%5lu %10.3e %10.3e %10.3e %10.3e chi_sq = %7.3f size = %.3f\n", iter, 
                gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), 
                gsl_vector_get (s->x, 2), gsl_vector_get (s->x, 3), 
                s->fval, size);
    } while (status == GSL_CONTINUE && iter < 100);

    output_final_rho_r_theta_vs_rho_r_theta_std(ws_r0_best, ws_xsi_best, ws_beta2_best, ws_beta4_best);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
    return;
}

void MakeDensity::output_final_rho_r_vs_rho_r_std(double ws_r0_in, double ws_xsi_in)
{
    ofstream results("rho_r_vs_rho_r_std.dat");
    calculate_density(ws_r0_in, ws_xsi_in);
    for(int ir = 0; ir < n_r; ir++)
    {
        rho_r[ir] = 0.0;
        rho_r_std[ir] = 0.0;
    }
    for(int ir = 0; ir < n_r; ir++)
        for(int iphi = 0; iphi < n_phi; iphi++)
            for(int itheta = 0; itheta < n_theta; itheta++)
            {
                 rho_r[ir] += rho[ir][iphi][itheta]*phi_weight[iphi]*cos_theta_weight[itheta];
                 rho_r_std[ir] += test_nucleus->standard_woods_saxon_distribution(r[ir], phi[iphi], cos_theta[itheta])*phi_weight[iphi]*cos_theta_weight[itheta];
            }
    results << "# ws_r0 = " << ws_r0_in << " ws_xsi = " << ws_xsi_in << endl;
    results << "# r [fm]   rho_r [1/fm^3]   rho_r_std [1/fm^3]" << endl;
    for(int ir = 0; ir < n_r; ir++)
        results << scientific << setw(20) << setprecision(6) << r[ir] << "  " 
                << scientific << setw(20) << setprecision(6) << rho_r[ir] << "  " 
                << scientific << setw(20) << setprecision(6) << rho_r_std[ir]
                << endl;
    results.close();
    return;
}

void MakeDensity::output_final_rho_r_theta_vs_rho_r_theta_std(double ws_r0_in, double ws_xsi_in, double ws_beta2_in, double ws_beta4_in)
{
    ofstream results("rho_r_theta_vs_rho_r_theta_std.dat");
    calculate_density(ws_r0_in, ws_xsi_in, ws_beta2_in, ws_beta4_in);
    for(int ir = 0; ir < n_r; ir++)
    {
        for(int itheta = 0; itheta < n_theta; itheta++)
        {
            rho_r_theta[ir][itheta] = 0.0;
            rho_r_theta_std[ir][itheta] = 0.0;
        }
    }
    for(int ir = 0; ir < n_r; ir++)
        for(int iphi = 0; iphi < n_phi; iphi++)
            for(int itheta = 0; itheta < n_theta; itheta++)
            {
                 rho_r_theta[ir][itheta] += rho[ir][iphi][itheta]*phi_weight[iphi];
                 rho_r_theta_std[ir][itheta] += test_nucleus->standard_woods_saxon_distribution(r[ir], phi[iphi], cos_theta[itheta])*phi_weight[iphi];
            }
    results << "# r0 = " << ws_r0_in << " xsi = " << ws_xsi_in 
            << " beta2 = " << ws_beta2_in << " beta4 = " << ws_beta4_in << endl;
    results << "# r [fm]   rho_r_theta [1/fm^3]   rho_r_theta_std [1/fm^3]" << endl;
    int n_skip = 4;
    for(int ir = 0; ir < n_r; ir++)
    {
        results << scientific << setw(20) << setprecision(6) << r[ir] << "  ";
        for(int i = 0; i < (int)(n_theta/2./n_skip); i++)
        {
            results << scientific << setw(20) << setprecision(6) 
                    << rho_r_theta[ir][i*n_skip] << "  "
                    << scientific << setw(20) << setprecision(6) 
                    << rho_r_theta_std[ir][i*n_skip] << "  "; 
        }
        results << endl;
    }
    results.close();
    return;
}
