#include <iomanip>

#include "Regge96.h"
#include "ParameterReader.h"
#include "Nucleus.h"
#include "MakeDensity.h"
#include "gauss_quadrature.h"

using namespace std;

MakeDensity::MakeDensity(ParameterReader *paraRdr_in)
{
    paraRdr = paraRdr_in;
    n_event = 1000;

    // NN cross sections in mb
    double ecm = paraRdr->getVal("ecm");
    double sig = hadronxsec::totalXsection(ecm, 0);
    double sigel = hadronxsec::elasticXsection(sig, ecm, 0, 0);
    siginNN = sig - sigel;
    gauss_nucl_width = sqrt(siginNN/(8.*M_PI));

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
    for(int i = 0; i < n_r; i++)
    {
        rho[i] = new double* [n_phi];
        rho_r[i] = 0.0;
        rho_r_std[i] = 0.0;
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
    }
    delete[] rho;
    delete[] rho_r;
    delete[] rho_r_std;

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
    double prefactor = 1./(pow(gauss_nucl_width, 3)*pow(M_PI, 1.5));
    double sigma_sq = gauss_nucl_width*gauss_nucl_width;
    test_nucleus->set_woods_saxon_parameters(ws_r0, ws_xsi);
    for(int i = 0; i < n_event; i ++)
    {
        cout << "processing event " << i << endl;
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

double MakeDensity::calculate_rho_r(double ws_r0, double ws_xsi)
{
    calculate_density(ws_r0, ws_xsi);
    for(int ir = 0; ir < n_r; ir++)
        for(int iphi = 0; iphi < n_phi; iphi++)
            for(int itheta = 0; itheta < n_theta; itheta++)
            {
                 rho_r[ir] += rho[ir][iphi][itheta]*phi_weight[iphi]*cos_theta_weight[itheta];
                 rho_r_std[ir] += test_nucleus->standard_woods_saxon_distribution(r[ir], phi[iphi], cos_theta[itheta])*phi_weight[iphi]*cos_theta_weight[itheta];
            }
    return(0.0);
}

double MakeDensity::compare_rho_r_with_standard_ws()
{
    double chi_sq = 0.0;
    for(int ir = 0; ir < n_r; ir++)
    {
        chi_sq += pow((rho_r[ir] - rho_r_std[ir]), 2);
    }
    return(chi_sq);
}
