#include <cstdlib>
#include <cmath>

#include "Nucleus.h"
#include "ParameterReader.h"

using namespace std;

Nucleus::Nucleus(ParameterReader* paraRdr_in)
{
    paraRdr = paraRdr_in;
    atomic_num = paraRdr->getVal("atomic_num");
    atomic_mass = (double) atomic_num;
    deformed = paraRdr->getVal("deformed");

    // generic (crude) default parameterization of radius and surface thickness
    rho_0 = 0.17;
    r_0_std = 1.12*pow(atomic_mass, 1./3.)-0.86/pow(atomic_mass, 1./3.);
    xsi_std = 0.54;
    //taken from C.W.De Jager et al. Atom.Data.Nucl.Data Tabl.36, 495 (1987).
    if(atomic_num == 197)
    {
      r_0_std = 6.38;
      xsi_std = 0.535;
      rho_0 = 0.1695;
    }
    else if (atomic_num == 63)
    {
      r_0_std = 4.20641;
      xsi_std = 0.5877;
      rho_0 = 0.1686;
    }
    else if (atomic_num == 208)
    {
        r_0_std = 6.62;
        xsi_std = 0.546;
        rho_0 = 0.17;
    }
    else if (atomic_num == 238)
    {
        //Taken from P.Filip et al.,PRC80,054903(2009)
        r_0_std = 6.81;
        xsi_std = 0.54;
        rho_0 = 0.17;
    }

    r_0 = r_0_std;
    xsi = xsi_std;
    beta2 = beta2_std;
    beta4 = beta4_std;

    ws_max = rho_0/(1. + exp(-r_0/xsi));

    r_max = r_0 + 5.0*xsi;
}

Nucleus::~Nucleus()
{

}

void Nucleus::set_woods_saxon_parameters(double r0_in, double xsi_in)
{
    r_0 = r0_in;
    xsi = xsi_in;
}

void Nucleus::get_nucleon_corrdinate(double& x, double& y, double& z)
{
    double r;
    double cos_theta;
    double phi;

    // uniform sampling in a sphere within r = r_max
    do
    {
        r = r_max*pow(drand48(), 1.0/3.0);
        cos_theta = 1.0-2.0*drand48();
        phi = 2*M_PI*drand48();
    } while(drand48()*ws_max > woods_saxon_distribution(r, phi, cos_theta));
    
    double sin_theta = sqrt(1. - cos_theta*cos_theta);
    x = r*sin_theta*cos(phi);
    y = r*sin_theta*sin(phi);
    z = r*cos_theta;

    return;
}

void Nucleus::generate_nucleus()
{

}

double Nucleus::woods_saxon_distribution(double r, double phi, double cos_theta)
{
    double density;

    if(deformed)
    {
        double y20 = spherical_harmonics(2, cos_theta);
        double y40 = spherical_harmonics(4, cos_theta);
        double r_0_deformed =r_0*(1.0 + beta2*y20 + beta4*y40);
        density = rho_0/(1.0 + exp((r - r_0_deformed)/xsi));
    }
    else
        density = rho_0/(1.0 + exp((r - r_0)/xsi));
    return(density);
}

double Nucleus::standard_woods_saxon_distribution(double r, double phi, double cos_theta)
{
    double density;

    if(deformed)
    {
        double y20 = spherical_harmonics(2, cos_theta);
        double y40 = spherical_harmonics(4, cos_theta);
        double r_0_deformed =r_0_std*(1.0 + beta2_std*y20 + beta4_std*y40);
        density = rho_0/(1.0 + exp((r - r_0_deformed)/xsi_std));
    }
    else
        density = rho_0/(1.0 + exp((r - r_0_std)/xsi_std));
    return(density);
}

double Nucleus::spherical_harmonics(int l, double cos_theta)
{
  //Currently assuming m=0 and available for Y_{20} and Y_{40}

  double ylm = 0.0;

  if(l == 2){

    ylm = 3.0*cos_theta*cos_theta-1.0;
    ylm *= 0.31539156525252005; //pow(5.0/16.0/M_PI,0.5);

  }else if (l == 4){

    ylm = 35.0*pow(cos_theta, 4);
    ylm -= 30.0*cos_theta*cos_theta;
    ylm += 3.0;
    ylm *= 0.10578554691520431; //3.0/16.0/pow(M_PI,0.5);

  }else{
    cerr << "Not available in Nucleus::spherical_harmonics" << endl;
  }
  
  return ylm;
}
