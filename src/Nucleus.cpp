#include <cstdlib>
#include <cmath>

#include "Nucleus.h"
#include "ParameterReader.h"
#include "gauss_quadrature.h"

using namespace std;

Nucleus::Nucleus(ParameterReader* paraRdr_in)
{
    paraRdr = paraRdr_in;
    atomic_num = paraRdr->getVal("atomic_num");
    atomic_mass = (double) atomic_num;
    deformed = paraRdr->getVal("deformed");
    r_min = paraRdr->getVal("r_min");

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
        rho_0 = 0.16;
    }
    else if (atomic_num == 238)
    {
        //Taken from P.Filip et al.,PRC80,054903(2009)
        r_0_std = 6.81;
        xsi_std = 0.54;
        rho_0 = 0.166;
    }

    if(deformed)
    {
        // Parameters taken from Moller et al., Atomic Data and Nuclear Data
        // Tables, 59, 185-381, (1995).
        if(atomic_num == 197)  //(Au)
        {
            beta2_std = -0.13;
            beta4_std = -0.03;
        }
        else if(atomic_num == 63)  //(Cu)
        {
            beta2_std = 0.162;
            beta4_std = 0.006;
        }
        else if(atomic_num == 238)  //(U)
        {
            beta2_std = 0.28;
            beta4_std = 0.093;
        }
        else
        {
            cerr << "Mass number not available for reparametrization" << endl;
        }
    }
    else
    {
        beta2_std = 0.0; beta4_std = 0.0;
    }


    r_0 = r_0_std;
    xsi = xsi_std;
    beta2 = beta2_std;
    beta4 = beta4_std;

    //determine rho_0 from atomic number
    double r_max = paraRdr->getVal("r_max");
    int n_r = paraRdr->getVal("n_r");
    int n_phi = paraRdr->getVal("n_phi");
    int n_theta = paraRdr->getVal("n_theta");

    double *r = new double [n_r];
    double *r_weight = new double [n_r];
    gauss_quadrature(n_r, 1, 0.0, 0.0, 0.0, r_max, r, r_weight);
    double *phi = new double [n_phi];
    double *phi_weight = new double [n_phi];
    gauss_quadrature(n_phi, 1, 0.0, 0.0, 0.0, 2*M_PI, phi, phi_weight);
    double *cos_theta = new double [n_theta];
    double *cos_theta_weight = new double [n_theta];
    gauss_quadrature(n_theta, 1, 0.0, 0.0, -1.0, 1.0, cos_theta, cos_theta_weight);

    double temp = 0.0;
    for(int ir = 0; ir < n_r; ir++)
        for(int iphi = 0; iphi < n_phi; iphi++)
            for(int itheta = 0; itheta < n_theta; itheta++)
                 temp += standard_woods_saxon_distribution(r[ir], phi[iphi], cos_theta[itheta])*r[ir]*r[ir]*r_weight[ir]*phi_weight[iphi]*cos_theta_weight[itheta];

    rho_0 = rho_0*atomic_num/temp;

    delete [] r;
    delete [] r_weight;
    delete [] phi;
    delete [] phi_weight;
    delete [] cos_theta;
    delete [] cos_theta_weight;

    ws_max = rho_0;

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

void Nucleus::set_woods_saxon_parameters(double r0_in, double xsi_in, double beta2_in, double beta4_in)
{
    r_0 = r0_in;
    xsi = xsi_in;
    beta2 = beta2_in;
    beta4 = beta4_in;
}

void Nucleus::sample_nucleon_corrdinate(double& x, double& y, double& z)
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
    nucleus.erase(nucleus.begin(), nucleus.end());
    double xcm = 0.0;
    double ycm = 0.0;
    double zcm = 0.0;
    for(int i_nucleon=0; i_nucleon < atomic_num; i_nucleon++)
    {
        double x,y,z;
        int icon=0;
        do
        {
            sample_nucleon_corrdinate(x,y,z);
            icon=0;
            for(int i = 0; i<(int)nucleus.size(); i++)
            {
                double x1=nucleus[i]->getX();
                double y1=nucleus[i]->getY();
                double z1=nucleus[i]->getZ();
                double r2 = (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1);
                if(r2 < r_min)
                {
                    icon=1;
                    break;
                }
            }
        } while(icon==1);

        xcm +=x;
        ycm +=y;
        zcm +=z;
        nucleus.push_back(new Particle(x,y,z));
    }

    for(int i_nucleon = 0; i_nucleon < atomic_num; i_nucleon++)
    { 
    // shift center of nucleus
        double x = nucleus[i_nucleon]->getX() - xcm/atomic_num;
        double y = nucleus[i_nucleon]->getY() - ycm/atomic_num;
        double z = nucleus[i_nucleon]->getZ() - zcm/atomic_num;
        nucleus[i_nucleon]->setX(x);
        nucleus[i_nucleon]->setY(y);
        nucleus[i_nucleon]->setZ(z);
    }

}

void Nucleus::get_nucleon_position(int index, double& x, double& y, double& z)
{
    x = nucleus[index]->getX();
    y = nucleus[index]->getY();
    z = nucleus[index]->getZ();
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
