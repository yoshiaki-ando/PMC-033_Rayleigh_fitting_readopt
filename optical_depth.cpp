#include <iostream>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <Vector3d.h>
#include <vector_gauss_quadrature.h>
#include "Across_pts.h"
#include "Date.h"
#include "Geocoordinate.h"
#include "Mie_scattering.h"
#include "pmc.h"
#include "Msis.h"
#include "pmc_simulation.h"


double mie_beta(AndoLab::Vector3d <double> r, void *parameters){

  PMC mp = *(PMC*)parameters;

  return mp.dist(r) * mp.bn();
}

double Mie_opt_depth(
    AndoLab::Vector3d <double> &r2, AndoLab::Vector3d <double> &r1,
    const int num_pmc, PMC *pmc){

  AndoLab::Vector3d <double> d = r2 - r1;
  const double distance = d.abs();
  const int Num = int(distance / PMC_Integral_Interval) + 1;
  AndoLab::Vector3d <double> dr = d/double(Num);

  double ans = 0.0;
  for(int i = 0; i < Num; i++){
    for(int npmc = 0; npmc < num_pmc; npmc++){
      ans += vector_gauss_quadrature(mie_beta, (void*)(pmc+npmc), r2+double(i)*dr, r2+double(i+1)*dr);
    }
  }

  return ans;
}

double Optical_depth(
    const double lambda, Date date,
    AndoLab::Vector3d <double> &r2, AndoLab::Vector3d <double> &r1,
    const int num_pmc, PMC *pmc){

  /* r2 -> r1 への光学的深さ */

  Parameters param;
  param.msis.set_date(date);  /* MSISE */
  param.wavelength = lambda;

  /* レイリー散乱による光学的深さ
   * (ガウス積分で全区間を積分してしまう→問題ない) */
//  const double distance { (r2-r1).abs() };
//  const int N_inte { int( distance / Integral_Interval ) + 1 };
//  AndoLab::Vector3d <double> dr = (r2-r1) / double(N_inte);
//  double delta { 0.0 };
//  for(int i = 0; i < N_inte; i++){
//    delta += gauss_integration_vector(beta, (void*)(&param),
//        r1 + double(i+1)*dr, r1 + double(i)*dr, gauss_z, gauss_w);
//  }
  double delta = vector_gauss_quadrature(beta, (void*)(&param), r2, r1);


  /* Mie散乱の寄与 */
  Region pmc_region = search_pmc_region(r2, r1);
  for(int i = 0; i < pmc_region.num; i++){
    delta += Mie_opt_depth(pmc_region.r_from[i], pmc_region.r_to[i],
        num_pmc, pmc);
  }
  return delta;

}
