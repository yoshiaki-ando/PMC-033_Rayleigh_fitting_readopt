#include <iostream>
#include "pmc_simulation.h"
#include "Across_pts.h"

void calculate_intensity(
    Date date, const int Num_PMC, PMC *pmc,
    const double lambda,
    const double Altitude_min, const double dAlt, const int N_alt,
    const double alt_low, const double alt_high,
    const double alpha,
    double *intensity){

  Msis msis(date);

  for(int iAlt = 0; iAlt < N_alt; iAlt++){

    std::cout << "Altitude : " << iAlt << " / " << N_alt << std::endl;
    intensity[iAlt] = 0.0;

    double alt = Altitude_min + iAlt*dAlt;

    if ( (alt < alt_low) || (alt > alt_high) ) continue;

    /* 高度alt, 北極からの角度αのリム点 */
    AndoLab::Vector3d <double> r_limb = limb_point(alt, alpha);
    /* 衛星からの視線と、大気圏上界との交点 */
    AndoLab::Vector3d <double> *Pts_atmos
    = Across_point_atmosphere(r_limb);

    /* PMC領域上面との交点 */
    AndoLab::Vector3d <double> Pts_upper_pmc[2], Pts_lower_pmc[2];
    if ( alt < Upper_PMC ){
      Across_pts acrs(Pts_atmos[0], Pts_atmos[1], Upper_PMC);
      acrs.copy(Pts_upper_pmc);
    }
    /* PMC領域下面との交点 */
    if ( alt < Lower_PMC ){
      Across_pts acrs(Pts_atmos[0], Pts_atmos[1], Lower_PMC);
      acrs.copy(Pts_lower_pmc);
    }

    /* 散乱角
     * Pts_atmos[0] - Pts_atmos[1] は、衛星の方向
     * -r_s は太陽光の入射方向
     */
    double th =
        angle_between( Pts_atmos[0] - Pts_atmos[1], -1.0*r_s(date));

    /* Mie散乱のための前処理(2) */
    for(int j_pmc = 0; j_pmc < Num_PMC; j_pmc++){
      pmc[j_pmc].calc_normalized_beta_th(th);
    }

    if ( alt < Lower_PMC ){

      double delta2r = 0.0;
      /* Pts_atmos[0]     -> Pts_upper_pmc[0] …散乱点はレイリー */
      intensity[iAlt] +=
          intensity_integral(Pts_atmos[0], Pts_upper_pmc[0],
          date, lambda,
          th, msis, delta2r, Num_PMC, pmc, Integral_Interval);

      /* Pts_upper_pmc[0] -> Pts_lower_pmc[0] */
      intensity[iAlt] +=
          intensity_integral(Pts_upper_pmc[0], Pts_lower_pmc[0],
          date, lambda,
          th, msis, delta2r, Num_PMC, pmc, PMC_Integral_Interval);

      /* Pts_lower_pmc[0] -> Pts_lower_pmc[1] */
      intensity[iAlt] +=
          intensity_integral(Pts_lower_pmc[0], Pts_lower_pmc[1],
          date, lambda,
          th, msis, delta2r, Num_PMC, pmc, Integral_Interval);

      /* Pts_lower_pmc[1] -> Pts_upper_pmc[1] */
      intensity[iAlt] +=
          intensity_integral(Pts_lower_pmc[1], Pts_upper_pmc[1],
          date, lambda,
          th, msis, delta2r, Num_PMC, pmc, PMC_Integral_Interval);

      /* Pts_upper_pmc[1] -> Pts_atmos[1]     */
      intensity[iAlt] +=
          intensity_integral(Pts_upper_pmc[1], Pts_atmos[1],
          date, lambda,
          th, msis, delta2r, Num_PMC, pmc, Integral_Interval);

    } else if ( alt < Upper_PMC ){ /* PMC内 */
      double delta2r = 0.0;
      /* Pts_atmos[0]     -> Pts_upper_pmc[0] …散乱点はレイリー */
      intensity[iAlt] +=
          intensity_integral(Pts_atmos[0], Pts_upper_pmc[0],
          date, lambda,
          th, msis, delta2r, Num_PMC, pmc, Integral_Interval);

      /* Pts_upper_pmc[0] -> Pts_lower_pmc[1] */
      intensity[iAlt] +=
          intensity_integral(Pts_upper_pmc[0], Pts_upper_pmc[1],
          date, lambda,
          th, msis, delta2r, Num_PMC, pmc, PMC_Integral_Interval);

      /* Pts_upper_pmc[1] -> Pts_atmos[1]     */
      intensity[iAlt] +=
          intensity_integral(Pts_upper_pmc[1], Pts_atmos[1],
          date, lambda,
          th, msis, delta2r, Num_PMC, pmc, Integral_Interval);

    } else {

      double delta2r = 0.0;
      intensity[iAlt] =
          intensity_integral(Pts_atmos[0], Pts_atmos[1],
          date, lambda,
          th, msis, delta2r, Num_PMC, pmc, Integral_Interval);

    }

  }

}
