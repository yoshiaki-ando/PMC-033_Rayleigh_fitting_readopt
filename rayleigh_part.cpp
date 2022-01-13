#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include <Vector3d.h>
#include "Date.h"
#include "Msis.h"
#include "pmc_simulation.h"
#include "albedo.h"
#include "optdepth.h"

bool is_illuminated(AndoLab::Vector3d <double> r, Date date){
  /* 大気中の点 r に太陽光が入射しているか */

  double rQ = AndoLab::abs( r_s(date) * (r * r_s(date)) );

  return ( r%r_s(date) >= 0.0 ) ||
      ( (r%r_s(date) < 0.0) && (rQ >= Radius_of_Earth) );
}


bool is_in_shadow_zone(AndoLab::Vector3d <double> r, Date date){
  /* 大気中の点 r が影領域に入っているか */
  return !is_illuminated(r,date);
}


AndoLab::Vector3d <double> upper_atmosphere_to_solar(
    AndoLab::Vector3d <double> r, Date date){
  /* r の点から太陽方向へ移動して大気圏界面に達する点 */

  return cross_point_at_altitude(r, r_s(date).n(),
      Altitude_of_Atmosphere );
}


bool is_inside_PMC_region(AndoLab::Vector3d <double> r){
  return ( (r.abs() > R_lower_pmc) && (r.abs() < R_upper_pmc) );
}


double intensity_integral(
    AndoLab::Vector3d <double> &r_from,
    AndoLab::Vector3d <double> &r_to,
    Date &date,
    const double lambda,
    const double th, /* 散乱角 */
    Msis &msis,
    double &delta,
    const int num_pmc,
    PMC *pmc,
    const double Interval){

  double I { 0.0 };

  /* 寄与を計算するベクトル */
  AndoLab::Vector3d <double> length = r_to - r_from;

  /* 分割点数。Interval毎 ===>> 要調査 */
  const int Num = int( length.abs() / Interval + 1 );
  AndoLab::Vector3d <double> dr = length / double(Num); /* 分割距離 */

  AndoLab::Vector3d <double> pre_rp = r_from;

  /* 反射のための各種設定 */
  AndoLab::Vector3d <double> yv( 0.0, 1.0, 0.0 );
  AndoLab::Vector3d <double> zv( 0.0, 0.0, 1.0 );
  AndoLab::Vector3d <double> r_gs = (r_from - r_to).n();

  for(int i_dr = 0; i_dr < Num; i_dr++){
    AndoLab::Vector3d <double> r_p = r_from + (i_dr+0.5)*dr;

    /* 視線方向の光学的深さは、視線から遠い方向へと計算してゆけば
     * 先の計算が使える */
    delta += Optical_depth(lambda, date, pre_rp, r_p, num_pmc, pmc);

    if ( is_illuminated(r_p, date) ){
      msis.set_position( Geocoordinate(r_p) );

      /* r の点がシャドウ領域でなければ
       * 1. 太陽方向の大気圏界面の点 r0 を見つける
       * 2. exp( - delta(r0,r) - delta(r,r2geo) ) * σ(θ) */
      AndoLab::Vector3d <double> r_i =
          upper_atmosphere_to_solar(r_p, date);

      /* 散乱係数 */
      double total_beta = sigma(th, lambda, msis) * msis.N();
      if ( is_inside_PMC_region( r_p ) ){
        for(int npmc = 0; npmc < num_pmc; npmc++){
          total_beta += pmc[npmc].bn_th() * pmc[npmc].dist( r_p );
        }
      }

      /* 反射 */
      double reflection = 0.0; /* 単一分子当たりの反射光 */
      const double alpha_max = std::acos( R0 / r_p.r() );
      const double Dalpha = alpha_max / M_alpha;
      const double theta_s = r_p.theta();
      const double phi_s = r_p.phi();

      /* 反射に関する光学的深さを読み出す */
      double *optdepth = new double [M_alpha * N_beta * 2];
      ifs_optdepth.read( (char*)optdepth, sizeof(double) * M_alpha*N_beta*2 );

      /** 光のあたる反射面数を求める */

      for(int m = 0; m < M_alpha; m++){
        const double alpha = (m + 0.5) * Dalpha;

        for(int n = 0; n < N_beta; n++){
          AndoLab::Vector3d <double> rp_R_mn(R0, alpha, (n+0.5)*Dbeta, AndoLab::coordinate::Spherical);
          AndoLab::Vector3d <double> r_R_mn = ( rp_R_mn.rotate(theta_s, yv) ).rotate(phi_s, zv); /* 反射点 */

          if ( r_R_mn%r_s(date) <= 0.0 ){
            continue;
          }

          const double Smn = R0*R0 * std::sin(alpha) * Dalpha * Dbeta;

          AndoLab::Vector3d <double> r_Rs_mn = r_p - r_R_mn; /* 反射点→散乱点ベクトル */
          AndoLab::Vector3d <double> r_Ri_mn =
              upper_atmosphere_to_solar(r_R_mn, date); /* 反射点に入射する光の大気圏入射点 */
          /* 反射点到達光 */
          double optdepth1 = optdepth[m*N_beta*2 + n*2 + 0];
          double I_r_mn = std::exp( - optdepth1 );

          double Ip_S_mn = Reflection_Coefficient * Smn / M_PI / r_Rs_mn.abs() / r_Rs_mn.abs()
            * I_r_mn * r_R_mn.n()%r_s(date) * r_Rs_mn.n()%r_R_mn.n(); /* 無損失散乱光 */
            //                      * I_r_mn * r_R_mn.n()%r_p.n() * r_Rs_mn.n()%r_R_mn.n(); /* 無損失散乱光 */
          /* 損失散乱光 */
          double optdepth2 = optdepth[m*N_beta*2 + n*2 + 1];
          double I_S_mn = Ip_S_mn * std::exp( - optdepth2 );
          double theta_s_mn = AndoLab::angle_between( r_Rs_mn, r_gs  ); /* 反射光散乱角 */
//        ofs_logmn << i_dr << " " << m << " " << n << " " << I_S_mn * sigma(theta_s_mn, lambda, msis) << std::endl;
          reflection += I_S_mn * sigma(theta_s_mn, lambda, msis);
        }
      }

      /* 反射の検討(ここまで) */

      I += std::exp(-delta) *
          (std::exp( - Optical_depth(lambda, date, r_i, r_p, num_pmc, pmc) ) * total_beta
              + reflection * msis.N() )
              * dr.abs();
    }
    pre_rp = r_p;
  }

  delta += Optical_depth(lambda, date, pre_rp, r_to, num_pmc, pmc );

  return I;
}
