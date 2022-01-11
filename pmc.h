/*
 * pmc_parameter.h
 *
 *  Created on: 2019/12/25
 *      Author: ando
 *
 * PMCの分布(空間、及び粒径)を決めるクラス
 * 空間はガウシアン、粒形は対数正規分布
 */

#ifndef PMC_H_
#define PMC_H_

#include <complex>

#include <Vector3d.h>

#include "Mie_scattering.h"


constexpr double radius_min { 10.0e-9 };
constexpr double radius_max { 260.0e-9 };

class PMC{
private:
  AndoLab::Vector3d <double> pR_center; /* [m,m,m] 分布の中心 */
  double pR_center_magnitude; /* [m] 分布の中心の地球中間からの距離 */

  double pSigma_z, pSigma_r;  /* [m] 分布の高度・水平方向の標準偏差 */

  double pR0;    /* [m] 粒径分布のモード半径 */
  double pSigma; /* 粒径分布の標準偏差みたいなもの */

  double pN0; /* [m^-3] 中心における数密度 */

  /* 散乱係数を求めるためのガウス積分の分点と重み */
  double gauss_w[8], gauss_z[8];

  double pNormalized_Beta;    /* 全散乱係数 */
  double pNormalized_Beta_th; /* 散乱係数 */

  Mie_scattering mie; /* Mie散乱の断面積, 3波長分 */

  void calc_normalized_beta(void);
  void set_gauss_legendre(void);

public:
  PMC(void);
  PMC(
      double Wavelength, /* [m], 散乱を求める際の波長 */
      std::complex <double> Refractive_Index, /* 散乱体の屈折率 */
      AndoLab::Vector3d <double> Center_position, /* PMCの空間分布の中心 */
      /* 空間分布の標準偏差, 高度方向と水平方向 */
      double SD_of_altitude_for_spatial_distribution,
      double SD_of_horizontal_direction_for_spatial_distribution,
      /* PMCの粒径分布
       * 粒径分布は「対数正規分布」と仮定する。
       * 1/(√(2π) logσ r) * exp( - log(r/r0)^2 / 2 / (logσ)^2 )
       * これのr0とσを与えると粒径分布 N_a(r) [m^{-1}] が決まる
       */
      double Mode_radius_r0,
      double SD_for_radius_distribution,
      /* 空間分布中心における全粒子数 */
      double Number_density
      );

//  void set_spatial_sd(double sigma_r, double sigma_z);
//  void center_position(FDTD::Vector <double> r_center){
//    pR_center = r_center;
//    pR_center_magnitude = pR_center.abs();
//  };
//  void set_radius_distribution(double r, double sigma);

  /* アクセサ */
  void set(
      double Wavelength, /* [m], 散乱を求める際の波長 */
      std::complex <double> Refractive_Index, /* 散乱体の屈折率 */
      AndoLab::Vector3d <double> Center_position, /* PMCの空間分布の中心 */
      /* 空間分布の標準偏差, 高度方向と水平方向 */
      double SD_of_altitude_for_spatial_distribution,
      double SD_of_horizontal_direction_for_spatial_distribution,
      /* PMCの粒径分布
       * 粒径分布は「対数正規分布」と仮定する。
       * 1/(√(2π) logσ r) * exp( - log(r/r0)^2 / 2 / (logσ)^2 )
       * これのr0とσを与えると粒径分布 N_a(r) [m^{-1}] が決まる
       */
      double Mode_radius_r0,
      double SD_for_radius_distribution,
      /* 空間分布中心における全粒子数 */
      double Number_density
      );

  AndoLab::Vector3d <double> rc(void){ return pR_center; }

  double sig_r(void){ return pSigma_r; }
  double sig_z(void){ return pSigma_z; }

  double r0(void){ return pR0; };
  double sigma(void){ return pSigma; };
  double N0(void){ return pN0; };

  double scs_t(double Radius){ return mie.scs_t(Radius); };
  void calc_normalized_beta_th(double theta);

  /* normalized total scattering coefficient β_t / N_0 */
  double bn(void){ return pNormalized_Beta; };

  /* normalized scattering coefficient β(θ) / N_0 */
  double bn_th(void){ return pNormalized_Beta_th; };
  double dist(AndoLab::Vector3d <double> r);

};

#endif /* PMC_H_ */
