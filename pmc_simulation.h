/*
 * pmc_rayleigh.h
 *
 *  Created on: 2019/07/09
 *      Author: ando
 */

#ifndef PMC_SIMULATION_H_
#define PMC_SIMULATION_H_

#include <iostream>
#include <fstream>

#include <Vector3d.h>
#include "Date.h"
#include "Mathparameter.h"
#include "Geoparameter.h"
#include "Geocoordinate.h"
#include "Msis.h"
#include "pmc.h"
#include "Region.h"

/* 観測データ (高度) */
constexpr double Obsrvd_data_Lowest_Altitude { 0.0 };
constexpr double Obsrvd_data_Highest_Altitude { 100.0e3 };
constexpr double Obsrvd_data_Step_Altitude { 1.0e3 };
constexpr int Number_of_Obsrvd_data_Altitude
{ int( (Obsrvd_data_Highest_Altitude - Obsrvd_data_Lowest_Altitude)
    / Obsrvd_data_Step_Altitude ) + 1 };

/* 観測データ (緯度経度) */
constexpr double Obsrvd_data_Lowest_Latitude { 38.0 };  /* [deg] */
constexpr double Obsrvd_data_Highest_Latitude { 81.0 }; /* [deg] */
constexpr double Obsrvd_data_Step_Latitude { 1.0 }; /* [deg] */
constexpr int Number_of_Obsrvd_data_Latitude
{ int( (Obsrvd_data_Highest_Latitude - Obsrvd_data_Lowest_Latitude)
    / Obsrvd_data_Step_Latitude + 1 + 0.5 ) * 2 };

/* 使用する波長 */
constexpr int Num_Lambda { 3 };
constexpr double Lambda[Num_Lambda] { 470e-9, 510e-9, 640e-9 }; /* 470nm, 510nm, 640nm */

/* 衛星視線上で大気圏上界面間を、分けて積分する区間 */
constexpr double Integral_Interval { 200.0e3 };
constexpr double PMC_Integral_Interval { 1.0e3 };

/* 観測ファイルの情報 */
constexpr double Lowest_latitude { 38.0 };
constexpr double Highest_latitude { 81.0 };
constexpr double Step_latitude { 1.0 };
constexpr int NUM_ALPHA
{ 2*int( (Highest_latitude - Lowest_latitude) / Step_latitude + 1 ) };

/*********************
 * 設定するパラメタ *
 *********************/

extern int Day_of_Year;
extern int Minute_of_Day;

/* 計算する高度 */
constexpr double Altitude_min { 10.e3 };
constexpr double Altitude_max { 100.e3 };
constexpr int N_alt { 90 };
constexpr double dAlt { (Altitude_max - Altitude_min) / N_alt };


/* PMC をフィッティングする領域 */
constexpr double Lower_Fitting { 60.0e3 };
constexpr double Upper_Fitting { 90.0e3 };
constexpr int idx_Lower_Fitting_PMC
{ int( (Lower_Fitting - Altitude_min)/dAlt + 0.5 ) };
constexpr int idx_Upper_Fitting_PMC
{ int( (Upper_Fitting - Altitude_min)/dAlt + 0.5 ) };

/* rayleighをフィッティングするための高度の上下限 */
constexpr double Alt_low_forFitting { 20.0e3 };
constexpr double Alt_high_forFitting { 43.0e3 };


/*******************************
 * 設定するパラメタ(ここまで) *
 *******************************/


/* パラメタ(派生) */

constexpr double R_lower_pmc { Radius_of_Earth + Lower_PMC };
constexpr double R_upper_pmc { Radius_of_Earth + Upper_PMC };

/* パラメタ(派生) ここまで */

/* インデックス */
constexpr int IDX_LATITUDE { 0 };
constexpr int IDX_LONGITUDE { 1 };

extern double* gauss_z;
extern double* gauss_w;

extern std::ofstream ofs_log;
extern int ga_number;

class Parameters{
private:
public:
  double wavelength;
  Msis msis;
};

/* 与えられた２点 r1, r2間の、与えられた時刻における光学的深さδを計算する
 * r1, r2 はひまわり座標
 */
double Optical_depth(
    const double Wavelength,
    Date,
    AndoLab::Vector3d <double> &r2,
    AndoLab::Vector3d <double> &r1,
    const int Number_of_PMC,
    PMC*
    );

/* 散乱断面積 */
double sigma(const double Theta, const double Wavelength, Msis );

/* 全散乱断面積 */
double sigma_t(const double Wavelength, Msis );

/* 体積散乱係数 */
double beta(AndoLab::Vector3d <double> r, void*);

/* 太陽を示すベクトル */
AndoLab::Vector3d <double> r_s(Date);

/* ひまわりから見て、北極から角度αで、高度Altitudeの球面とのリム点
 * αが正なら東経140°以下側、負なら東経140°以上 */
AndoLab::Vector3d <double> limb_point(double Altitude_in_m, double Angle_from_North_Pole_in_rad);

/* ひまわりから r の点を見た時、大気圏上空(高度100km)との交点2点を計算 */
AndoLab::Vector3d <double> *Across_point_atmosphere(AndoLab::Vector3d <double> r);

/* Dateの時刻における r Integral_Intervalの点でのレイリー散乱
 * 太陽から大気圏に入り、rの点に至るまでの減衰し、
 * rの点で散乱し、r2geo(大気圏界面上)に至るまで減衰した光強度
 */
double Rayleigh_scattering(
    AndoLab::Vector3d <double> r,
    AndoLab::Vector3d <double> r2geo,
    Date);

/* rの点から d 方向に進んで高度 H [m] となる点 */
AndoLab::Vector3d <double> cross_point_at_altitude(
    AndoLab::Vector3d <double> r,
    AndoLab::Vector3d <double> d,
    const double H);

AndoLab::Vector3d <double> *Across_point_upper_pmc(AndoLab::Vector3d <double>);
AndoLab::Vector3d <double> *Across_point_lower_pmc(AndoLab::Vector3d <double>);

bool is_illuminated(AndoLab::Vector3d <double> r, Date date);
AndoLab::Vector3d <double> upper_atmosphere_to_solar(
    AndoLab::Vector3d <double> r, Date date);

double intensity_integral(
    AndoLab::Vector3d <double> &r_from,
    AndoLab::Vector3d <double> &r_to,
    Date &date,
    const double Wavelength,
    const double th, /* 散乱各 */
    Msis &msis,
    double &delta,
    const int Number_of_PMC,
    PMC*,
    const double Interval_for_integration
    );


void output_pmc(AndoLab::Vector3d <double> Center,
    double SD_in_horizon,
    double SD_in_vertical);
AndoLab::Vector3d <double> rotate_alpha(
    AndoLab::Vector3d <double> r, const double alpha);
AndoLab::Vector3d <double> rotate_theta(
    AndoLab::Vector3d <double> r, const double theta);

/* r2->r1 の領域で、PMCの領域を通過する領域数と始点・終点ペアの探索 */
Region search_pmc_region(
    AndoLab::Vector3d <double> &r2,
    AndoLab::Vector3d <double> &r1
    );

void calculate_intensity(
    Date Day_of_Year,
    const int Number_of_PMC,
    PMC *pmc,
    const double Wavelength,
    const double Lower_Altitude,
    const double Step_Altitude,
    const int Number_of_Altitude,
    const double Lower_Altitude_to_calculate, /* 計算する高度下限 */
    const double Upper_Altitude_to_calculate, /* 計算する高度上限 */
    const double alpha,
    double *intensity             /* 返り値。各高度の光強度 */
    );

/* x軸上の衛星ひまわりから見たとき、点r の yz平面上への射影 */
AndoLab::Vector3d <double> projection_on_yz(AndoLab::Vector3d <double> r);

void get_observation_data(
    const int Number_of_Altitude,
    const double Lowest_Altitude,
    const double Step_of_Altitude,
    double ***Observed_data,
    double **Observed_Latitude_and_Longitude
    );
void get_rayleigh(
    Date,
    const int Number_of_angle,
    const double *alpha,
    const int Number_of_Altitude,
    const double Lowest_Altitude,
    const double Step_of_Altitude,
    double ***Calculated_Rayleigh_scattering
    );
void fitting_rayleigh(
    const int Number_of_angle,
    const double *alpha,
    const int Number_of_Altitude,
    const double Lowest_Altitude,
    const double Step_of_Altitude,
    double ***Calculated_Rayleigh_scattering,
    double ***Observed_data,
    int *Index_for_fitting_in_observed_data,
    double **Dark_Level,
    double **Coefficients
    );

/* 観測データのうち、フィッティングに使用する緯度経度インデックス配列を設定 */
void set_fitting_latlon(
    int *Index_array_of_Latitude_and_Longitude_for_fitting,
    double **Observed_Latitude_and_Longitude,
    const int Number_of_alpha,
    double **fitted_Latitude_and_Longitude
    );

void read_pmc_param( /* PMCのパラメタをファイルから読み込む */
    int &Number_of_PMC, /* PMC分布の数 */
    double** &pmc_param /* 各PMC分布のパラメタ */
    );

inline double lat2theta(double latitude){
  return (90 - latitude)*Deg2rad;
}
inline double phi(double theta){
  return acos( Radius_of_Earth / Rgeo / sin(theta) );
}

#endif /* PMC_SIMULATION_H_ */
