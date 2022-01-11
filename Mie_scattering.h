/*
 * mie_scattering.h
 *
 *  Created on: 2019/07/01
 *      Author: ando
 */

#ifndef MIE_SCATTERING_H_
#define MIE_SCATTERING_H_

#include <complex>
#include <string>

class Mie_scattering{
private:
  double pWavelength; /* [m] 波長 */
  double p_a; /* [m] 散乱球の半径 */
  std::complex <double> p_m; /* [-] 屈折率 */
  double p_ka; /* [rad] 真空中の波数k0 × 半径 a */

  int pN; /* 級数の項数 */
public:
  Mie_scattering(void){};
  Mie_scattering(double lambda, double a, std::complex <double> m):
    pWavelength( lambda ), p_a( a ), pN( 20 ){
    p_m = m;
    p_ka = 2.0*M_PI/lambda * p_a;
  }

  void set(double lambda, double a, std::complex <double> m){
    pWavelength = lambda;
    p_a = a;
    pN = 20;
    p_m = m;
    p_ka = 2.0*M_PI/lambda * p_a;
  }

  void set_series_number(int N){ pN = N; }
  std::complex <double> cn(int n);
  std::complex <double> dn(int n);
  std::complex <double> S1n(int n, double th);
  std::complex <double> S2n(int n, double th);
  std::complex <double> S1(double th);
  std::complex <double> S2(double th);
  double i1(double th); /* i1 */
  double i2(double th); /* i2 */

  /* 散乱断面積(無偏光) σ(θ) = λ^2 / (4π^2) * (i1 + i2)/2 */
  double scs(double th);
  double scs(double a, double th);

  /* 全散乱断面積(無偏光、損失あり) */
  double scs_t(void);
  double scs_t(
      double a /* [m] 散乱級の半径 */
      );
};

double lambda2f(double lambda); /* lambda [m] ==> f [Hz]の変換 */
double lambda2omega(double lambda); /* lambda [m] ==> ω [rad/sec]の変換 */

/* 陪ルジャンドル関数 P_n^1(cosθ) / sinθ */
double P1n_over_sin(int n, double theta);

/* 陪ルジャンドル関数の微分 d P_n^1(cosθ) / dθ */
double dP1n_dth(int n, double th);

void set_options( /* オプションによるパラメタ設定 */
    int argc, char **argv,
    double &Radius_of_scatterer,
    double &Wavelength,
    std::complex <double> &Refractive_index,
    int &Number_of_points_for_scattering_angle,
    double &Theta_min, double &Theta_max,
    std::string &Output_filename
    );

std::complex <double> interpolated_m(double wavelength_in_m);

/* パラメタの出力 */
void output_param(
    const double Radius_of_sphere_in_m,
    const double Wavelength_in_m,
    std::complex <double> Reflactive_index,
    const int Number_of_sampling_for_scattering_angle,
    const double Minimum_value_of_Theta,
    const double Maximum_value_of_Theta);

#endif /* MIE_SCATTERING_H_ */
