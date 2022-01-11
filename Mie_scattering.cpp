#include <iostream>
#include <cmath>
#include <complex>

#include <complex_bessel.h>

#include <physical_constants.h>

#include "Mie_scattering.h"

/* 波長から周波数への変換 */
double lambda2f(double lambda){
  constexpr double c0 = C0;
  return ( c0 / lambda );
}

/* 波長から角周波数への変換 */
double lambda2omega(double lambda){
  return 2.0*M_PI*lambda2f(lambda);
}

/* d{ jn(z) }/dz = (n/z) * jn(z) - j_{n+1}(z) */
std::complex <double> djn(int n, std::complex <double> z){
  return (
      double(n)/z
      * sp_bessel::sph_besselJ(double(n), std::complex <double> (z))
  -   sp_bessel::sph_besselJ(double(n+1), std::complex <double> (z))
  );
}

/* d{ hn(z) }/dz = (n/z) * hn(z) - h_{n+1}(z) */
std::complex <double> dhn(int n, std::complex <double> z){
  return (
      double(n)/z
      * sp_bessel::sph_hankelH2(double(n), std::complex <double> (z))
  -   sp_bessel::sph_hankelH2(double(n+1), std::complex <double> (z))
  );
}

/* ψn = z jn(z) */
std::complex <double> psi_n(int n, std::complex <double> z){
  return ( z * sp_bessel::sph_besselJ( double(n), z ) );
}

/* ξn = z hn(z) */
std::complex <double> xi_n(int n, std::complex <double> z){
  return ( z * sp_bessel::sph_hankelH2( double(n), z ) );
}

/* (ψn)' = {z jn(z)}' = jn(z) + z jn(z)' */
std::complex <double> dpsi_n(int n, std::complex <double> z){
  return ( sp_bessel::sph_besselJ( double(n), z ) + z * djn(n,z) );
}

/* (ξn)' = {z hn(z)}' = hn(z) + z hn(z)' */
std::complex <double> dxi_n(int n, std::complex <double> z){
  return ( sp_bessel::sph_hankelH2( double(n), z ) + z * dhn(n,z) );
}

/* 係数 cn */
std::complex <double> Mie_scattering::cn(int n){
  std::complex <double> ka = std::complex <double> (p_ka);
  std::complex <double> k2a = p_m*ka;

  std::complex <double> value
  = - (psi_n(n, k2a) * dpsi_n(n, ka) - p_m*psi_n(n,ka) * dpsi_n(n,k2a))
  / (psi_n(n, k2a) * dxi_n(n, ka) - p_m*xi_n(n,ka) * dpsi_n(n,k2a));

  return value;
}

/* 係数 dn */
std::complex <double> Mie_scattering::dn(int n){
  std::complex <double> ka = std::complex <double> (p_ka);
  std::complex <double> k2a = p_m*ka;

  std::complex <double> value
  = - (p_m*psi_n(n, k2a) * dpsi_n(n, ka) - psi_n(n,ka) * dpsi_n(n,k2a))
  / (p_m*psi_n(n, k2a) * dxi_n(n, ka) - xi_n(n,ka) * dpsi_n(n,k2a));

  return value;
}

std::complex <double> Mie_scattering::S1n(int n, double th){
  return (
      (2*n+1)/double(n*(n+1))
      * ( cn(n)*dP1n_dth(n,th) + dn(n)*P1n_over_sin(n,th) )
      );
}

std::complex <double> Mie_scattering::S2n(int n, double th){
  return (
      (2*n+1)/double(n*(n+1))
      * ( cn(n)*P1n_over_sin(n,th) + dn(n)*dP1n_dth(n,th) )
      );
}

/* S1 */
std::complex <double> Mie_scattering::S1(double th){
  std::complex <double> sum( 0.0, 0.0 );
  for(int n = 1; n <= pN; n++){
    sum += S1n(n, th);
  }
  return sum;
}

/* S2 */
std::complex <double> Mie_scattering::S2(double th){
  std::complex <double> sum( 0.0, 0.0 );
  for(int n = 1; n <= pN; n++){
    sum += S2n(n, th);
  }
  return sum;
}

/* i1 */
double Mie_scattering::i1(double th){
  double a_s1 = std::abs( S1(th) );
  return a_s1*a_s1;
}

/* i2 */
double Mie_scattering::i2(double th){
  double a_s2 = std::abs( S2(th) );
  return a_s2*a_s2;
}

/* 散乱断面積(無偏光) σ(θ) = λ^2 / (4π^2) * (i1 + i2)/2 */
double Mie_scattering::scs(double th){
  return (
      pWavelength*pWavelength / 4.0 / M_PI / M_PI
      * ( i1(th) + i2(th) )/2.0
      );
}

double Mie_scattering::scs(double a, double th){
  p_a = a;
  p_ka = 2.0*M_PI/pWavelength * p_a;
  return (
      pWavelength*pWavelength / 4.0 / M_PI / M_PI
      * ( i1(th) + i2(th) )/2.0
      );
}

/* 全散乱断面積(無偏光、損失あり) */
double Mie_scattering::scs_t(void){
  double sum = 0.0;
  for(int n = 1; n <= pN; n++){
    sum += (2.0*n + 1.) * real( cn(n) + dn(n) );
  }
  return -pWavelength*pWavelength / 2.0 / M_PI * sum;
}

/* 全散乱断面積(無偏光、損失あり)、半径可変 */
double Mie_scattering::scs_t(double a){
  p_a = a;
  p_ka = 2.0*M_PI/pWavelength * p_a;

//  double sum = 0.0;
//  for(int n = 1; n <= pN; n++){
//    sum += (2.0*n + 1.) * real( cn(n) + dn(n) );
//  }
//  return -pWavelength*pWavelength / 2.0 / M_PI * sum;
  return scs_t();
}


