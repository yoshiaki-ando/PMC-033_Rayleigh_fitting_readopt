#include <Vector3d.h>
#include "Geocoordinate.h"
#include "Msis.h"
#include "pmc_simulation.h"

/* 波長λ[m]のときの標準大気の屈折率 ー １
 *
 * from [1]
 *
 */
double ns_1(double lambda){
  lambda *= 1.0e6; /* convert from [m] to [μm]  */
  const double lambda_m2 = 1.0/lambda/lambda; /* λ^{-2} */

  /* 文献[1] 式(1)の右辺 */
  double index
  = 6432.8 + 2949810.0 / (146.0 - lambda_m2)
  + 25540. / (41.0 - lambda_m2);

  return index / 1.0e8;
}

/* King因子 from Tomasi */
double f(double lambda){
  lambda *= 1e6; /* λ[μm] */
  double lambda_2 = 1.0/lambda/lambda;
  double lambda_4 = lambda_2 * lambda_2;
  return ( (0.78084 * (1.034 + 3.17e-4*lambda_2)
      + 0.20946 * (1.096 + 1.385e-3*lambda_2 + 1.448e-4*lambda_4)
      + 0.00934) / 0.99964 );
}

/* depolarization factor */
double rho(const double lambda){
  return 6.*( f(lambda) - 1.0 ) / ( 7.0*f(lambda) + 3.0 );
}

double gam(const double lambda){
  return rho(lambda) / ( 2.0 - rho(lambda) );
}

double Phase_func(const double Theta, const double lambda){
  return 3./4. * ( 1. + 3.*gam(lambda) + (1. - gam(lambda))
      * std::cos(Theta)*std::cos(Theta) )
      / ( 1.0 + 2.0 * gam(lambda) );

}

double sigma(const double th, const double lambda, Msis msis){
  return sigma_t(lambda,msis)/4.0/M_PI * Phase_func(th, lambda);
}

double sigma_t(const double lambda, Msis msis){

  constexpr double Pressure_of_standard_air = 101325.; /* [Pa] */
  constexpr double Temperature_of_standard_air = (273.15 + 15.); /* [K] */
  const double lambda2 = lambda*lambda;
  const double lambda4 = lambda2*lambda2;

  const double n_1 = (msis.p() / msis.T())
      / (Pressure_of_standard_air / Temperature_of_standard_air)
      * ns_1(lambda);
  const double n = n_1 + 1.0;

  double sigma = 24.0 * M_PI * M_PI * M_PI
      / lambda4 / msis.N() / msis.N()
      * n_1 * n_1 * (n + 1.0) * (n + 1.0)
      / (n*n + 2.0) / (n*n + 2.0)
      * f(lambda);

  return sigma;
}

double beta(AndoLab::Vector3d <double> r, void* p){
  Parameters *param = (Parameters*)p;
  param->msis.set_position( Geocoordinate ( r ) );

  return sigma_t( param->wavelength, param->msis ) * param->msis.N();
}
