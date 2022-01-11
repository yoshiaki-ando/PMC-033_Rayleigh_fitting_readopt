/*
 * 陪ルジャンドル関数に関するルーチン
 *
 * 1. d P_n^1(cosθ) / dθ
 *    = cosθ { d P_n(cosθ)/d(cosθ) } + n(n+1) P_n(cosθ)
 *
 * 2. P_n^1(cosθ) / sinθ
 *    = - d { P_n(x) } / dx, θ = 0, π
 *
 */
#include <iostream>
#include <cmath>
#include <boost/math/special_functions/legendre.hpp>

#include "Mie_scattering.h"

double dP1n_dth(int n, double th){
  /*
   * 陪ルジャンドル関数の微分 d P_n^1(cosθ) / dθ
   *
   * d P_n^1(cosθ) / dθ = cosθ { d P_n(cosθ)/d(cosθ) }
   *                        + n(n+1) P_n(cosθ)
   *
   */
  double v;

  if ( th > M_PI ) th = M_PI;
  v = -std::cos(th) * P1n_over_sin(n, th)
      - n*(n+1) * boost::math::legendre_p(n, std::cos(th));

  return v;
}

double P1n_over_sin(int n, double theta){
  double v = 0.0;

  // エラー処理
  if ( (theta < 0.0) || (theta > M_PI) ){
    std::cout << "Error (P1n_over_sin) theta should be [0: pi]." << std::endl;
    std::exit(1);

  } else if ( (theta == 0) || (theta == M_PI) ){
    // lim      P_n^1(cosθ)     d P_n
    // θ→0,π ------------ = - ----- (1, -1)
    //            sinθ           dx
    v = - boost::math::legendre_p_prime(n, std::cos(theta));

  } else {
    v = boost::math::legendre_p(n, 1, std::cos(theta)) / std::sin(theta);
  }

  return v;
}
