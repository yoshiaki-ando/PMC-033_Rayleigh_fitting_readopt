#include <cmath>
#include <Vector3d.h>
#include "pmc_simulation.h"

/* rの点から、d方向へ進んで高度 H (原点からの距離 = 地球の半径 + H)
 * との交点を求める
 */
AndoLab::Vector3d <double> cross_point_at_altitude(
    AndoLab::Vector3d <double> r,
    AndoLab::Vector3d <double> d, /* 単位ベクトルとすること */
    const double H /* [m] 高度 */
    ){

  const double R_H = Radius_of_Earth + H; /* 原点からの距離 */
  double k = -1.0*r%d + std::sqrt(
      (r%d)*(r%d) - r.abs()*r.abs() + R_H*R_H );
  return r + k*d;
}

/* x軸上の衛星ひまわりから見たとき、点r の yz平面上への射影 */
AndoLab::Vector3d <double> projection_on_yz(AndoLab::Vector3d <double> r){

  AndoLab::Vector3d <double> r_geo { Rgeo, 0.0, 0.0 }; /* ひまわりの位置 */

  /* r0 = R r^(θ, φ)
   * d^ = (r0 - Rgeo)/|r0 - Rgeo|
   * (Rgeo + t d^)・x^ = 0 となるtを探す
   * Rgeo.x + t d^.x = 0 ←→ t = - Rgeo.x / d^.x;
   * (Rgeo + t d^)・z^ / |Rgeo + t d^| = cosα
   **/

  AndoLab::Vector3d <double> dn = (r - r_geo).n();
  double t = - r_geo.x() / dn.x();
  return r_geo + t * dn;
}

