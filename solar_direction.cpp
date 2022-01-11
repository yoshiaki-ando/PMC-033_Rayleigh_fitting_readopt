#include <iostream>
#include <cmath>

#include <Vector3d.h>
#include "Date.h"
#include "pmc.h"
#include "Geoparameter.h"
#include "Mathparameter.h"

AndoLab::Vector3d <double> r_s(Date d){ /* Date d はUTで与えられる */
  constexpr double inclination_of_axis { 23.43 }; /* [deg] 公転の面からの地軸の傾き */
  constexpr double gamma { inclination_of_axis * Deg2rad };

  constexpr double minutes_of_a_day { 24. * 60 }; /* 1日の分数 */
  /*
   * ひまわりの経度で南中するUT
   * 0°E では UT12時で南中、15°E では UT11時で南中、
   * 同様に東経15°ごとに UT12時より 1時間ずつ早くなる
   */
  constexpr double minute_offset { (12. - Longitude_of_himawari/15.0) * 60.0 }; /* 2h37m12s */

  /*
   * 指定時刻における太陽の経度の角度φs
   *  1分あたり、1/(1日の分24*60) * 2π[rad]だけ
   * 反対(-φ方向)に周る。
   *
   * Formula A (Azimuthal angle of the sun)
   */
  const double phi_s
  = -2.0*M_PI * (d.minute_day() - minute_offset) / minutes_of_a_day;

  //  std::cout << phi_s / M_PI * 180.0 << std::endl;
  /*
   * 夏至からの日数、1年で規格化する
   */
  const double s = (d.doy_from_solstice() + d.minute_day()/minutes_of_a_day) / 365.;

  /*
   * Formula B (Polar angle of the sun)
   */
  const double cos_s_sin_g = std::cos(2.0*M_PI*s) * std::sin( gamma );
  const double th = std::atan2(
      std::sqrt( 1 - cos_s_sin_g * cos_s_sin_g ), cos_s_sin_g
      );

  constexpr double lat { 0.0 }; /* すべての緯度経度で同じベクトルとなる */
  constexpr double lon { Longitude_of_himawari };
  const double Ts { d.minute_day()/60.0 + 9.0 }; /* 時刻[h]（中央標準時）*/
  const double J { d.doy() + 0.5 }; /* day of year + 0.5 */
  constexpr double w { 2.0 * M_PI / 365.0 };
  const double delta = (0.33281 - 22.984*cos(w*J) - 0.34990*cos(2*w*J) - 0.13980*cos(3*w*J)
  + 3.7872*sin(w*J) + 0.03250*sin(2*w*J) + 0.07187*sin(3*w*J)) * Deg2rad; /* 太陽赤緯[rad] */

  const double e = 0.0072*cos(w*J) - 0.0528*cos(2*w*J) - 0.0012*cos(3*w*J)
  - 0.1229*sin(w*J) - 0.1565*sin(2*w*J) - 0.0041*sin(3*w*J); /* 均時差[h] */
  const double t = (15.0*(Ts + (lon - 135.0)/15.0 + e) - 180.0) *Deg2rad; /* 時角[rad] */
  const double h = asin(sin(lat*Deg2rad)*sin(delta) + cos(lat*Deg2rad)*cos(delta)*cos(t)); /* 太陽高度角[rad] */
  double th_h = M_PI/2.0 - h;//太陽天頂角[rad]

  const double sinA = cos(delta) * sin(t) / cos(h);
  const double cosA = (sin(h) * sin(lat*Deg2rad) - sin(delta)) / (cos(h)*cos(lat*Deg2rad));
  double A = atan2(sinA,cosA);//方位角、南=0、西=90
  // if(A > 2.0*M_PI){
  //   A = A - 2.0*M_PI;
  // }
  // if(A < 0.0){
  //   A = A + 2.0*M_PI;
  // }

  AndoLab::Vector3d <double> r;
  r = AndoLab::Vector3d <double> (1.0, (90.0 - lat)*Deg2rad, (lon - Longitude_of_himawari)*Deg2rad, AndoLab::coordinate::Spherical );

  /*太陽へ向かう単位ベクトル*/
  return (r.rotate(th_h, r.phi_vector())).rotate(A,-1.0*r);
  //return AndoLab::Vector3d <double> (1.0, th, phi_s, AndoLab::coordinate::Spherical );
}

/*
 * '21 9/3: E140.8°に修正。コメントを追加
 */
