/*
 * PMC分布(2次元を出力するプログラム)
 */
#include <iostream>
#include <fstream>
#include <cmath>

#include <Vector3d.h>

#include "pmc.h"
#include "Mathparameter.h"
#include "Geoparameter.h"
#include "pmc_simulation.h"

void output_pmc(PMC pmc){

  Geocoordinate g_c( pmc.rc() );
  std::cout << "N" << g_c.latitude() << ", E"
      << g_c.longitude() << std::endl;

  /* 地軸からの角度。本来、渡されるものであるが、
   * これは一時的なプログラムであるため、ここで定義 */
  constexpr double ALPHA { -4.0 * Deg2rad };

  constexpr int NUM_H { 300 };
  constexpr int NUM_Z { 50 };

  const double z0 { pmc.rc().r() }; /* 高度(+地球の半径)の中心 */
  const double z_min { z0 - 4.0*pmc.sig_z() };
  const double z_max { z0 + 4.0*pmc.sig_z() };
  const double dz { (z_max - z_min)/NUM_Z };

  const double cos_b { z0 / Rgeo };
  const double th0 { acos( cos_b ) };

  const double sig_th { pmc.sig_r() / z0 }; /* 地球の中心からの角度のσ */

  const double th_min { th0 - 4.*sig_th };
  const double th_max { th0 + 4.*sig_th };
  const double dth { (th_max - th_min) / NUM_H };

  std::ofstream ofs("pmc_dist_p900.dat");
  for(int i_z = 0; i_z <= NUM_Z; i_z++){
    double alt { z_min + i_z * dz };

    for(int j_th = 0; j_th <= NUM_H; j_th++){
      double th { th_min + j_th * dth };

      AndoLab::Vector3d <double> r0( alt*cos(th), 0.0, alt*sin(th) );
      AndoLab::Vector3d <double> r = rotate_alpha(r0, ALPHA);

      ofs << r0.x() << " " << r0.z() << " "
          << pmc.dist(r) << "\n";
    }
    ofs << "\n";

  }

}
