#include <iostream>
#include <fstream>

#include <Vector3d.h>

#include "Mathparameter.h"
#include "Geoparameter.h"
#include "pmc_simulation.h"
#include "optdepth.h"

std::ifstream ifs_optdepth;

void get_rayleigh(
    Date date,
    const int Num_alpha,
    const double *alpha,
    const int N_alt,
    const double Altitude_min,
    const double dAlt,
    double ***intensity
    ){

  /* 数密度0のPMCを設定 */

  /* PMCの設定 */
  AndoLab::Vector3d <double> r_pmc_center
  = rotate_alpha( rotate_theta( limb_point(83.e3, 0.0*Deg2rad),
      0.e3/(83.e3+Radius_of_Earth) ),
      -4.0*Deg2rad );

  const std::complex <double>
  Reflactive_index_of_Ice { 1.3126, 8.036e-10 }; /* 要出典 */

  PMC pmc1;

  ifs_optdepth.open( ("data/" + str_optdepth).c_str(), std::ios::in|std::ios::binary);

  for(int j_lambda = 0; j_lambda < 3; j_lambda++){
    pmc1.set(Lambda[j_lambda], Reflactive_index_of_Ice,
        r_pmc_center, 0.6e3, 300.e3, 80.0e-9, 1.4, 0.0 );

    for(int i_alpha = 0; i_alpha < Num_alpha; i_alpha++){
      std::cout << j_lambda << " : "
          << i_alpha << " / " << Num_alpha << std::endl;

      calculate_intensity(date, 1, &pmc1, Lambda[j_lambda],
          Altitude_min, dAlt, N_alt,
          Altitude_min, Altitude_max,
          alpha[i_alpha], intensity[j_lambda][i_alpha]);

    }
  }

  ifs_optdepth.close();

  for(int i = 0; i < 3; i++){
    std::ofstream ofs1( ("data/rayleigh_" + std::to_string(i) + ".dat").c_str() );

    for(int j = 0; j < Num_alpha; j++){
      for(int k = 0; k < N_alt; k++){
        ofs1 << alpha[j] << " " << Altitude_min + k*dAlt << " "
            << intensity[i][j][k] << "\n";
        //          std::cout << alpha[j] << " " << Altitude_min + k*dAlt << " "
        //              << intensity[i][j][k] << "\n";
      }
      ofs1 << "\n";
    }
    ofs1.close();
  }
}
