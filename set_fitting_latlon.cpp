/*
 * set_fitting_latlon.cpp
 *
 *  Created on: 2021/02/15
 *      Author: ando
 */
#include "pmc_simulation.h"

void set_fitting_latlon(int *idx_latlon, double **Obsrvd_latlon,
    const int num_alpha, double **fitted_latlon){

  for(int i = 0; i < num_alpha; i++){

    /* 観測の緯度経度から、フィッティングに使う緯度経度を探す */
    for(int j_obs = 0; j_obs < Number_of_Obsrvd_data_Latitude; j_obs++){
      if ( (std::abs(Obsrvd_latlon[IDX_LATITUDE][j_obs] - fitted_latlon[IDX_LATITUDE][i]) < 0.1) && /* 緯度 */
          (std::abs(Obsrvd_latlon[IDX_LONGITUDE][j_obs] - fitted_latlon[IDX_LONGITUDE][i]) < 0.25) ){ /* 経度 */
        idx_latlon[i] = j_obs;
        break;
      }
    }

//    std::cout << i << " " << idx_latlon[i] << std::endl;
  }

}

