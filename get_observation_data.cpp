#include <iostream>
#include <fstream>
#include <string>

#include <Command.h>
#include "pmc_simulation.h"

void get_observation_data(
    const int N_alt,
    const double Altitude_min,
    const double dAlt,
    double ***intensity,
    double **latlon
    ){

  std::string line;

  /* 観測データにはあるが、使用しない部分のオフセット */
  constexpr int Offset_for_not_altitude_data { 2 }; /* 緯度・経度情報 */
  const int Offset_altitude /* 0〜Alt_min [km]までは使用しない */
  { Offset_for_not_altitude_data +
    int((Altitude_min - Obsrvd_data_Lowest_Altitude)/dAlt) };

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    std::ifstream ifs( ("data/h08_b0" + std::to_string(j_lambda+1) + "_s01s02_20160709_210000.txt").c_str() );
    std::getline(ifs, line); /* 1行目はただ height */
    std::getline(ifs, line); /* 2行目は lat lon 0 〜 100 */

    for(int i_alpha = 0; i_alpha < Number_of_Obsrvd_data_Latitude; i_alpha++){
      std::getline(ifs, line); /* 3行目以降は lat lon 0 〜 100の輝度 */
      AndoLab::Command cmd(line, ' ');

      /* 緯度・経度情報を記憶 */
      latlon[IDX_LATITUDE][Number_of_Obsrvd_data_Latitude - 1 - i_alpha] = cmd.get_d(0);
      latlon[IDX_LONGITUDE][Number_of_Obsrvd_data_Latitude - 1 - i_alpha] = cmd.get_d(1);

      for(int k = 0; k < N_alt; k++){
        /* 角度αが負はアメリカ側、正がヨーロッパ側、観測データは逆なので、逆転して保存 */
        intensity[j_lambda][Number_of_Obsrvd_data_Latitude - 1 - i_alpha][k] =
            cmd.get_d(Offset_altitude + k);
//        std::cout << k << " " << Offset_altitude + k << " " << cmd.get_d(Offset_altitude + k) << std::endl;
      }
    }
    ifs.close();
  }


}
