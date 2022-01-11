/*
 *** 変更履歴 ***
 * 2020/07/07: (fitting_rayleigh) 観測データに合わせたフィッティング
   - 線形最小二乗フィッティング、およびlogをとった非線形最小二乗フィッティング
   [ToDo] 非線形フィッティングは現在「30回反復を固定」
 *
 */
#include <iostream>
#include <fstream>
#include "pmc_simulation.h"

/* 高度から観測データの高度インデックスへの変換 */
int idx_alt(const int Number_of_Altitude,
    const double Lowest_Altitude,
    const double Step_of_Altitude,
    const double Altitude){

  int v;
  double didx = (Altitude - Lowest_Altitude) / Step_of_Altitude;

  if ( didx > Number_of_Altitude ){
    v = -1;
  } else {
    v = int( didx + 0.5 );
  }

  return v;
}

/* フィッティングのための高度インデックスの上下限 */
int idx_Alt_low_forFitting, idx_Alt_high_forFitting;

/* 非線形フィッティングのための関数 */
/* 最小二乗誤差 E のパラメタ a での微分 dE/da = F(a) */
double Fa(const double a, const double b,
          const double *yi, const double *Yi){
  double v { 0.0 };
  for(int idx_a = idx_Alt_low_forFitting; idx_a < idx_Alt_high_forFitting; idx_a++){
    v += ( std::log(a*Yi[idx_a] + b) - std::log(yi[idx_a]) )
      * Yi[idx_a] / (a*Yi[idx_a] + b);
  }

  return v;
}

/* dF(a)/da */
double dFa_da(const double a, const double b,
          const double *yi, const double *Yi){
  double v { 0.0 };
  for(int idx_a = idx_Alt_low_forFitting; idx_a < idx_Alt_high_forFitting; idx_a++){
    double t { Yi[idx_a] / (a*Yi[idx_a] + b) };
    v += t*t* ( 1.0 + std::log(a*Yi[idx_a] + b) - std::log(yi[idx_a]) );
  }
  return v;
}

void fitting_rayleigh(
    const int Number_of_angle,
    const double *alpha,
    const int Number_of_Altitude,
    const double Lowest_Altitude,
    const double Step_of_Altitude,
    double ***Calculated_Rayleigh_scattering,
    double ***Observed_data,
    int *idx_fitting_data,
    double **DL,
    double **C
    ){
      
  /*
   * ダークレベルDLの決定
   * 90-99kmの単純平均
   */
  constexpr int Lower_Altitude_for_DarkLevel { 90 };
  constexpr int Upper_Altitude_for_DarkLevel { 99 };
  constexpr int Number_of_Altitude_for_DarkLevel
  { Upper_Altitude_for_DarkLevel - Lower_Altitude_for_DarkLevel + 1 };

  for(int Lambda = 0; Lambda < Num_Lambda; Lambda++){
    for(int na = 0; na < Number_of_angle; na++){
      double avg { 0.0 };
      for(int alt = Lower_Altitude_for_DarkLevel;
          alt <= Upper_Altitude_for_DarkLevel; alt++){
        int idx_a =
            idx_alt(Number_of_Altitude, Lowest_Altitude, Step_of_Altitude,
                1e3*alt);
        std::cout << Lambda << ": " << alt << ", "
            << Observed_data[Lambda][idx_fitting_data[na]][idx_a] << std::endl;
        avg += Observed_data[Lambda][idx_fitting_data[na]][idx_a];

      }
      DL[Lambda][na] = avg / Number_of_Altitude_for_DarkLevel;
      std::cout << Lambda << ", " << DL[Lambda][na] << std::endl;
    }
  }

  /* フィッティング係数 */
  idx_Alt_low_forFitting =
      idx_alt(Number_of_Altitude, Lowest_Altitude, Step_of_Altitude,
          Alt_low_forFitting);
  idx_Alt_high_forFitting =
      idx_alt(Number_of_Altitude, Lowest_Altitude, Step_of_Altitude,
          Alt_high_forFitting);
  for(int Lambda = 0; Lambda < Num_Lambda; Lambda++){
    for(int na = 0; na < Number_of_angle; na++){

      /* 線形なフィッティング */
      double nume { 0.0 };
      double sum_Yi { 0.0 };
      for(int idx_a = idx_Alt_low_forFitting; idx_a < idx_Alt_high_forFitting; idx_a++){
//        std::cout << "> " << idx_a << " " << Calculated_Rayleigh_scattering[Lambda][na][idx_a] << " & "
//            << Observed_data[Lambda][idx_fitting_data[na]][idx_a] << std::endl;
        sum_Yi += Calculated_Rayleigh_scattering[Lambda][na][idx_a] *
            Calculated_Rayleigh_scattering[Lambda][na][idx_a];
        nume += (Observed_data[Lambda][idx_fitting_data[na]][idx_a] - DL[Lambda][na]) *
            Calculated_Rayleigh_scattering[Lambda][na][idx_a];
      }
      double alpha { nume / sum_Yi };
//      std::cout << nume << " / " << sum_Yi << std::endl;
//      std::cout << "Alpha : " << alpha << " -> ";

      /* 非線形フィッティング */
      for(int i = 0; i < 30; i++){
        alpha = alpha - Fa(alpha, DL[Lambda][na],
            Observed_data[Lambda][idx_fitting_data[na]],
            Calculated_Rayleigh_scattering[Lambda][na])
                  / dFa_da(alpha, DL[Lambda][na],
                      Observed_data[Lambda][idx_fitting_data[na]],
                      Calculated_Rayleigh_scattering[Lambda][na]);
      }
//      std::cout << alpha << std::endl;
      C[Lambda][na] = alpha;

    }
  }

}
