/*
 * read_pmc_param.cpp
 * PMCをマニュアルに置く位置などのパラメタを読み込む
 *
 * 読み込むファイルは data/pmc.prm とする
 * 1行目に PMC分布の数
 * 2行目以降にPMC分布の数の行だけ
 * 高度[km] 水平方向シフト[km] 高度方向σz[km] 水平方向σh[km] 分布中央の数密度N0[1/cm^3] モード半径r0[nm]
 * を入れる
 *
 *  Created on: 2021/02/22
 *      Author: ando
 */
#include <fstream>
#include <string>

#include <Command.h>
#include <memory_allocate.h>

/* 単位変換 */
constexpr double KM2M { 1e3 };
constexpr double CMSQ2MSQ { 1e6 };
constexpr double NM2M { 1e-9 };

/* 各パラメタの単位変換 */
constexpr double CNVT[6] = { KM2M, KM2M, KM2M, KM2M, CMSQ2MSQ, NM2M };

void read_pmc_param(int &Number_of_PMC, double** &pmc_param){

  std::string line;
  std::ifstream ifs_param("data/pmc.prm");
  std::getline(ifs_param, line);
  AndoLab::Command cmd(line, ' ');

  Number_of_PMC = cmd.get_i(0);

  pmc_param = AndoLab::allocate_memory2d(Number_of_PMC, 6, 0.0);
  for(int i_pmc = 0; i_pmc < Number_of_PMC; i_pmc++){
    std::getline(ifs_param, line);
    AndoLab::Command cmd2(line, ' ');

    for(int i = 0; i < 6; i++){
      pmc_param[i_pmc][i] = cmd2.get_d(i) * CNVT[i];
    }
  }
  ifs_param.close();
}
