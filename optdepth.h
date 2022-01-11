/*
 * optdepth.h
 *
 * 光学的深さを再計算しないように、計算したものを保存
 * または、保存したものを読み込む
 *
 *  Created on: 2022/01/09
 *      Author: ando
 */

#ifndef OPTDEPTH_H_
#define OPTDEPTH_H_

#include <iostream>
#include <fstream>
#include <string>

extern std::ofstream ofs_optdepth;
extern std::ifstream ifs_optdepth;
extern std::string str_optdepth;


#endif /* OPTDEPTH_H_ */
