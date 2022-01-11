/*
 * 日付を予め与えておき、ひまわり座標を与えれば
 * 数密度、温度、気圧が得られる
 */

#include <iostream>
#include <physical_constants.h>
#include <nrlmsise-00.h>
#include "Date.h"
#include "Msis.h"

void Msis::set_date(Date date){
  /* スイッチの9が-1なら必要 */
  for(int i = 0; i < 7; i++){
    Ap.a[i] = 100.0; /* どこかのサイトで探す */
  }

  /* デフォルトのスイッチ */
  Flags.switches[0] = 0;
  for(int i = 1; i < 24; i++){
    Flags.switches[i] = 1;
  }

  Input.year = 2019;
  Input.f107A = 150.; /* F10.7フラックスの81日間平均。150.でよい */
  Input.f107 = 150.; /* 前日のF10.7フラックス。150.でよい */
  Input.ap = 4.0; /* 磁気指数。4.0でよい */
  Input.ap_a = &Ap;

  Input.doy = date.doy();
  Input.sec = date.minute_day() * 60;
}

void Msis::set_position(double lat, double longi, double alt){
  Input.g_lat = lat;
  Input.g_long = longi;
  Input.alt = alt * 1.0e-3;
  Input.lst = Input.sec/3600 + Input.g_long/15;

  gtd7(&Input, &Flags, &Output);
}

void Msis::set_position(Geocoordinate geo){
  set_position(geo.latitude(), geo.longitude(), geo.altitude());
}

double Msis::N(void){
  double n = 0.0;
  constexpr double cm_cube2m_cube = 1.0e6;

  for(int i = 0; i < 5; i++){
    n += Output.d[i];
  }

  return (n + Output.d[6] + Output.d[7])*cm_cube2m_cube;
}

double Msis::p(void){
  constexpr double kB = BOLTZMANN_CONSTANT;
  return N()*kB*T();
}

