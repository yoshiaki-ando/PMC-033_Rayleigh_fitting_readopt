/*
 * msis.h
 *
 *  Created on: 2019/07/09
 *      Author: ando
 */

#ifndef MSIS_H_
#define MSIS_H_

#include <nrlmsise-00.h>
#include "Date.h"
#include "Geocoordinate.h"

class Msis{
private:
  struct ap_array Ap; /* Ap指数 */
  struct nrlmsise_input Input; /* 入力データ */
  struct nrlmsise_flags Flags; /* 各種スイッチ */
  struct nrlmsise_output Output;  /* 出力データ */
  double pLatitude, pLongitude, pAltitude; /* 緯度、経度、高度 */

public:
  Msis(void){};
  Msis(Date d){ set_date(d); };
  void set_date(Date);
  double latitude(void){ return pLatitude; }
  double longitude(void){ return pLongitude; }
  double altitude(void){ return pAltitude; }

  void set_position(
      double latitude,
      double longitude,
      double altitude
      );
  void set_position(Geocoordinate);

  double T(void){ return Output.t[1]; } /* [K] 温度 */
  double N(void);                        /* [m^{-3}] 数密度 */
  double p(void);                        /* [Pa] 気圧 */
};



#endif /* MSIS_H_ */
