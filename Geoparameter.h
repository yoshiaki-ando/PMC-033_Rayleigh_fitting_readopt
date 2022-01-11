/*
 * Geoparameter.h
 *
 *  Created on: 2021/01/14
 *      Author: ando
 */

#ifndef GEOPARAMETER_H_
#define GEOPARAMETER_H_

#include <cmath>

constexpr double Radius_of_Earth { 6370.e3 };
constexpr double Altitude_of_GEO { 36000.e3 };
constexpr double Altitude_of_Atmosphere { 1000.e3 };
constexpr double Lower_PMC { 78.0e3 };
constexpr double Upper_PMC { 88.0e3 };

constexpr double Rgeo { Altitude_of_GEO + Radius_of_Earth };
constexpr double R_Atmosphere { Radius_of_Earth + Altitude_of_Atmosphere };

constexpr double RADIUS_OF_EARTH { 6370.0e3 };
constexpr double ALTITUDE_OF_GEO { 36000.0e3 };
constexpr double ALTITUDE_OF_ATMOSPHERE { 100.0e3 };
constexpr double Longitude_of_himawari { 140.7 };

#endif /* GEOPARAMETER_H_ */
