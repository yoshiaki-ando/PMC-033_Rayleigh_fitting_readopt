#ifndef ALBEDO_H_
#define ALBEDO_H_

/* アルベド考慮の定数 */
constexpr double Reflection_Coefficient { 0.3 }; /* アルベド */
constexpr int M_alpha { 30 }; /* 散乱点からの角度の分割数 */
constexpr int N_beta { 30 }; /* 散乱点を中心とする周の角度の分割数 */
constexpr double Dbeta { 2.0*M_PI / N_beta };
constexpr double R0 { Radius_of_Earth };

#endif
