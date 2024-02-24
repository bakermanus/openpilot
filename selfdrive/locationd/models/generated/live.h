#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_5615539772670674435);
void live_err_fun(double *nom_x, double *delta_x, double *out_3184297168680991965);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_1346894790680346190);
void live_H_mod_fun(double *state, double *out_3305990753123808482);
void live_f_fun(double *state, double dt, double *out_5565335988102264817);
void live_F_fun(double *state, double dt, double *out_6137194908109777179);
void live_h_4(double *state, double *unused, double *out_4631659118297702821);
void live_H_4(double *state, double *unused, double *out_2040864674671359739);
void live_h_9(double *state, double *unused, double *out_4601435483314280008);
void live_H_9(double *state, double *unused, double *out_9118660463773744407);
void live_h_10(double *state, double *unused, double *out_916556768011054270);
void live_H_10(double *state, double *unused, double *out_7557164162780314037);
void live_h_12(double *state, double *unused, double *out_4510911896990896196);
void live_H_12(double *state, double *unused, double *out_7060321082703321534);
void live_h_35(double *state, double *unused, double *out_6073292320071698116);
void live_H_35(double *state, double *unused, double *out_5407526732043967115);
void live_h_32(double *state, double *unused, double *out_8888027542773416496);
void live_H_32(double *state, double *unused, double *out_191883680823115118);
void live_h_13(double *state, double *unused, double *out_6130404106778203215);
void live_H_13(double *state, double *unused, double *out_5152835145042968857);
void live_h_14(double *state, double *unused, double *out_4601435483314280008);
void live_H_14(double *state, double *unused, double *out_9118660463773744407);
void live_h_33(double *state, double *unused, double *out_4134089893039182451);
void live_H_33(double *state, double *unused, double *out_8558083736682824719);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}