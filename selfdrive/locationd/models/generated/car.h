#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_8021895057910774678);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_162742718340394609);
void car_H_mod_fun(double *state, double *out_4674772262529250125);
void car_f_fun(double *state, double dt, double *out_8882169823147127836);
void car_F_fun(double *state, double dt, double *out_342562984584705687);
void car_h_25(double *state, double *unused, double *out_7899195939648737914);
void car_H_25(double *state, double *unused, double *out_5948428563763182938);
void car_h_24(double *state, double *unused, double *out_1654583961423930527);
void car_H_24(double *state, double *unused, double *out_8940389911144045542);
void car_h_30(double *state, double *unused, double *out_8174390001933243803);
void car_H_30(double *state, double *unused, double *out_8466761522270431565);
void car_h_26(double *state, double *unused, double *out_1382147028301512024);
void car_H_26(double *state, double *unused, double *out_2206925244889126714);
void car_h_27(double *state, double *unused, double *out_1500154422234888768);
void car_H_27(double *state, double *unused, double *out_7756388480255176834);
void car_h_29(double *state, double *unused, double *out_5270680804115462168);
void car_H_29(double *state, double *unused, double *out_8976992866584823749);
void car_h_28(double *state, double *unused, double *out_6501827689895734413);
void car_H_28(double *state, double *unused, double *out_3894593849515293175);
void car_h_31(double *state, double *unused, double *out_4261746411394715564);
void car_H_31(double *state, double *unused, double *out_5979074525640143366);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}