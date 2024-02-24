#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8021895057910774678) {
   out_8021895057910774678[0] = delta_x[0] + nom_x[0];
   out_8021895057910774678[1] = delta_x[1] + nom_x[1];
   out_8021895057910774678[2] = delta_x[2] + nom_x[2];
   out_8021895057910774678[3] = delta_x[3] + nom_x[3];
   out_8021895057910774678[4] = delta_x[4] + nom_x[4];
   out_8021895057910774678[5] = delta_x[5] + nom_x[5];
   out_8021895057910774678[6] = delta_x[6] + nom_x[6];
   out_8021895057910774678[7] = delta_x[7] + nom_x[7];
   out_8021895057910774678[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_162742718340394609) {
   out_162742718340394609[0] = -nom_x[0] + true_x[0];
   out_162742718340394609[1] = -nom_x[1] + true_x[1];
   out_162742718340394609[2] = -nom_x[2] + true_x[2];
   out_162742718340394609[3] = -nom_x[3] + true_x[3];
   out_162742718340394609[4] = -nom_x[4] + true_x[4];
   out_162742718340394609[5] = -nom_x[5] + true_x[5];
   out_162742718340394609[6] = -nom_x[6] + true_x[6];
   out_162742718340394609[7] = -nom_x[7] + true_x[7];
   out_162742718340394609[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4674772262529250125) {
   out_4674772262529250125[0] = 1.0;
   out_4674772262529250125[1] = 0;
   out_4674772262529250125[2] = 0;
   out_4674772262529250125[3] = 0;
   out_4674772262529250125[4] = 0;
   out_4674772262529250125[5] = 0;
   out_4674772262529250125[6] = 0;
   out_4674772262529250125[7] = 0;
   out_4674772262529250125[8] = 0;
   out_4674772262529250125[9] = 0;
   out_4674772262529250125[10] = 1.0;
   out_4674772262529250125[11] = 0;
   out_4674772262529250125[12] = 0;
   out_4674772262529250125[13] = 0;
   out_4674772262529250125[14] = 0;
   out_4674772262529250125[15] = 0;
   out_4674772262529250125[16] = 0;
   out_4674772262529250125[17] = 0;
   out_4674772262529250125[18] = 0;
   out_4674772262529250125[19] = 0;
   out_4674772262529250125[20] = 1.0;
   out_4674772262529250125[21] = 0;
   out_4674772262529250125[22] = 0;
   out_4674772262529250125[23] = 0;
   out_4674772262529250125[24] = 0;
   out_4674772262529250125[25] = 0;
   out_4674772262529250125[26] = 0;
   out_4674772262529250125[27] = 0;
   out_4674772262529250125[28] = 0;
   out_4674772262529250125[29] = 0;
   out_4674772262529250125[30] = 1.0;
   out_4674772262529250125[31] = 0;
   out_4674772262529250125[32] = 0;
   out_4674772262529250125[33] = 0;
   out_4674772262529250125[34] = 0;
   out_4674772262529250125[35] = 0;
   out_4674772262529250125[36] = 0;
   out_4674772262529250125[37] = 0;
   out_4674772262529250125[38] = 0;
   out_4674772262529250125[39] = 0;
   out_4674772262529250125[40] = 1.0;
   out_4674772262529250125[41] = 0;
   out_4674772262529250125[42] = 0;
   out_4674772262529250125[43] = 0;
   out_4674772262529250125[44] = 0;
   out_4674772262529250125[45] = 0;
   out_4674772262529250125[46] = 0;
   out_4674772262529250125[47] = 0;
   out_4674772262529250125[48] = 0;
   out_4674772262529250125[49] = 0;
   out_4674772262529250125[50] = 1.0;
   out_4674772262529250125[51] = 0;
   out_4674772262529250125[52] = 0;
   out_4674772262529250125[53] = 0;
   out_4674772262529250125[54] = 0;
   out_4674772262529250125[55] = 0;
   out_4674772262529250125[56] = 0;
   out_4674772262529250125[57] = 0;
   out_4674772262529250125[58] = 0;
   out_4674772262529250125[59] = 0;
   out_4674772262529250125[60] = 1.0;
   out_4674772262529250125[61] = 0;
   out_4674772262529250125[62] = 0;
   out_4674772262529250125[63] = 0;
   out_4674772262529250125[64] = 0;
   out_4674772262529250125[65] = 0;
   out_4674772262529250125[66] = 0;
   out_4674772262529250125[67] = 0;
   out_4674772262529250125[68] = 0;
   out_4674772262529250125[69] = 0;
   out_4674772262529250125[70] = 1.0;
   out_4674772262529250125[71] = 0;
   out_4674772262529250125[72] = 0;
   out_4674772262529250125[73] = 0;
   out_4674772262529250125[74] = 0;
   out_4674772262529250125[75] = 0;
   out_4674772262529250125[76] = 0;
   out_4674772262529250125[77] = 0;
   out_4674772262529250125[78] = 0;
   out_4674772262529250125[79] = 0;
   out_4674772262529250125[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_8882169823147127836) {
   out_8882169823147127836[0] = state[0];
   out_8882169823147127836[1] = state[1];
   out_8882169823147127836[2] = state[2];
   out_8882169823147127836[3] = state[3];
   out_8882169823147127836[4] = state[4];
   out_8882169823147127836[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8882169823147127836[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8882169823147127836[7] = state[7];
   out_8882169823147127836[8] = state[8];
}
void F_fun(double *state, double dt, double *out_342562984584705687) {
   out_342562984584705687[0] = 1;
   out_342562984584705687[1] = 0;
   out_342562984584705687[2] = 0;
   out_342562984584705687[3] = 0;
   out_342562984584705687[4] = 0;
   out_342562984584705687[5] = 0;
   out_342562984584705687[6] = 0;
   out_342562984584705687[7] = 0;
   out_342562984584705687[8] = 0;
   out_342562984584705687[9] = 0;
   out_342562984584705687[10] = 1;
   out_342562984584705687[11] = 0;
   out_342562984584705687[12] = 0;
   out_342562984584705687[13] = 0;
   out_342562984584705687[14] = 0;
   out_342562984584705687[15] = 0;
   out_342562984584705687[16] = 0;
   out_342562984584705687[17] = 0;
   out_342562984584705687[18] = 0;
   out_342562984584705687[19] = 0;
   out_342562984584705687[20] = 1;
   out_342562984584705687[21] = 0;
   out_342562984584705687[22] = 0;
   out_342562984584705687[23] = 0;
   out_342562984584705687[24] = 0;
   out_342562984584705687[25] = 0;
   out_342562984584705687[26] = 0;
   out_342562984584705687[27] = 0;
   out_342562984584705687[28] = 0;
   out_342562984584705687[29] = 0;
   out_342562984584705687[30] = 1;
   out_342562984584705687[31] = 0;
   out_342562984584705687[32] = 0;
   out_342562984584705687[33] = 0;
   out_342562984584705687[34] = 0;
   out_342562984584705687[35] = 0;
   out_342562984584705687[36] = 0;
   out_342562984584705687[37] = 0;
   out_342562984584705687[38] = 0;
   out_342562984584705687[39] = 0;
   out_342562984584705687[40] = 1;
   out_342562984584705687[41] = 0;
   out_342562984584705687[42] = 0;
   out_342562984584705687[43] = 0;
   out_342562984584705687[44] = 0;
   out_342562984584705687[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_342562984584705687[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_342562984584705687[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_342562984584705687[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_342562984584705687[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_342562984584705687[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_342562984584705687[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_342562984584705687[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_342562984584705687[53] = -9.8000000000000007*dt;
   out_342562984584705687[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_342562984584705687[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_342562984584705687[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_342562984584705687[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_342562984584705687[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_342562984584705687[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_342562984584705687[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_342562984584705687[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_342562984584705687[62] = 0;
   out_342562984584705687[63] = 0;
   out_342562984584705687[64] = 0;
   out_342562984584705687[65] = 0;
   out_342562984584705687[66] = 0;
   out_342562984584705687[67] = 0;
   out_342562984584705687[68] = 0;
   out_342562984584705687[69] = 0;
   out_342562984584705687[70] = 1;
   out_342562984584705687[71] = 0;
   out_342562984584705687[72] = 0;
   out_342562984584705687[73] = 0;
   out_342562984584705687[74] = 0;
   out_342562984584705687[75] = 0;
   out_342562984584705687[76] = 0;
   out_342562984584705687[77] = 0;
   out_342562984584705687[78] = 0;
   out_342562984584705687[79] = 0;
   out_342562984584705687[80] = 1;
}
void h_25(double *state, double *unused, double *out_7899195939648737914) {
   out_7899195939648737914[0] = state[6];
}
void H_25(double *state, double *unused, double *out_5948428563763182938) {
   out_5948428563763182938[0] = 0;
   out_5948428563763182938[1] = 0;
   out_5948428563763182938[2] = 0;
   out_5948428563763182938[3] = 0;
   out_5948428563763182938[4] = 0;
   out_5948428563763182938[5] = 0;
   out_5948428563763182938[6] = 1;
   out_5948428563763182938[7] = 0;
   out_5948428563763182938[8] = 0;
}
void h_24(double *state, double *unused, double *out_1654583961423930527) {
   out_1654583961423930527[0] = state[4];
   out_1654583961423930527[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8940389911144045542) {
   out_8940389911144045542[0] = 0;
   out_8940389911144045542[1] = 0;
   out_8940389911144045542[2] = 0;
   out_8940389911144045542[3] = 0;
   out_8940389911144045542[4] = 1;
   out_8940389911144045542[5] = 0;
   out_8940389911144045542[6] = 0;
   out_8940389911144045542[7] = 0;
   out_8940389911144045542[8] = 0;
   out_8940389911144045542[9] = 0;
   out_8940389911144045542[10] = 0;
   out_8940389911144045542[11] = 0;
   out_8940389911144045542[12] = 0;
   out_8940389911144045542[13] = 0;
   out_8940389911144045542[14] = 1;
   out_8940389911144045542[15] = 0;
   out_8940389911144045542[16] = 0;
   out_8940389911144045542[17] = 0;
}
void h_30(double *state, double *unused, double *out_8174390001933243803) {
   out_8174390001933243803[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8466761522270431565) {
   out_8466761522270431565[0] = 0;
   out_8466761522270431565[1] = 0;
   out_8466761522270431565[2] = 0;
   out_8466761522270431565[3] = 0;
   out_8466761522270431565[4] = 1;
   out_8466761522270431565[5] = 0;
   out_8466761522270431565[6] = 0;
   out_8466761522270431565[7] = 0;
   out_8466761522270431565[8] = 0;
}
void h_26(double *state, double *unused, double *out_1382147028301512024) {
   out_1382147028301512024[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2206925244889126714) {
   out_2206925244889126714[0] = 0;
   out_2206925244889126714[1] = 0;
   out_2206925244889126714[2] = 0;
   out_2206925244889126714[3] = 0;
   out_2206925244889126714[4] = 0;
   out_2206925244889126714[5] = 0;
   out_2206925244889126714[6] = 0;
   out_2206925244889126714[7] = 1;
   out_2206925244889126714[8] = 0;
}
void h_27(double *state, double *unused, double *out_1500154422234888768) {
   out_1500154422234888768[0] = state[3];
}
void H_27(double *state, double *unused, double *out_7756388480255176834) {
   out_7756388480255176834[0] = 0;
   out_7756388480255176834[1] = 0;
   out_7756388480255176834[2] = 0;
   out_7756388480255176834[3] = 1;
   out_7756388480255176834[4] = 0;
   out_7756388480255176834[5] = 0;
   out_7756388480255176834[6] = 0;
   out_7756388480255176834[7] = 0;
   out_7756388480255176834[8] = 0;
}
void h_29(double *state, double *unused, double *out_5270680804115462168) {
   out_5270680804115462168[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8976992866584823749) {
   out_8976992866584823749[0] = 0;
   out_8976992866584823749[1] = 1;
   out_8976992866584823749[2] = 0;
   out_8976992866584823749[3] = 0;
   out_8976992866584823749[4] = 0;
   out_8976992866584823749[5] = 0;
   out_8976992866584823749[6] = 0;
   out_8976992866584823749[7] = 0;
   out_8976992866584823749[8] = 0;
}
void h_28(double *state, double *unused, double *out_6501827689895734413) {
   out_6501827689895734413[0] = state[0];
}
void H_28(double *state, double *unused, double *out_3894593849515293175) {
   out_3894593849515293175[0] = 1;
   out_3894593849515293175[1] = 0;
   out_3894593849515293175[2] = 0;
   out_3894593849515293175[3] = 0;
   out_3894593849515293175[4] = 0;
   out_3894593849515293175[5] = 0;
   out_3894593849515293175[6] = 0;
   out_3894593849515293175[7] = 0;
   out_3894593849515293175[8] = 0;
}
void h_31(double *state, double *unused, double *out_4261746411394715564) {
   out_4261746411394715564[0] = state[8];
}
void H_31(double *state, double *unused, double *out_5979074525640143366) {
   out_5979074525640143366[0] = 0;
   out_5979074525640143366[1] = 0;
   out_5979074525640143366[2] = 0;
   out_5979074525640143366[3] = 0;
   out_5979074525640143366[4] = 0;
   out_5979074525640143366[5] = 0;
   out_5979074525640143366[6] = 0;
   out_5979074525640143366[7] = 0;
   out_5979074525640143366[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_8021895057910774678) {
  err_fun(nom_x, delta_x, out_8021895057910774678);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_162742718340394609) {
  inv_err_fun(nom_x, true_x, out_162742718340394609);
}
void car_H_mod_fun(double *state, double *out_4674772262529250125) {
  H_mod_fun(state, out_4674772262529250125);
}
void car_f_fun(double *state, double dt, double *out_8882169823147127836) {
  f_fun(state,  dt, out_8882169823147127836);
}
void car_F_fun(double *state, double dt, double *out_342562984584705687) {
  F_fun(state,  dt, out_342562984584705687);
}
void car_h_25(double *state, double *unused, double *out_7899195939648737914) {
  h_25(state, unused, out_7899195939648737914);
}
void car_H_25(double *state, double *unused, double *out_5948428563763182938) {
  H_25(state, unused, out_5948428563763182938);
}
void car_h_24(double *state, double *unused, double *out_1654583961423930527) {
  h_24(state, unused, out_1654583961423930527);
}
void car_H_24(double *state, double *unused, double *out_8940389911144045542) {
  H_24(state, unused, out_8940389911144045542);
}
void car_h_30(double *state, double *unused, double *out_8174390001933243803) {
  h_30(state, unused, out_8174390001933243803);
}
void car_H_30(double *state, double *unused, double *out_8466761522270431565) {
  H_30(state, unused, out_8466761522270431565);
}
void car_h_26(double *state, double *unused, double *out_1382147028301512024) {
  h_26(state, unused, out_1382147028301512024);
}
void car_H_26(double *state, double *unused, double *out_2206925244889126714) {
  H_26(state, unused, out_2206925244889126714);
}
void car_h_27(double *state, double *unused, double *out_1500154422234888768) {
  h_27(state, unused, out_1500154422234888768);
}
void car_H_27(double *state, double *unused, double *out_7756388480255176834) {
  H_27(state, unused, out_7756388480255176834);
}
void car_h_29(double *state, double *unused, double *out_5270680804115462168) {
  h_29(state, unused, out_5270680804115462168);
}
void car_H_29(double *state, double *unused, double *out_8976992866584823749) {
  H_29(state, unused, out_8976992866584823749);
}
void car_h_28(double *state, double *unused, double *out_6501827689895734413) {
  h_28(state, unused, out_6501827689895734413);
}
void car_H_28(double *state, double *unused, double *out_3894593849515293175) {
  H_28(state, unused, out_3894593849515293175);
}
void car_h_31(double *state, double *unused, double *out_4261746411394715564) {
  h_31(state, unused, out_4261746411394715564);
}
void car_H_31(double *state, double *unused, double *out_5979074525640143366) {
  H_31(state, unused, out_5979074525640143366);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
