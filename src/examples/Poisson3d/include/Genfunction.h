#ifndef DENDRO_5_0_GENFUNCTIONS_H
#define DENDRO_5_0_GENFUNCTIONS_H

#include "functional"
#include <math.h>

void set_grid_funs(std::function<double(double)> gx2x, std::function<double(double)> gy2y, std::function<double(double)> gz2z);

void genfunction_0(double x, double y, double z, double* var);
void genfunction_1(double x, double y, double z, double* var);
#endif //DENDRO_5_0_GENFUNCTIONS_H
