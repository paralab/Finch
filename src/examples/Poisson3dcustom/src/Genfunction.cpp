//This file was generated by Finch.

/*

*/
std::function<void(double,double,double,double*)> genfunction_0 = [gridX_to_X, gridY_to_Y, gridZ_to_Z](const double x,const double y,const double z,double *var){
var[0] = -14 * M_PI * M_PI * sin(3 * M_PI * gridX_to_X(x)) * sin(2 * M_PI * gridY_to_Y(y)) * sin(M_PI * gridZ_to_Z(z));
};