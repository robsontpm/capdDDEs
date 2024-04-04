#include "nonrig-setup.h"

// this will generate the code of all the classes and subroutines

template class capd::ddes::MackeyGlass<double, double>;
template class capd::ddeshelper::NonrigorousHelper<Eq, 1, capd::DMatrix, capd::DVector>;
