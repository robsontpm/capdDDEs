/*
 * This should instantiate the library (generate and compile all code from the classes below).
 * It should make compilation of other files faster (if there are no changes in ExemplaryDDE (==Eq) class)
 */

#include "setup.h"

template class ExemplaryDDE<double, double>;
template class ExemplaryDDE<capd::interval, capd::interval>;

