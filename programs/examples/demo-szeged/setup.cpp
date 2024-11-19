/*
 * This should instantiate the library (generate and compile all code from the classes below).
 * It should make compilation of other files faster (if there are no changes in ElNino (==Eq) class)
 */

#include "setup.h"

template class DemoSzeged<double, double>;
template class DemoSzeged<capd::interval, capd::interval>;
template class capd::ddes::NonrigorousSetup<DEq>;
template class capd::ddes::RigorousSetup<IEq>;

