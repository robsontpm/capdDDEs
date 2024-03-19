/*
 * This should instantiate the library (generate and compile all code from the classes below).
 * It should make compilation of other files faster (if there are no changes in ElNino (==Eq) class)
 */

#include "setup-nonrig.h"

/** this is __definition__ of a new type, with all methods and such generated automatically. */
template class ElNino<double, double>;

/** this is __definition__ of a new type, with all methods and such generated automatically. */
template class capd::ddeshelper::NonrigorousHelper<Eq>;

///** specialized version of the template from setup-nonrig.h. It is a concrete function, so it must be in .cpp. */
//template<> inline double get_pi(){ return 3.14159265358979323846; }
