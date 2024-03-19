/*
 * this declare classes to be used in the example. It is used both
 * in the programs and in the setup.cpp file to generate neccessary code
 * from the library.
 */

#ifndef DEMO_ELNINIO_NONRIG_SETUP_H_
#define DEMO_ELNINIO_NONRIG_SETUP_H_

#include "equation.h"

/** this is just a __declaration__ of a new type. The __definition__ is in the setup-nonrig.cpp file. */
typedef ElNino<double, double> Eq;
/** this is just a __declaration__ of a new type. The __definition__ is in the setup-nonrig.cpp file. */
typedef capd::ddeshelper::NonrigorousHelper<Eq> Setup;

#endif /* DEMO_ELNINIO_NONRIG_SETUP_H_ */
