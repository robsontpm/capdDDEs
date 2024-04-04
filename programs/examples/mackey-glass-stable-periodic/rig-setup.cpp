#include "rig-setup.h"

// this will generate the code of all the classes and subroutines

template class capd::intervals::Interval<double, capd::rounding::DoubleRounding>;
template class capd::ddes::MackeyGlass<Interval, Interval>;
template class capd::ddeshelper::RigorousHelper<Eq, 1>;
template class capd::ddeshelper::DDECompareHelper<RigSetup::Vector>;
