#include "rig-setup.h"

// this will generate the code of all the classes and subroutines

template class capd::intervals::Interval<double, capd::rounding::DoubleRounding>;
template class capd::ddes::MackeyGlass<Interval, Interval>;
template class capd::ddeshelper::RigorousHelper<IEq, 1>;
template class capd::ddeshelper::DDECompareHelper<IDDEs::Vector>;
