/*
 * DDECommon.cpp
 *
 *  Created on: Feb 12, 2024
 *      Author: robson
 */

#include <capd/ddes/DDECommon.h>

namespace capd{
namespace ddes{

int closestSmallerInt(capd::interval const & value){ return closestSmallerInt(value.leftBound()); }

void helper_dump_line(std::istream & in){
	std::string dump; std::getline(in, dump);
}

void helper_dump_badge(std::istream & in){
	std::string dump; in >> dump;
}

double ecloseStep(double h){ return h; }

std::string showEnclosedInterval(double a, double b){ std::ostringstream oss; oss << "[" << a << ", " << b << ")"; return oss.str(); }


} // namespace ddes
} // namespace capd
