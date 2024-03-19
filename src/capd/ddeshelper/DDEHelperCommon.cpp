/*
 * DDEHelperCommon.cpp
 *
 *  Created on: Feb 12, 2024
 *      Author: robson
 */

#include <capd/ddeshelper/DDEHelperCommon.h>

namespace capd {
namespace ddeshelper {

void runSystemCommand(std::string cmd){
	// TODO: (FUTURE) - make sure that everything works... and system is available on all systems (should be, according to docs)
	// TODO: (NOT URGENT) - create define flag to turn this on or off. Might be some problem, as other routines (like drawing) might use this...
	// TODO: (NOT URGENT) - in that situation maybe just throw exception when using routines that depend on this??
	if (system(cmd.c_str())) { /* just to shut up the compiler warnings... */ }
}


} // namespace ddeshelper
} // namespace capd
