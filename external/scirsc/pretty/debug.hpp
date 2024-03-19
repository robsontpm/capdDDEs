/*
 * pretty_debug.h
 *
 * This file is for simple and pretty formated output for debug purposes
 * I know that debugger can do some decent work on debbuging, but in
 * particular applications it may be quite helpful to output human-readable
 * debug information for end users, e.g. as in my case, when I provide
 * library for rigorous scientific computation of enclosures of true solutions
 * to delay differential equations (infinite dimensional space of silutions).
 * For a person who tries to make a computer assisted proof it is important
 * to analyze steps of the algorithm and its output to get the idea of the
 * bootlenecks and possible numerical problems. For this, I have created this simple
 * code, but it is more general so it can be used in any program.
 *
 * For examples on use see my programs for computer assisted proofs i DDEs.
 *
 *  Created on: Apr 21, 2016
 *      Author: robson
 */

#ifndef SCIRSC_PRETTY_DEBUG_HPP_
#define SCIRSC_PRETTY_DEBUG_HPP_

#include <scirsc/pretty/debug.h>

namespace scirsc{
namespace pretty{

PrettyDebug& PrettyDebug::operator<<(const Modifier& mod){
	return mod(*this);
}

template<typename AnyType>
PrettyDebug& PrettyDebug::operator<<(const AnyType& item){
	std::ostringstream os; os.precision(precision); os << item;
	std::string output = os.str();
	bool newline = (output.back() == '\n');
	std::istringstream is(output);
	std::ostringstream os_formatted;
	std::string to;

	std::getline(is, to, '\n');
	os_formatted << leftover << to;
	leftover = "\n" + preamble();
	while(std::getline(is, to, '\n')){
		os_formatted << leftover << to;
	}
	leftover = (newline ? "\n" + preamble() : "");

	to = os_formatted.str();
	if (currentLevel <= importantLevel) important << to;
	if (currentLevel <= displayLevel) all << to;

	return *this;
	}

} // namespace pretty
} // namespace scirsc

#endif /* SCIRSC_PRETTY_DEBUG_HPP_ */
