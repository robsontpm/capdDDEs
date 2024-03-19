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

#ifndef SCIRSC_PRETTY_DEBUG_H_
#define SCIRSC_PRETTY_DEBUG_H_

#include <iostream>
#include <iomanip>
#include <stdarg.h>
#include <vector>
#include <string>
#include <sstream>

const int PRETTY_MAX_TABS = 256;

#ifndef PRETTY_IMPORTANT_LEVEL
#define PRETTY_IMPORTANT_LEVEL -1
#endif

#ifndef PRETTY_DISPLAY_LEVEL
#define PRETTY_DISPLAY_LEVEL 4
#endif

#ifndef PRETTY_TAB_SIZE
#define PRETTY_TAB_SIZE 4
#endif


namespace scirsc{
namespace pretty{

// declare PrettyDebug for Modifier to see it
class PrettyDebug;

// definition must be here to make PrettyDebug see it :)
class Modifier{
public:
	virtual PrettyDebug& operator()(PrettyDebug& debug) const = 0;
	virtual ~Modifier() {}
};

/**
 * This class is for pretty formated debug output (e.g. history of the integration
 * It can output debug data depending on current settings to two separate
 * destinations: all and important
 *
 * all is meant to capture all messages with level < displayLevel, while
 * important has its own limit of importantLevel.
 *
 * REMARK: yes, i know, this can be generalized to any number of outputs...
 *         but for now I will stick to this for naming convinience.
 *
 * You can alter the state of the object with Modifier class
 * there are several predefined modifiers in the pretty namespace.
 *
 * The class prepends [ something ] in front of input
 * The first something is set in constructor (default = ':: debug ::'
 * The size of this first something is important as it sets out the
 * width of [       ] for other items. Whenever you supply a new name
 * as something it will be stretched (with spaces) or catenated to
 * match the width of this field. So choose wisely :)
 */
class PrettyDebug{
protected:
	std::string preamble(){
		//std::ostringstream lev; lev << currentLevel;
		return "[ " /* + lev.str() + ": " */ + levelName.back() + " ] " + tabs;
	}
public:
	void leftoverPreamble() { leftover = "\n" + preamble(); }
	std::ostream& all;
	std::ostream& important;

	// tabs config
	int tabSize;
	int currentTab;
	char tabs[PRETTY_MAX_TABS];

	// error reporting level
	int importantLevel;
	int displayLevel;
	int currentLevel;

	// default precision
	int precision;

	// naming
	std::vector<std::string> levelName;

	// what to print next if << called
	std::string leftover;

	// we construct with default parameters, and the user can modify
	PrettyDebug(std::ostream& all,
				std::ostream& important,
				std::string name = ":: debug ::",
				int displayLevel = 4,
				int importantLevel = -1,
				int tabSize = 4,
				int precision = 16):
			all(all),
			important(important),
			tabSize(tabSize),
			currentTab(0),
			importantLevel(importantLevel),
			displayLevel(displayLevel),
			currentLevel(0),
			precision(precision){
		for (int i = 0; i < PRETTY_MAX_TABS; i++) tabs[i] = ' ';
		tabs[currentTab] = '\0';
		levelName.push_back(name);
		leftover = preamble();
	};

	PrettyDebug& operator<<(const Modifier& mod);

	template<typename AnyType>
	PrettyDebug& operator<<(const AnyType& item);

	~PrettyDebug(){
		all << "\n";
		important << "\n";
	}

	int setImportantLevel(int level){ int tmp = importantLevel; importantLevel = level; return tmp; }
	int setDisplayLevel(int level){ int tmp = displayLevel; displayLevel = level; return tmp; }
};

// ***************************************************************************
// * Below we have modifier classes that alters the behaviour ****************
// * of the pretty debug                                      ****************
// ***************************************************************************

// It is worth notting what is the idea, why we define __modifier and
// modifier variables after definition of a class.
// __modifier is of given specific class, while modifier is of generic Modifier
// class. This way, I can define
//    PrettyDebug& operator<<(const Modifier& mod)
// in PrettyDebug, so that ass Modifiers are used by this. Then I define
// generic << operator, which redirects everything directly to underlying streams
// otherwise I would need to override friend operator << in each of the Modifier
// subclasses or define << operator in PrettyPrint.
//
// There are two possible definitions of Modifier& modifier.
// One is direct, that is Modifier& modifier = __modifier
// this is used when modifier does not have state and/or parameters
// the second is by a function that returns Modifier&. In that
// case passing of arguments is possible. Those arguments are passed
// to __modifier via some setup() function and the reference
// to __modifier is returned. It should be used at once, as
// subsequent call to modifier() will destroy state of __modifier
// This construction is simple and no big overhaul on the space / time
// is imposed, but the modifier() must be used with caution.

// TODO: make it not inline (methods below)

/**
 * This class makes current debug level less important,
 * thus it may make it cease to be printed
 */
class LessImportant : public Modifier {
public:
	PrettyDebug& operator()(PrettyDebug& debug) const {
		debug.currentLevel++;
		return debug;
	}
} __lessImportant;
Modifier& lessImportant = __lessImportant;

/**
 * This class makes current debug level more important,
 * so if it is not printed it may start to get printed
 */
class MoreImportant : public Modifier {
public:
	PrettyDebug& operator()(PrettyDebug& debug) const {
		if (debug.currentLevel > 0)
			debug.currentLevel--;
		return debug;
	}
} __moreImportant;
Modifier& moreImportant = __moreImportant;

/**
 * Tells debug to print more spaces in front of output
 * Never exceeds maximum defined in the library
 */
class MoreTab : public Modifier {
public:
	PrettyDebug& operator()(PrettyDebug& debug) const {
		debug.tabs[debug.currentTab] = ' ';
		if (debug.currentTab + debug.tabSize < PRETTY_MAX_TABS)
			debug.currentTab += debug.tabSize;
		debug.tabs[debug.currentTab] = '\0';
		debug.leftoverPreamble();
		return debug;
	}
} __moreTab;
Modifier& moreTab = __moreTab;

/**
 * Tells debug to print less spaces in front of output
 * Never drops below 0.
 */
class LessTab : public Modifier {
public:
	PrettyDebug& operator()(PrettyDebug& debug) const {
		debug.tabs[debug.currentTab] = ' ';
		debug.currentTab -= debug.currentTab > debug.tabSize ? debug.tabSize : debug.currentTab;
		debug.tabs[debug.currentTab] = '\0';
		debug.leftoverPreamble();
		return debug;
	}
} __lessTab;
Modifier& lessTab = __lessTab;

/**
 * Put a name for a given selection, so that this name appears in front of
 * the debug output lines from now on (until popName used)
 */
class PutName : public Modifier {
public:
	std::string item;
	PrettyDebug& operator()(PrettyDebug& debug) const {
		std::string str = item.substr(0, debug.levelName[0].length());	// remove overflow characters
		while (str.length() < debug.levelName[0].length())	// if neccessary
			str += ' ';										// 	  add spaces
		debug.levelName.push_back(str); 					// remember string
		debug << "Putting: " << debug.levelName.back() << "\n";
		debug.leftoverPreamble();
		return debug;
	}
	Modifier& setup(std::string str){
		this->item = str;
		return *this;
	}
} __putName;
Modifier& putName(std::string name){
	return __putName.setup(name);
}

/**
 * Pop current name and revert to the previous name
 * Allways stops at default name defined in debug contructor
 */
class popName : public Modifier {
public:
	PrettyDebug& operator()(PrettyDebug& debug) const {
		debug << "Popping: " << debug.levelName.back() << "\n";
		if (debug.levelName.size() > 1)
			debug.levelName.pop_back();
		debug.leftoverPreamble();
		return debug;
	}
} __popName;
Modifier& popName = __popName;

/**
 * Combines several modifiers into one
 */
class ModCombined : public Modifier {
public:
	typedef Modifier* ListType;
	PrettyDebug& operator()(PrettyDebug& debug) const {
		for (auto it = mods.begin(); it != mods.end(); it++)
			(**it)(debug);
		return debug;
	}
	std::vector<ListType> mods;

	ModCombined(int n, ...){
	    va_list ap;
	    va_start(ap, n);
		for (int i = 0; i < n; i++)
			mods.push_back(va_arg(ap, ListType));
		va_end(ap);
	}
	template<typename IteratorType>
	ModCombined(IteratorType begin, IteratorType end){
		for (IteratorType it = begin; it != end; it++)
			mods.push_back(*it);
	}
	ModCombined& append(Modifier* mod){
		mods.push_back(mod);
		return *this;
	}
} __descend(2, &moreTab, &lessImportant), __ascend(2, &lessTab, &moreImportant);
Modifier& descend = __descend;
Modifier& ascend = __ascend;

// for enter and exit modifiers
PutName __enterSetName; // auxiliary
ModCombined __enter(3, &moreTab, &lessImportant, &__enterSetName);
ModCombined __leave(3, &lessTab, &moreImportant, &__popName);
Modifier& enter(std::string functionName){
	__enterSetName.setup(functionName);
	return __enter;
}
Modifier& leave = __leave;

} // namespace pretty
} // namespace scirsc

#ifndef PRETTY_GLOBAL_DEBUG
	#define PRETTY_GLOBAL_DEBUG " :: global :: "
#else
	#define PRETTY_GLOBAL_DEBUG_ON
#endif

#ifdef PRETTY_GLOBAL_DEBUG_ON
	scirsc::pretty::PrettyDebug __globpd(std::cerr, std::cout, PRETTY_GLOBAL_DEBUG, PRETTY_DISPLAY_LEVEL, PRETTY_IMPORTANT_LEVEL, PRETTY_TAB_SIZE);
	#define glodeb(data) __globpd << data << "\n";
#else
	#define glodeb(data)
#endif

#endif /* SCIRSC_PRETTY_DEBUG_H_ */
