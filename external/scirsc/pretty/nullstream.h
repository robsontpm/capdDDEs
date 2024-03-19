/*
 * nullstream.h
 *
 *  Created on: 22-01-2013
 *      Author: Robson
 */

#ifndef SCIRSC_PRETTY_NULLSTREAM_H_
#define SCIRSC_PRETTY_NULLSTREAM_H_

namespace scirsc{
namespace pretty{

/**
 * this class ignores all manipulators and input
 *
 * It is safe to use this to disable output for example in
 * PrettyDebug in all or important...
 */
class nullstream{
public:
	template<typename T>
	nullstream& operator<<(const T& item){
		return *this;
	}
};

} // pretty
} // scirsc

#endif /* SCIRSC_PRETTY_NULLSTREAM_H_ */
