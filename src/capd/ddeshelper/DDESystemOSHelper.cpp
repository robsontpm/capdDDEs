/*
 * DDESystemOSHelper.cpp
 */

#include <capd/ddeshelper/DDEHelperCommon.h>
#include <capd/ddeshelper/DDESystemOSHelper.h>

namespace capd {
namespace ddeshelper {

/** */
template<>
bool ArgumentParser::parse(const std::string &param, capd::interval &out,
		const char *const help_s) {
	startNewParam(param);
	if (help_s) (*this) << help_s << "\n";
	(*this) << "NOTE 1: loaded argument will be rigorous, even if the ends of the given interval are not representable numbers!" << "\n";
	(*this) << "NOTE 2: it works with normal 'double' values too, i.e. 1.23 will load a small interval around 1.23!" << "\n";
	(*this) << "[default: " << std::setprecision(16) << out << "]";
	bool found = false;
	std::string rep;
	interval v = 0;
	for (ArgcType i = 0; i < argc; ++i){
		auto test = conditionalExtractValue(argv[i], param, rep);
		if (!test) continue;
		try { // process the data:
			if (rep[0] == '['){ // we have an interval
				auto coma_pos = rep.find(',');
				if (coma_pos == std::string::npos) out = interval(rep, "(DOES NOT HAVE COMMA!)");
				else {
					auto comma = rep.begin() + coma_pos;
					std::istringstream ileft(std::string(rep.begin()+1, comma));	// + to skip [
					std::string left; ileft >> left;    // remove extra spaces if any.
					std::istringstream iright(std::string(comma + 1, rep.end()-1)); // -1 for closing ]
					std::string right; iright >> right; // remove extra spaces if any.
					v = interval(left, right);
				}
			}else{ // we have double (probably);
				v = interval(rep, rep);
			}
			test = false;
		} catch (std::exception& ex){
			test = false;
		}
		found = found or test;
		out = v; // it is safe now to pass it to output. All conversions succeded.
	}
	return found;
}

} // namespace ddeshelper
} // namespace capd
