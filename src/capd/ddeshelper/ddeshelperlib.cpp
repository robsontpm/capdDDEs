/*
 * ddeshelperlib.cpp
 *
 *  Created on: Feb 12, 2024
 *      Author: robson
 */

#include <capd/ddeshelper/ddeshelperlib.h>

namespace capd{
namespace ddeshelper{

bool startsWith(std::string const& a, std::string const& b){
	return 0 == a.compare(0, b.length(), b);
}

bool conditionalExtractStream(std::string const& line, std::string const& key, std::istringstream& output){
	if (startsWith(line, key)){
		output.str(line.substr(key.length(), line.length())); output.clear();
		return true;
	}
	return false;
}

bool conditionalExtractStream(
		std::string const& line, std::string const& key,
		std::istringstream& output,
		std::map<std::string, std::string>& dict
){
	if (conditionalExtractStream(line, key, output)){
		dict[key] = line.substr(key.length(), line.length());
		return true;
	}
	return false;
}

std::string makeCommandLine(int argc, char** argv){
	std::ostringstream oss;
	oss << argv[0];
	for (int i = 1; i < argc; ++i)
		oss << " \"" << argv[i] << "\"";
	return oss.str();
}

bool conditionalCheckEqual(std::string const& line, std::string const& key){
	return line == key;
}

bool conditionalCheckEqual(
		std::string const& line, std::string const& key,
		std::map<std::string, std::string>& dict
){
	if (dict.find(key) == dict.end())
		dict[key] = "false";
	if (conditionalCheckEqual(line, key)){
		dict[key] = "true";
		return true;
	}
	return false;
}

void mkdir_p(std::string const& dirpath){
	std::ostringstream cmd;
	// TODO: sanitize input?
	cmd << "mkdir -p '" << dirpath << "'";
	int dumpres = system(cmd.str().c_str());
	if (dumpres) { /* just to make compiler shut up */ };
}

} // namespace ddeshelper
} // namespace capd
