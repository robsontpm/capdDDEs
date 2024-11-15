/*
 * ddeshelperlib.h
 *
 * This file constitutes part of DDEs rigorous integration framework developed
 * in PhD Thesis under supervision of prof. Piotr Zgliczynski:
 *
 * 		"Rigorous Integration of Delay Differential Equations", Jagiellonian University, 2015
 *
 * When using in scientific work please consider citing my articles (preffered),
 * PhD thesis and/or my webpage. For the publications describing the source code
 * please refer to http://scirsc.org/p/papers. Most notable paper up to date is:
 *
 *     Szczelina, R.; Zgliczyński, P.; "Algorithm for rigorous integration of
 *     Delay Differential Equations and the computer-assisted proof of periodic
 *     orbits in the Mackey-Glass equation", http://dx.doi.org/10.1007/s10208-017-9369-5
 *     Foundations of Computational Mathematics (2018), Vol. 18, Iss 6, Pages 1299--1332
 *
 * This work would not be possible without aid and expertise of people involved in
 * CAPD developing library (Computer Assisted Proofs in Dynamics).
 * Please refer to http://capd.ii.uj.edu.pl and consider citing also this library
 * when using those codes in any scientific work.
 *
 * Author: Robert Szczelina, PhD
 * Faculty of Mathematics and Computer Science, Jagiellonian University AND
 * (former) Małopolska Center of Biotechnology, Jagiellonian University
 * email: 	robert.szczelina@uj.edu.pl
 * www: 	scirsc.org
 *
 * This source code is provided under GNU GPL license
 * (v.2 or whatever compatible with CAPD license)
 */

#ifndef _CAPD_DDES_DDESHELPERSLIB_H_
#define _CAPD_DDES_DDESHELPERSLIB_H_

//#ifdef _CAPD_MPCAPD_H_
//#include <capd/ddeshelper/DDEMultiPrecHelper.h>
//#endif

#include <capd/ddeshelper/DDEHelperCommon.h>
#include <capd/ddeshelper/DDEHelperDrawing.hpp>
#include <capd/ddeshelper/DDECompareHelper.h>
#include <capd/ddeshelper/DDECoordinateFrameHelper.hpp>
// user should decide by himself
//#include <capd/ddeshelper/DDEHelperNonrigorous.hpp>
//#include <capd/ddeshelper/DDEHelperRigorous.hpp>


// TODO: (NOT IMPORTANT): Move to some other file.
namespace capd{
namespace ddeshelper{

std::string makeCommandLine(int argc, char** argv);

////////////////////////////////
template<typename VT>
void strToVector(std::string &data, VT& output){
	if (data[0] == '{'){
		std::istringstream iss(data);
		if (data[1] == '['){ 	// Interval
			iss >> output;
		}else{					// Double
			VT tmp(output.dimension());
			iss >> tmp;
			output = VT(tmp);
		}
	}else{
		capd::ddeshelper::readBinary(data, output);
	}
}
////////////////////////////////

bool startsWith(const std::string &a, const std::string &b);

bool conditionalExtractStream(const std::string &line, const std::string &key,
		std::istringstream &output);

bool conditionalExtractStream(
		const std::string &line, const std::string &key,
		std::istringstream& output,
		std::map<std::string, std::string>& dict
);

////////////////////////////////
template<typename T>
bool conditionalExtractValue(const std::string &line, const std::string &key,
		T &output) {
	if (startsWith(line, key)){
		std::istringstream iss(line.substr(key.length(), line.length()));
		iss >> output;
		return true;
	}
	return false;
}
////////////////////////////////

////////////////////////////////
template<typename T>
bool conditionalExtractValue(
		const std::string &line, const std::string &key,
		T& output, std::map<std::string, std::string>& dict
){
	// save default value if not yet stored!
	if (dict.find(key) == dict.end()){
		std::ostringstream data; data.precision(16);
		data << output;
		dict[key] = data.str();
	}
	if (conditionalExtractValue(line, key, output)){
		// store extracted value as string representation!
		std::ostringstream data; data.precision(16);
		data << output;
		dict[key] = data.str();
		return true;
	}
	return false;
}
////////////////////////////////

bool conditionalCheckEqual(const std::string &line, const std::string &key);

bool conditionalCheckEqual(
		const std::string &line, const std::string &key,
		std::map<std::string, std::string>& dict);

void mkdir_p(const std::string &dirpath);


/**
 * Works on two solutions (collections of jets + value at t0) of the same kind.
 * Makes one solution into another by dumping higher order terms in jets
 *
 * ***BEWARE***: (1) it copies the value at t0, (2) the output jet controls
 * how many jets are copied. If not enough jets in v then ERROR
 * (seg. fault most probably)!
 */
template<typename SolType>
void copyReduce(const SolType &v, SolType &u) {
	u.setValueAtCurrent(v.getValueAtCurrent());
	auto jetu = u.begin();
	auto jetv = v.begin();
	for (; jetu != u.end(); ++jetu, ++jetv){
		auto coeffu = (*jetu)->begin();
		auto coeffv = (*jetv)->begin();
		for (; coeffu != (*jetu)->end(); ++coeffu, ++coeffv){
			(*coeffu) = (*coeffv);
		}
	}
}

/**
 * This is for a nice parsing of arguments from command line
 *
 * This is not so optimal ('m' * 'argc' * 'string comparison' operations)
 * Where m is number of "parse" parameters. It could be done better, but I decided to keep
 * it simple and produce decent help string automatically.
 *
 * For now it only supports params in the following form: param=value and param
 * You can do --param=value and --param if you wish
 * For now it does not support multiple forms of params, i.e. -f and --file meaning the same, itp.
 *
 * Example use:
 *
 *     ```
 * 	   capd::ddeshelper::ArgumentParser args(argc, argv);
 *	   std::string outdir = ".";
 *	   std::string prefix = "attractor";
 *	   int p = 32;
 *	   int n = 3;
 *	   std::string plot_mode = "poincare";
 *	   args.parse("outdir=", outdir);
 *	   args.parse("prefix=", prefix, "A prefix of all the filenames");
 *	   args.parse("p=", p);
 *	   args.parse("n=", n, "order of the method");
 *	   args.parse("plot=", plot_mode, {std::string("poincare"), std::string("xpx"), std::string("phasespace")}, "Kind of plot to produce and this is a very long explanation what happens in that parameter lorem ipsum sit dolor");
 *	   int itest = 5;
 *	   std::string stest;
 *	   args.parse("stest=", stest, {std::string("poincare"), std::string("xpx"), std::string("phasespace")} );
 *	   args.parse("itest=", itest, {-1,0,1} );
 *	   bool flag = args.parse(std::string("flag"));
 *	   ```
 *
 *	TODO: (not important) add support for variation in params, eg. args.parse({"-h", "--help"}), itp. should be easy
 *	TODO: (not important) extract those from the project and send to github as a small library
 *	TODO: (not important) add option to handle separated values, ie. "--param data"
 *	TODO: (not important) add option to include help header and help footer (before and after params)
 *	TODO: (not important) move implementation outside class to a .cpp/.hpp file.
 */
class ArgumentParser {
	const unsigned int MAX_LENGTH = 78;
	typedef char** ArgvType;
	typedef int ArgcType;
	ArgcType& argc;
	ArgvType& argv;
	std::ostringstream help_oss;
	std::string current_pad;
	unsigned int current_col;

	/** starts new param in the help, prepares padding for a nice output */
	void startNewParam(const std::string &param) {
		help_oss << "\n";
		std::string extra_pad = "   ";
		std::string new_pad = "";
		for (auto c : param)
			new_pad += " ";
		new_pad += extra_pad;
		while (new_pad.length() < current_pad.length()){
			extra_pad += " ";
			new_pad += " ";
		}
		// not very likely, but...
		if (new_pad.length() >= MAX_LENGTH){
			new_pad = extra_pad;
			help_oss << "\n" << new_pad;
		}
		current_pad = new_pad;
		help_oss << param << extra_pad;
		current_col = current_pad.length();
	}

	/** outputs with a nice pad */
	template<typename T>
	friend ArgumentParser& operator<<(ArgumentParser &out, const T &item) {
		std::ostringstream tmp; tmp << item;
		char stop = 1;
		std::string tmps = tmp.str() + stop;
		std::string word = "";
		for (auto& c : tmps){
			// std::cout << "c|" << c << "|word|" << word << std::endl;
			if (word.length() + out.current_pad.length() >= out.MAX_LENGTH){
				// word is too long to put in a one line, split...
				out.help_oss << word; word = "";
				out.current_col = out.current_pad.length();
			}else if (out.current_col >= out.MAX_LENGTH){
				out.help_oss << "\n" << out.current_pad;
				out.current_col = out.current_pad.length();
			} else if (c == ' ' || c == '\t') {
				out.help_oss << word; word = "";
				out.help_oss << " "; ++out.current_col;
				continue;
			} else if (c == '\n' || c == '\r') {
				out.help_oss << word; word = "";
				out.help_oss << "\n" << out.current_pad;
				out.current_col = out.current_pad.length();
				continue;
			} else if (c == stop) {
				out.help_oss << word;
				continue;
			}
			word += c; ++out.current_col;
		}
		return out;
	}
public:
	/**
	 * this should be called first after entering main(int argc, char* argv[]).
	 * Parameters are obvious from the above.
	 */
	ArgumentParser(ArgcType& argc, ArgvType& argv, int help_line_length=78):
		argc(argc), argv(argv), current_col(0), MAX_LENGTH(help_line_length) {

	}

	/**
	 * returns the command line formatted.
	 *
	 * The default sep makes it such that it is a well formatted
	 * bash command, but with comments in subsequent lines,
	 * so that if you write:
	 * cout << "# " << args.getCommandLine() << endl;
	 * you end up with a fully commented out program execution command.
	 *
	 * If you use sep = " \\\n    ", then you will have uncommented lines
	 * If you use sep = " ", then you will get command line in one line.
	 */
	std::string getCommandLine(std::string sep=" \\\n#    "){
		std::ostringstream out;
		out << argv[0];
		for (ArgcType i = 1; i < argc; ++i)
			out << sep << argv[i];
		return out.str();
	}

	/**
	 * checks with parse() for -h, --h, --help, and returns
	 * true if anyone of those found. If you use this function,
	 * please do not parse for help yourself.
	 */
	bool isHelpRequested() {
		return parse("-h") || parse("--help") || parse("--h");
	}

	/** Returns the help string with all parameters formatted */
	std::string getHelp() {
		return help_oss.str();
	}

	/**
	 * parsing a flag, i.e. a parameter without value, such as -h, --help, itp.
	 * help_s is the extra message to include in help string.
	 */
	bool parse(const std::string &param, const std::string &help_s="") {
		startNewParam(param);
		if (help_s != "") (*this) << help_s << "\n";
		(*this) << "[this is a flag parameter (present/absent)]";
		std::string dump;
		for (ArgcType i = 0; i < argc; ++i)
			if (argv[i] == param)
				return true;
		return false;
	}

	/**
	 * Basic parsing template - it takes a prefix (param), and then
	 * tries to parse with istringstream the type T to out.
	 * The value given in out will be used as default value shown in help.
	 * If parsing fails, the value in out will stay not changed.
	 * help_s is the extra information to ut into help string for this paramater.
	 *
	 * It returns true if the requested parameter was in the command line.
	 */
	template<typename T>
	bool parse(const std::string &param, T &out,
			const std::string &help_s = "") {
		startNewParam(param);
		if (help_s != "") (*this) << help_s << "\n";
		(*this) << "[default: " << out << "]";
		bool found = false;
		for (ArgcType i = 0; i < argc; ++i)
			found = found or conditionalExtractValue(argv[i], param, out);
		return found;
	}

	/** same as the other tempated parse, but accepts old-fashion strings */
	template<typename T>
	bool parse(const std::string &param, T &out, const char *const help_s) {
		return parse(param, out, std::string(help_s));
	}

	/** this verison is almost the same as basic paring, but has a list of allowed values */
	template<typename T>
	bool parse(const std::string &param, T &out,
			std::initializer_list<T> allowed_values,
			const std::string &help_s="") {
		startNewParam(param);
		if (help_s != "") (*this) << help_s << "\n";
		(*this) << "[default: " << out << "] ";
		(*this) << "[allowed: "; std::string sep = " ";
		for (auto v : allowed_values){
			(*this) << sep << v;
			sep = ", ";
		}
		(*this) << "] ";
		T tmp_out;
		bool found = false;
		for (ArgcType i = 0; i < argc; ++i){
			if (!conditionalExtractValue(argv[i], param, tmp_out)) continue;
			bool allowed = false;
			for (auto v : allowed_values){
				if (v == tmp_out){
					allowed = true;
					break;
				}
			}
			if (allowed){
				out = tmp_out;
				found = true;
			}
		}
		return found;
	}

	/**
	 * this can check parsed value against given predicate (e.g. lambda [](T item) -> bool)
	 * since I cannot get docs automatically from predicate, you should provide
	 * information for user in help_s explaining what are the accepted values.
	 */
	template<typename T, typename F>
	bool parse(const std::string &param, T &out,
			F predicate, const std::string &help_s="") {
		startNewParam(param);
		if (help_s != "") (*this) << help_s << "\n";
		(*this) << "[default: " << out << "] ";
		(*this) << "[note: has predicate, check docs] ";
		T tmp_out;
		bool found = false;
		for (ArgcType i = 0; i < argc; ++i){
			if (!conditionalExtractValue(argv[i], param, tmp_out)) continue;
			if (F(tmp_out)){
				out = tmp_out;
				found = true;
			}
		}
		return found;
	}
};

} // namespace ddeshelper
} // namespace capd

#endif /* _CAPD_DDES_DDESHELPERSLIB_H_ */
