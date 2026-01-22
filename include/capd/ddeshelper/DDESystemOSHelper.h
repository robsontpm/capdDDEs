/*
 * DDESystemOSHelper.h
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

#ifndef _CAPD_DDESYSTEMOSHELPER_H_
#define _CAPD_DDESYSTEMOSHELPER_H_

#include <ctime>
#include <iomanip>
#include <capd/ddeshelper/DDEHelperCommon.h>

namespace capd {
namespace ddeshelper {

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

std::string makeCommandLine(int argc, char** argv);

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

// TODO: (NOT URGENT): this is because the files are all messed up. I need to cleanup those helpers.
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
 *
 *	NOTE: a lot of the methods does not have const modifier. I assume this object is rather volatile and it has only one instance.
 */
class ArgumentParser {
	const unsigned int MAX_LENGTH = 78;
	typedef char** ArgvType;
	typedef int ArgcType;
	typedef std::string TokenType;
	typedef std::vector<TokenType> TokenListType;
	ArgcType argc;
	ArgvType argv;
	TokenListType tokens;
	std::ostringstream help_oss;
	std::ostringstream parsing_log;
	std::string current_pad;
	unsigned int current_col;
	std::string parse_time = "";

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
	/**
	 * converts argv to tokens
	 * TODO: NOTE: this might be used now to convert consecutive pairs into tokens, like in the old programs. Rethink!
	 */
	void argv2Tokens(){
		parsing_log << "argv/argc:" << argc << "\n";
		for (int i = 0; i < argc; ++i) {
			parsing_log << "token: " << argv[i] << "\n";
			tokens.push_back(std::string(argv[i]));
		}
	}
public:
	/**
	 * this should be called first after entering main(int argc, char* argv[]).
	 * Parameters are obvious from the above.
	 */
	ArgumentParser(ArgcType argc, ArgvType argv, int help_line_length=78):
		argc(argc), argv(argv), current_col(0), MAX_LENGTH(help_line_length) {
		std::time_t time = std::time(0);
		const int M = 30; char time_c_str[M+1] = {};
		std::strftime(time_c_str, M, "%c", std::localtime(&time));
		parse_time = std::string(time_c_str);
		argv2Tokens();
	}
	/** repopulate data from the given argc and argv. Probably not used that much. */
	void setup(ArgcType argc, ArgvType argv, bool remove_old=false) {
		this->argv = argv;
		this->argc = argc;
		if (remove_old) tokens.clear();
		argv2Tokens();
	}
	/**
	 * setup data from the given stream. The stream is read in whole,
	 * lines that starts with # are skipped, and each line becomes one token.
	 */
	void setup(std::istream& input, bool remove_old=false) {
		if (remove_old) tokens.clear();
	    std::string line;
	    while (std::getline(input, line)) {
	    	parsing_log << "line: " << line << "\n";
	    	// erase the front white spaces
	        line.erase(0, line.find_first_not_of(" \t"));
	        // if the line is just a comment or empty, do nothing
	        if (line.empty() || line[0] == '#') continue;
	        // now we have, that line has some chars not being white and for sure don't start with #
	        std::istringstream iss(line);
	        std::getline(iss, line, '#'); // this will remove potential comment at the end of the line
	        // we have removed comment, but we might have some trailing white spaces.
	        // but we know that there must be some non-white chars in the front, so we
	        // do not need to worry about find_last_not_of returning npos.
	        line.erase(line.find_last_not_of(" \t") + 1);
	        // line is in good shape now!
	        parsing_log << "token: " << line << "\n";
	        tokens.push_back(line);
	    }
	}
	/**
	 * setup data from the given filepath (read through ifstream).
	 * @see setup(std::ifstream&, bool)
	 */
	void setupFromFile(std::string& filepath, bool remove_old=false) {
		parsing_log << "file: " << filepath << "\n";
	    std::ifstream f(filepath);
	    setup(f, remove_old);
	    f.close();
	}
	/**
	 * Checks if the given parameter is present in the current tikens list (e.g. from argv)
	 * if so, it loads the file given by this param, and then repopulates the
	 * token list as given in @see setup(std::ifstream&, bool).
	 * If overwrite_argv=false (default), then argv params takes precedence over those from files
	 * if it is true, then the config file params takes precedence.
	 */
	void checkConfigFile(
				std::string param="--config=",
				std::string help_s="The configuration file for extra parameters.",
				bool overwrite_argv=false){
		std::string filepath;
		if (this->parse(param, filepath, help_s)){
			this->setupFromFile(filepath);
			// repopulate the argv in the end of the list,
			// so that they takes precendence
			if (!overwrite_argv) argv2Tokens();
		}
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
	/** returns the representation of the time at what the argparser was created (~start of the program) */
	std::string getStartTime(){
		return parse_time;
	}
	/** returns a nice description of the parse, to be put in the output of files for future reference. */
	std::string getStoryMessage(std::string sep="\\\n#    "){
		return "This program was run on " + parse_time + ", with the following command: " + sep + getCommandLine(sep);
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
	std::string getHelp() const {
		return help_oss.str();
	}

	/** Returns the log of all parsed tokens. Might be useful for debug or for history. */
	std::string getParsingLog() const {
		return parsing_log.str();
	}

	/**
	 * parsing a flag, i.e. a parameter without value, such as -h, --help, itp.
	 * help_s is the extra message to include in help string.
	 */
	bool parse(const std::string &param, const char *const help_s=0) {
		startNewParam(param);
		if (help_s) (*this) << help_s << "\n";
		(*this) << "[this is a flag parameter (present/absent)]";
		std::string dump;
//		for (ArgcType i = 0; i < argc; ++i)
//			if (argv[i] == param)
//				return true;
		for (auto& item: tokens)
			if (item == param)
				return true;
		return false;
	}

	/** same as the other templated parse, but accepts strings */
	template<typename T>
	bool parse(const std::string &param, const std::string& help_s) {
		return parse(param, help_s.c_str());
	}

	/**
	 * Basic parsing template - it takes a prefix (param), and then
	 * tries to parse with istringstream the type T to out.
	 * The value given in out will be used as default value shown in help.
	 * If parsing fails, the value in out will stay not changed.
	 * help_s is the extra information to put into help string for this paramater.
	 *
	 * It returns true if the requested parameter was in the command line.
	 */
	template<typename T>
	bool parse(const std::string &param, T &out,
			const char *const help_s=0) {
		startNewParam(param);
		if (help_s) (*this) << help_s << "\n";
		(*this) << "[default: " << std::setprecision(16) << out << "]";
		bool found = false;
//		for (ArgcType i = 0; i < argc; ++i)
//			found = conditionalExtractValue(argv[i], param, out) or found;
		for (auto& item: tokens)
			found = conditionalExtractValue(item, param, out) or found;
		return found;
	}

	/** same as the other templated parse, but accepts strings */
	template<typename T>
	bool parse(const std::string &param, T &out, const std::string &help_s) {
		return parse(param, out, help_s.c_str());
	}

	/** this verison is almost the same as basic paring, but has a list of allowed values */
	template<typename T>
	bool parse(const std::string &param, T &out,
			std::initializer_list<T> allowed_values,
			const char *const help_s=0) {
		startNewParam(param);
		if (help_s) (*this) << help_s << "\n";
		(*this) << "[default: " << out << "] ";
		(*this) << "[allowed: "; std::string sep = " ";
		for (auto v : allowed_values){
			(*this) << sep << v;
			sep = ", ";
		}
		(*this) << "] ";
		T tmp_out;
		bool found = false;
//		for (ArgcType i = 0; i < argc; ++i){
//			if (!conditionalExtractValue(argv[i], param, tmp_out)) continue;
		for (auto& item: tokens){
			if (!conditionalExtractValue(item, param, tmp_out)) continue;
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

	/** to accept strings as help */
	template<typename T>
	bool parse(const std::string &param, T &out,
			std::initializer_list<T> allowed_values,
			const std::string &help_s="") {
		return parse(param, out, allowed_values, help_s.c_str());
	}

	/**
	 * this can check parsed value against given predicate (e.g. lambda [](T item) -> bool)
	 * since I cannot get docs automatically from predicate, you should provide
	 * information for user in help_s explaining what are the accepted values.
	 */
	template<typename T, typename F>
	bool parse(const std::string &param, T &out,
			F predicate, const char *const help_s=0) {
		startNewParam(param);
		if (help_s) (*this) << help_s << "\n";
		(*this) << "[default: " << out << "] ";
		(*this) << "[note: has predicate, check docs] ";
		T tmp_out;
		bool found = false;
//		for (ArgcType i = 0; i < argc; ++i){
//			if (!conditionalExtractValue(argv[i], param, tmp_out)) continue;
		for (auto& item: tokens){
			if (!conditionalExtractValue(item, param, tmp_out)) continue;
			if (F(tmp_out)){
				out = tmp_out;
				found = true;
			}
		}
		return found;
	}

	/** to accept strings as help */
	template<typename T, typename F>
	bool parse(const std::string &param, T &out,
			F predicate, const std::string &help_s="") {
		return parse(param, out, predicate, help_s.c_str());
	}
};

/**
 * used in some helpers to define where the output goes.
 * WARNING: It does not check the correctness of those paths! You need to assure they are correct!
 * WARNING: It is not malicious attack safe! Use only when you trust the user and the input!
 *
 * dirPath - a base directory to where the output save
 * prefix - as name says - prefix of all generated files.
 */
class PathConfig {
public:
	std::string dirPath;
	std::string prefix;

	/** make a config path, default is ./ with no prefix */
	PathConfig(std::string dirPath="./", std::string prefix="", bool absolute_allowed=false):
			dirPath(dirPath), prefix(prefix) {
		sanitize(absolute_allowed);
	}
	/**
	 * Makes a config reading parameters wd={dirPath} and prefix={prefix} from command line parser.
	 * You can change the parameter names for parser in the optional arguments.
	 */
	PathConfig(ArgumentParser& args,
			bool absolute_allowed=false,
			std::string wd_param="wd=",
			std::string prefix_param="prefix=",
			std::string wd_default="./",
			std::string prefix_default="",
			std::string wd_help="A path to the working directory.",
			std::string prefix_help="A prefix to all of the generated files."
	): dirPath(wd_default), prefix(prefix_default){
		args.parse(wd_param, dirPath, wd_help);
		args.parse(prefix_param, prefix, prefix_help);
		sanitize(absolute_allowed);
	}

	/** make a new PathConfix with a longer prefix : "prefix-suffix" */
	PathConfig suffix(std::string suffix) const {
		return PathConfig(dirPath, (prefix == "" ? "" : "-") + suffix);
	}
	/** dirpath with prefix */
	std::string fullpath() const { return dirPath + "/" + prefix; }
	/** makes a new prefixed filepath for a given filename */
	std::string filepath(std::string name) const { return fullpath() + (prefix == "" ? "" : "-") + name; }
	/** makes just a prefixed filename (no path attached) */
	std::string filename(std::string name) const { return prefix + (prefix == "" ? "" : "-") + name; }

	/** assure there is really this folder - not very safe!!! */
	bool mkdir_p() const {
		if (dirPath != "/" && dirPath != "./" && dirPath != "" && dirPath != "."){
			capd::ddeshelper::runSystemCommand(std::string("mkdir -p '") + dirPath + "' ");
			return true;
		}else{
			return false;
		}
	}

	std::string dirpath() const  { return dirPath; };

private:
	/** WARNING: this is by no means industry standard, just some randommost common checks!! */
	void sanitize(bool absolute_allowed){
		if (dirPath != "" && dirPath.back() != '/')
			dirPath += "/";
		if (!absolute_allowed){
			if (dirPath.length() > 0 && dirPath[0] == '/')
				throw std::logic_error("PathConfig::__construct__(): absolute paths are not allowed: '" + dirPath + "'");
		}
	}
};

} // namespace ddeshelper
} // namespace capd

#endif /* _CAPD_DDESYSTEMOSHELPER_H_ */
