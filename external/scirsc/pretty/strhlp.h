#ifndef _CAPD_CUSTOM_STRHLP_H_
#define _CAPD_CUSTOM_STRHLP_H_

#include <iostream>
#include <sstream>
#include <math.h>

// TODO: add pretty namespaces

namespace scirsc{

/**
 * some usefull functions to manipulate strings and converting
 * numbers into strings with padding and such. 
 */
namespace strhlp{

/**
 * simple as that...
 */
bool startsWith(std::string item, std::string prefix){
	if (prefix.length() > item.length()){
		return false;
	}
	for (int i = 0; i < prefix.length(); i++){
		if (item[i] != prefix[i]){
			return false;
		}
	}
	return true;
}


/**
 * Prints formatted string:
 * Examples:
 *
 *  a) rjust(123) 			-> 	"123"
 *  b) rjust(123, 6)  		-> 	"   123"
 *  c) rjust(123, 6, '0')	-> 	"000123"
 *  d) rjust(12345678, 6)	-> 	"12345678"
 *
 * Works on double and int.
 * Do not align dots in doubles.
 *
 * This is not wery efficient...
 */
template<typename T>
std::string rjust(T data, int digits = -1, char pad = ' '){
	std::ostringstream os; if (digits > 0) os.precision(digits-5);
	os << data;
	std::string text = os.str();
	while (digits > 0 && text.length() < digits ) text = pad + text;
	return text;
}


/**
 * Prints formatted string:
 * Examples:
 *
 *  a) ljust(123) 		    -> 	"123"
 *  b) ljust(123, 6)  	    -> 	"123   "
 *  c) ljust(123, 6, '0')	-> 	"123000"
 *  d) ljust(12345678, 6)	-> 	"12345678"
 *
 * Works on double and int.
 * Do not align dots in doubles.
 *
 * This is not wery efficient...
 */
template<typename T>
std::string ljust(T data, int digits = -1, char pad = ' '){
	std::ostringstream os; if (digits > 0) os.precision(digits-5);
	os << data;
	std::string text = os.str();
	while (digits > 0 && text.length() < digits ) text = text + pad;
	return text;
}


/**
 * Backward compatibility
 */
template<typename T>
std::string str(T data, int digits = -1, char pad = ' '){
	return rjust(data, digits, pad);
}

/**
 * fixed / forced length
 */
template<typename T>
std::string str_fl(T data, int digits = -1, char pad = ' '){
	std::ostringstream os; if (digits > 0) os.precision(digits-3);
	os << data;
	std::string text = os.str();
	int tries = 1;
	while (digits > 0 && text.length() > digits && (digits - 3 - tries > 0)){
		std::ostringstream os; os.precision(digits - 3 - tries);
		os << data;
		text = os.str();
		++tries;
	}
	while (digits > 0 && text.length() < digits ) text = pad + text;
	return text;
}


/**
 * Fill with zeros
 */
template<typename T>
std::string zfill(T data, int digits = -1){
	return rjust(data, digits, '0');
}


/**
 *
 */
template<typename T>
std::string fixed(T data, int dl, int dr){
	std::ostringstream os; os.precision(dl + dr + 1);
	os << data;
	std::string l = os.str(), r;
	int i = l.find(".");
	r = (i == std::string::npos ? "0" : l.substr(i+1));
	l = l.substr(0, i);

	while (r.length() < dr) r += " ";
	while (l.length() < dl) l = " " + l;

	return (l + "." + r).substr(0, dl+dr+1);
}

/**
 * Pretty Print eigen vectors and eigenvalues
 */
template<typename OutputSpec, typename VectorSpec, typename MatrixSpec>
OutputSpec& prettyPrintEigenvaluesAndEigenvectors(
				OutputSpec& output,
				const VectorSpec& real_val,
				const VectorSpec& imag_val,
				const MatrixSpec& real_vec,
				const MatrixSpec& imag_vec){
	int dim = real_val.dimension();
	typedef typename VectorSpec::ScalarType ScalarType;
	VectorSpec modulus(dim);
	for (int i = 0; i < dim; i++){
		// modulus[i] = sqrt(power(real_val[i], 2) + power(imag_val[i], 2));   // clang complains here about power... ??? I rephrased it below, but should find the reason in the future...
		modulus[i] = sqrt(real_val[i] * real_val[i] + imag_val[i] * imag_val[i]);
		std::string s_modulus = strhlp::str(modulus[i], 19);
		std::string s_real = strhlp::str(real_val[i], 19);
		std::string s_imag = (imag_val[i] != 0. ? strhlp::str(imag_val[i], 19) + "i" : "                    ");
		output << (imag_val[i] != 0. ? "i" : "r") << " " << s_real << (imag_val[i] != 0. ? " + " : "   ") << s_imag << " [" << s_modulus << "] : " << real_vec.column(i) << " + " << imag_vec.column(i) << "i\n";
	}
	output << "Eigenvalues stats:\nreal part: ";
	for (int i = 0; i < dim; i++){
		std::string s_real = strhlp::str(real_val[i], 19);
		output << s_real << "|";
	} output << "\n";
	output << "imag part: ";
	for (int i = 0; i < dim; i++){
		std::string s_imag = (imag_val[i] != 0. ? strhlp::str(imag_val[i], 19) : "                   ");
		output << s_imag << "|";
	} output << "\n";
	output << "modulus  : ";
	for (int i = 0; i < dim; i++){
		std::string s_modulus = strhlp::str(modulus[i], 19);
		output << s_modulus << "|";
	} output << "\n";
	return output;
}

/**
 * Reurns full comand line as string - together with program invocation
 * useful when want to show what command was called
 * WARNING: must pass argc and argv arguments - it assumes that
 * argv has at least one element (argv[0] == program name)
 */
std::string commandLine(int argc, char* argv[]){
	std::ostringstream out;
	out << argv[0]; //
	for (int i = 1; i < argc; i++){
		out << " \"" << argv[i] << "\"";
	}
	return out.str();
}

/**
 * this function copies input to output in character (text) mode.
 * return output for further operations.
 */
std::ostream& copy(std::istream& in, std::ostream& out){
	int MAX_READ = 1024;
	char ch; char data[MAX_READ+1];
	while (in.get(ch)){
		out << ch; 					// output character (usually \n)
		in.get(data, MAX_READ);		// read that much characters or until delimiter (\n) -> delimiter stays in input!
		out << data;				// output data
	}
	return out;
}

/**
 * Reads lines from the input skipping those which start with '#' (comment)
 * Then return last line without hash or empty string (all lines started with #)
 */
std::string readlineWithHashskip(std::istream& in){
	std::string str;
	while (std::getline(in, str) && str[0] == '#'){} // skip comments
	if (str.length() > 0 && str[0] == '#') str = ""; // if last line also comment - return empty string
	return str;
}


/**
 * skip empty lines / characters - it removes characters ' ', '\t', '\n', '\r'
 * from current input.
 */
std::istream& emptyskip(std::istream& in){
	int c = in.peek();
	while (c == ' ' || c == '\t' || c == '\n' || c == '\r'){
		in.get(); c = in.peek();
	}
	return in;
}


/**
 * reads comment lines in input stream (those starts with #)
 *
 * hash => character to skip (default: #)
 * howMany == -1 => skip all hashes
 * howMany == n > 0 => skip all hashes
 */
std::istream& hashskip(std::istream& in, char hash = '#', int howMany = -1){
	emptyskip(in);
	int c = in.peek();
	std::string dump;
	int i = 0;
	while (c == hash && (i < howMany || howMany < 0)){
		i++;
		std::getline(in, dump);
		c = in.peek();
	}
	return in;
}


// TODO: why not working? (error: use of deleted function ‘std::basic_stringstream<char>::basic_stringstream(const std::basic_stringstream<char>&)’)
///**
// * Reads lines from the input skipping those which start with '#' (comment)
// * Then return last line as an std::istream without hash or empty string (all lines started with #)
// */
//std::istringstream hashskip(std::istream& in){
//	std::string str;
//	while (std::getline(in, str) && str[0] == '#'){} // skip comments
//	if (str.length() > 0 && str[0] == '#') str = ""; // if last line also comment - return empty string
//	std::istringstream is(str);
//	return is;
//}

/**
 * Pretty Prints Matrix
 */
template<typename MatrixSpec>
std::ostream& prettyPrintM(std::ostream& out, const MatrixSpec& C, std::string prefix = ""){
	for (int i = 1; i <= C.numberOfRows(); i++){
		out << prefix;
		for (int j = 1; j <= C.numberOfColumns(); j++){
			out << strhlp::fixed(C(i,j), 3, 8) << ", " << std::flush;
		}
		out << "\n" << std::flush;
	}
	return out;
}

/**
 * Pretty Prints Vector
 */
template<typename VectorSpec>
std::ostream& prettyPrintV(std::ostream& out, const VectorSpec& v, std::string prefix = ""){
	out << prefix;
	for (int i = 0; i < v.dimension(); i++){
		out << strhlp::str(v[i], 15) << ", "<< std::flush;
	}
	return out;
}

// for backward compatibility
template<typename MatrixSpec>
std::ostream& prettyPrint(std::ostream& out, const MatrixSpec& C, std::string prefix = ""){
	return prettyPrintM(out, C, prefix);
}


/**
 * returns latex ready rows for pretty formating fourier functions
 */
template<typename AScalarSpec, typename BScalarSpec>
std::string texFourierPair(AScalarSpec c, BScalarSpec s, int k, std::string T = "{2\\pi \\over T}", std::string esc = "\\"){
	std::ostringstream ss, cs;
	ss << s; cs << c;
	std::string sss = ss.str(), css = cs.str();
	char ssign = sss[0]; if (ssign != '-') ssign = '+';
	char csign = css[0]; if (csign != '-') csign = '+';
	if (sss[0] == '-' || sss[0] == '+') sss = sss.substr(1);
	if (css[0] == '-' || css[0] == '+') css = css.substr(1);
	std::ostringstream o;
	o << csign << " & " << css << " " << esc << "cdot " << esc << "cos" << esc << "left( " << T << " " << esc << "cdot " << k << " " << esc << "cdot t " << esc << "right) & " << ssign << " & " << sss << " " << esc << "cdot " << esc << "sin" << esc << "left( " << T << " " << esc << "cdot " << k << " " << esc << "cdot t " << esc << "right)";
	return o.str();
}


} // namespace strhlp

} // namespace scirsc

#endif /* _CAPD_CUSTOM_STRHLP_H_ */
