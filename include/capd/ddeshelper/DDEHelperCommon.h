/*
 * DDEHelperCommon.h
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

#ifndef _CAPD_DDEHELPERCOMMON_H_
#define _CAPD_DDEHELPERCOMMON_H_

#include <capd/ddes/ddeslib.h>

#include <fstream>

namespace capd {
namespace ddeshelper {


template<typename T>
std::string slice(T const& item, int len = 72){
	std::ostringstream oss;
	oss << item;
	return oss.str().substr(0, len);
}

/** wrapper for calling system commands, we use it for plotting */
void runSystemCommand(std::string cmd);

/**
 * used in some helpers to define where the output goes.
 * It does not check the correctness of those paths! You need to assure they are correct!
 *
 * dirPath - a base directory to where the output save
 * prefix - as name says - prefix of all generated files.
 */
class PathConfig {
public:
	std::string dirPath;
	std::string prefix;

	/** make a config path, default is ./ with no prefix */
	PathConfig(std::string dirPath=".", std::string prefix=""):
			dirPath(dirPath), prefix(prefix) {
		if (dirPath != "" && dirPath.back() != '/')
			dirPath += "/";
	}

	/** make a new PathConfix with a longer prefix : "prefix-suffix" */
	PathConfig suffix(std::string suffix) const {
		return PathConfig(dirPath, prefix + "-" + suffix);
	}
	/** dirpath with prefix */
	std::string fullpath() const { return dirPath + "/" + prefix; }
	/** makes a new prefixed filepath for a given filename */
	std::string filepath(std::string name) const { return dirPath + "/" + prefix + "-" + name; }
	/** makes just a prefixed filename (no path attached) */
	std::string filename(std::string name) const { return prefix + "-" + name; }

	/** assure there is really this folder - not very safe!!! */
	void mkdir_p() const {
		if (dirPath != "/" && dirPath != "./" && dirPath != ""){
			if (dirPath != ".")
				capd::ddeshelper::runSystemCommand(std::string("mkdir -p '") + dirPath + "' ");
		}else{
			throw std::logic_error("PathConfig::mkdir_p(): path is unsafe!");
		}
	}

	std::string dirpath() const  { return dirPath; };
};

// TODO: (NOT URGENT) move those to other files!
template<typename T>
union BinaryData {
	T value;
	char repr[sizeof(T)];

	BinaryData(){};
	BinaryData(T v): value(v) {};

	~BinaryData(){}
};
template<typename T>
std::streamsize binsize(const BinaryData<T>& d){
	return sizeof(T);
}

// TODO: write helpers for readText ?

/**
 * User must supply out which supports write and is set to binary!
 * It does not store Vector size, it might be used for both vectors and matrices
 *
 * TODO: (NOT URGENT) move those to other files!
 */
template<typename VectorType>
void saveBinary(std::ostream& out, const VectorType& vector){
	typedef typename VectorType::ScalarType ScalarType;
	for (auto it = vector.begin(); it != vector.end(); ++it){
		BinaryData<ScalarType> value(*it);
		out.write(value.repr, binsize(value));
	}
}
/**
 * User must supply out which supports write and is set to binary!
 * TODO: (NOT URGENT) move those to other files!
 */
template<typename ScalarType>
void saveBinaryScalar(std::ostream& out, const ScalarType& v){
	BinaryData<ScalarType> value(v);
	out.write(value.repr, binsize(value));
}
/**
 * User must supply out which supports write and is set to binary!
 * TODO: (NOT URGENT) move those to other files!
 * TODO: (NOT URGENT) there is a typo in the name, it should be Builtin or BuiltIn...
 */
template<typename BuildinType>
void saveBinaryBuildin(std::ostream& out, const BuildinType& v){
	BinaryData<BuildinType> value(v);
	out.write(value.repr, binsize(value));
}
/**
 * For simplicity...
 *
 * TODO: (NOT URGENT) move those to other files!
 */
template<typename VectorType>
void saveBinary(const std::string& filepath, const VectorType& vector){
	std::ofstream out(filepath, std::ios::binary);
	saveBinary(out, vector);
	out.close();
}

/**
 * User must supply in which supports read and is set to binary!
 * Vector must be of appropriate size
 * It works for both matrix and vector! (any Container!)
 *
 * TODO: (NOT URGENT) move those to other files!
 */
template<typename VectorType>
void readBinary(std::istream& in, VectorType& vector){
	typedef typename VectorType::ScalarType ScalarType;
	for (auto it = vector.begin(); it != vector.end(); ++it){
		BinaryData<ScalarType> value;
		in.read(value.repr, binsize(value));
		*it = value.value;
	}
}
/**
 * User must supply in which supports read and is set to binary!
 * TODO: (NOT URGENT) move those to other files!
 */
template<typename ScalarType>
void readBinaryScalar(std::istream& in, ScalarType& v){
	BinaryData<ScalarType> value;
	in.read(value.repr, binsize(value));
	v = value.value;
}
/**
 * User must supply in which supports read and is set to binary.
 * TODO: (NOT URGENT) move those to other files!
 */
template<typename BuildinType>
void readBinaryBuiltin(std::istream& in, BuildinType& v){
	BinaryData<BuildinType> value;
	in.read(value.repr, binsize(value));
	v = value.value;
}
/**
 * For simplicity...
 *
 * TODO: (NOT URGENT) move those to other files!
 */
template<typename VectorType>
void readBinary(const std::string& filepath, VectorType& vector){
	std::ifstream in(filepath, std::ios::binary);
	readBinary(in, vector);
	in.close();
}

/**
 * used to conventionally read data from files.
 * If a vector/matrix is given in line, then it will be read
 * Otherwise, if # or @ is the first character (comment), it will be discarded (and false returned)
 * Otherwise, it will try to read a vector/matrix from a file given in the line.
 * Returns true if item was successfully read
 *
 * TODO: (NOT URGENT) move those to other files!
 */
template<typename T>
bool relativeReadData(std::istream& in, T& item){
	std::string line; in >> line; // TODO: (FUTURE) use readline for better handling
	// check if line is explicitly given vector first
	if (line[0] == '{'){ std::istringstream iss(line); iss >> item; return true; }
	// check if it is a comment
	else if (line[0] == '#' || line[0] == '@') return false;
	else {
		std::ifstream for_text_input(line);
		if (for_text_input.peek() == '{'){ // this should be enough for our purpose
			for_text_input >> item;
			for_text_input.close();
		}else{
			for_text_input.close();
			// we read vector assuming it is a binary representation.
			readBinary(line, item);
		}
		return true;
	}
}

/**
 * this is a helper function that returns the combined dimension of several explicitely given vectors.
 * Usage:
 * 		auto d = commonDimension(v, u, w);			// any number of arguments.
 * 		auto d = commonDimension(v, u, w, x, y, z); // etc.
 * TODO: test it.
 */
template<typename VectorType, typename ...VectorTypes>
typename VectorType::size_type
commonDimension(VectorType const & first, VectorTypes const & ...others) {
    return first.dimension() + (sizeof...(others) ? commonDimension(others...) : 0);
}
/**
 * this is a helper function that puts all the vectors into out vectors,
 * sequentially: out = (first[0],..., first[first.dimension()-1], ...)
 * WARNING this does not resize out to accommodate all and only the elements of input vectors.
 * WARNING: if you want to get the right dimension, use the joinVectors(...) function.
 * Usage:
 *      capd::DVector v(1), u(2), w(3), out(8);
 * 		joinVectorsRaw(out, v, u, w); // out will still be dimension 8, and first 6 coordinates will be v, u, w.
 * TODO: test it.
 */
template<typename OutIterator, typename VectorType, typename ...VectorTypes>
void joinVectorsRaw(OutIterator& out, VectorType const & first, VectorTypes const & ...others) {
	for (auto vIt = first.begin(); vIt != first.end(); ++vIt, ++out)
		*out = *vIt;
	if (sizeof...(others))
		joinVectorsRaw(out, others...);
}
/**
 * this is a helper function that puts all the vectors into one vector,
 * sequentially: out = (first[0],..., first[first.dimension()-1], ...)
 * Usage:
 *      capd::DVector v(1), u(2), w(3);
 * 		auto out = joinVectorsRaw(out, v, u, w); // out will be of dimension 6
 * TODO: test it.
 */
template<typename VectorType, typename ...VectorTypes>
VectorType joinVectors(VectorType const & first, VectorTypes const & ...others) {
	typename VectorType::size_type dim = commonDimension(first, others...);
	VectorType result(dim);
	joinVectorsRaw(result.begin(), first, others...);
	return result;
}

} // namespace ddeshelper
} // namespace capd

#endif /* _CAPD_DDEHELPERCOMMON_H_ */
