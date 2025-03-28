/*
 * DDEHelperDrawing.cpp
 *
 *  Created on: Feb 12, 2024
 *      Author: robson
 */

#include <capd/ddes/DDECommon.h>
#include <capd/ddeshelper/DDEHelperDrawing.h>

namespace capd {
namespace ddeshelper {

std::pair<std::string, std::string> split_pathprefix(std::string const& s){
	auto pos = s.rfind("/");
	if (pos != std::string::npos){
		return std::make_pair(s.substr(0, pos+1), s.substr(pos+1));
	}else{
		return std::make_pair(std::string(""), s);
	}
}

std::pair<std::string, std::string> split_name_ext(std::string const& s){
	auto pos = s.rfind(".");
	if (pos != std::string::npos){
		return std::make_pair(s.substr(0, pos), s.substr(pos+1));
	}else{
		return std::make_pair(s, std::string(""));
	}
}

void to_dat(std::ostream& out, capd::interval const & v){
	capd::interval u = v, r;
	u.split(r);
	out << u.leftBound() << "  " << r.rightBound();
}

void to_dat(std::ostream& out, double const & v){
	out << v << " 0.0";
}

void plot_datfile(
		std::string const& pathprefix,
		std::vector<std::string> const& gpoptions,
		bool live,
		std::string filesuffix){
	auto dirandname = split_pathprefix(pathprefix);
	auto nameandext = split_name_ext(dirandname.second);
	auto dirpath = dirandname.first;
	auto datname = dirandname.second;
	if (nameandext.second == "") datname += ".dat";
	auto fileprefix = nameandext.first;
	#ifdef DDES_ALLOW_SYSTEM
	if (dirpath != "") capd::ddeshelper::runSystemCommand(std::string("mkdir -p " + dirpath));
	#endif
	std::string datpath = dirpath + datname;
	std::string pngname = fileprefix + filesuffix + ".png";
	std::string gpname = fileprefix + filesuffix + ".gp";
	std::string gppath = dirpath + gpname;
	std::ofstream outg(gppath);
	if (!live){
		outg << "set terminal png size 1600,1200" << std::endl;
		outg << "set output '" << pngname << "'" << std::endl;
	}
	outg << "plot \\\n";
	for (auto& item: gpoptions){
		outg << "    '" << datname << "' " << item << ", \\\n";
	}
	outg << "\n\n";
	outg.close();
	#ifdef DDES_ALLOW_SYSTEM
	std::ostringstream cmd; cmd << "cd '" << dirpath << "' && gnuplot " << (live ? "-p " : "") << "'" << gpname << "'";
	capd::ddeshelper::runSystemCommand(cmd.str());
	#endif
}

void xplot_many(
		std::string const& pathprefix,
		std::string plottype,
		std::vector<std::string> const& gpplots,
		bool live,
		std::string filesuffix){
	auto dirandname = split_pathprefix(pathprefix);
	auto dirpath = dirandname.first;
	auto fileprefix = dirandname.second;
	#ifdef DDES_ALLOW_SYSTEM
	if (dirpath != "") capd::ddeshelper::runSystemCommand(std::string("mkdir -p " + dirpath));
	#endif
	std::string pngname = fileprefix + filesuffix + ".png";
	std::string gpname = fileprefix + filesuffix + ".gp";
	std::string gppath = dirpath + gpname;
	std::ofstream outg(gppath);
	if (!live){
		outg << "set terminal png size 1600,1200" << std::endl;
		outg << "set output '" << pngname << "'" << std::endl;
	}
	if (gpplots.size() > 0){
		std::string sep = " ";
		outg << plottype;
		for (auto& item: gpplots){
			outg << sep << item;
			sep = ", \\\n";
		}
		outg << "\n\n";
	}
	outg.close();
	#ifdef DDES_ALLOW_SYSTEM
	std::ostringstream cmd; cmd << "cd '" << dirpath << "' && gnuplot " << (live ? "-p " : "") << "'" << gpname << "'";
	capd::ddeshelper::runSystemCommand(cmd.str());
	#endif
}

void plot_many(
		std::string const& pathprefix,
		std::vector<std::string> const& gpplots,
		bool live,
		std::string filesuffix){
	xplot_many(pathprefix, "plot", gpplots, live, filesuffix);
}

void splot_many(
		std::string const& pathprefix,
		std::vector<std::string> const& gpplots,
		bool live,
		std::string filesuffix){
	xplot_many(pathprefix, "splot", gpplots, live, filesuffix);
}

} // namespace ddeshelper
} // namespace capd
