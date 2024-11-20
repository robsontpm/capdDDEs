/*
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

#ifndef _CAPD_DDEHELPERDRAWING_H_
#define _CAPD_DDEHELPERDRAWING_H_

#include <capd/ddes/DDECommon.h>
#include <capd/ddeshelper/DDEHelperCommon.h>

#include <capd/ddes/DDEForwardTaylorCurvePiece.hpp>
#include <string>

namespace capd {
namespace ddeshelper {

/** splits "path/to/a/file.txt" into "path/to/a/" and "file.txt" */
std::pair<std::string, std::string> split_pathprefix(std::string const& s);

/** splits "file.txt" into "file" and "txt" */
std::pair<std::string, std::string> split_name_ext(std::string const& s);

/** saves as several double values separated by spaces. Special case for capd::intervals::Interval template */
template<typename T_Bound, typename T_Rnd>
void to_dat(std::ostream& out, capd::intervals::Interval<T_Bound, T_Rnd> const & v){
	typedef capd::intervals::Interval<T_Bound, T_Rnd> T_Intv;
	T_Intv u = v, r;
	u.split(r);
	out << u.leftBound() << " " << r.rightBound();
}

/** saves as several double values separated by spaces. special case for capd::interval (TODO: (NOT URGENT) add handling for any filib type)*/
void to_dat(std::ostream& out, capd::interval const & v);
/** saves as several double values separated by spaces. special case for normal double */
void to_dat(std::ostream& out, double const & v);

/** saves as several double values separated by spaces. Returns as string instead to sending to a stream. */
template<typename AnyTypeSpec>
std::string to_dat(AnyTypeSpec const & item){
	std::ostringstream oss; to_dat(oss, item); return oss.str();
}

/* saves as several double values separated by spaces. */
template<typename VectorSpec>
void to_dat(std::ostream& out, VectorSpec const & v){
	std::string prefix = "";
	for (typename VectorSpec::size_type i = 0; i < v.dimension(); i++){
		out << prefix; to_dat(out, v[i]);
		prefix = " ";
	}
}

/**
 * outputs to the stream the pairs (t, curve(t))
 * for all t = n * h in [t0, t1].
 * It requires curve to be defined on all [t0, t1].
 * It works best with TimePoints that are supporting a discrete
 * equidistant grid (see for example DiscreteGrid type) and
 * h to be compatible with the grid, but in principle it
 * should work with any representation of time.
 */
template<typename TimePointSpec, typename CurveSpec, typename StepSpec>
void value_to_gnuplot(
		std::ostream& out,
		TimePointSpec const &t0,
		TimePointSpec const &t1,
		StepSpec const & h,
		CurveSpec const & curve){
	to_dat(out, typename CurveSpec::RealType(t0)); out << " "; to_dat(out, curve.eval(t0)); out << "\n";
	StepSpec t = StepSpec(t0) + h;
	while (t < StepSpec(t1)){
		try{
			to_dat(out, typename CurveSpec::RealType(t)); out << " "; to_dat(out, curve.eval(t)); out << "\n";
		} catch(...){
			std::cout << "Draw Error at t = " << t << std::endl;
			// do nothing...
		}
		t += h;
	}
	to_dat(out, typename CurveSpec::RealType(t1)); out << " "; to_dat(out, curve.eval(t1)); out << "\n";
}

/**
 * outputs to the stream the pairs (t, curve(t))
 * for all t = n * h in [t0, t1].
 * It requires curve to be defined on all [t0, t1].
 * It works best with TimePoints that are supporting a discrete
 * equidistant grid (see for example DiscreteGrid type) and
 * h to be compatible with the grid, but in principle it
 * should work with any representation of time.
 *
 * It uses dt to shift evaluation, can be used to evaluate over whole intervals
 */
template<typename TimePointSpec, typename CurveSpec, typename StepSpec>
void value_to_gnuplot(
		std::ostream& out,
		TimePointSpec const &t0,
		TimePointSpec const &t1,
		StepSpec const & h,
		CurveSpec const & curve,
		typename CurveSpec::RealType dt){
	to_dat(out, typename CurveSpec::RealType(t0)); out << " "; to_dat(out, curve.getPiece(t0).evalAtDelta(dt)); out << "\n";
	StepSpec t = StepSpec(t0) + h;
	while (t < StepSpec(t1)){
		try{
			to_dat(out, typename CurveSpec::RealType(t)); out << " "; to_dat(out, curve.getPiece(t).evalAtDelta(dt)); out << "\n";
		} catch(...){
			std::cout << "Draw Error at t = " << t << std::endl;
			// do nothing...
		}
		t += h;
	}
//	if (dt == 0.) // TODO: (NOT URGENT?) rethink!
	to_dat(out, typename CurveSpec::RealType(t1)); out << " "; to_dat(out, curve.eval(t1)); out << "\n";
}

template<typename TimePointSpec, typename CurveSpec, typename StepSpec>
void plot_value(
		std::string const& pathprefix,
		TimePointSpec const &t0,
		TimePointSpec const &t1,
		StepSpec const & h,
		CurveSpec const & curve,
		bool live = true){
	auto dirandname = split_pathprefix(pathprefix);
	auto dirpath = dirandname.first;
	auto fileprefix = dirandname.second;
	#ifdef DDES_ALLOW_SYSTEM
	if (dirpath != "") capd::ddeshelper::runSystemCommand(std::string("mkdir -p " + dirpath));
	#endif
	std::string datpath = pathprefix + "ddes-plot.dat"; // TODO: (Not urgent) rethink if add suffix. Maybe check if prefix has filename in it or ends with / ?
	std::string datname = fileprefix + "ddes-plot.dat"; // TODO: (not urgent) check other places in drawing for such behaviour.
	std::string pngname = fileprefix + "ddes-plot.png";
	std::string gppath = pathprefix + "ddes-plot.gp";
	std::string gpname = fileprefix + "ddes-plot.gp";
	std::ofstream outf_sol(datpath);
	capd::ddeshelper::value_to_gnuplot(outf_sol, curve.pastTime(), curve.currentTime(), h, curve);
	outf_sol.close();
	std::ofstream outg(gppath);
	if (!live){
		outg << "set terminal png size 1600,1200" << std::endl;
		outg << "set output '" << pngname << "'" << std::endl;
	}
	outg << "plot '" << datname << "' using 1:($3-$4) with lines, '" << datname << "' using 1:($3+$4) with lines" << std::endl;
	outg.close();
	#ifdef DDES_ALLOW_SYSTEM
	std::ostringstream cmd; cmd << "cd '" << dirpath << "' && gnuplot " << (live ? "-p " : "") << "'" << gpname << "'";
	capd::ddeshelper::runSystemCommand(cmd.str());
	#endif
}

template<typename CurveSpec, typename StepSpec>
void plot_value(
		std::string const& pathprefix,
		StepSpec const & h,
		CurveSpec const & curve,
		bool live = true){
	typename CurveSpec::TimePointType t0 = curve.leftDomain();
	typename CurveSpec::TimePointType t1 = curve.rightDomain();
	plot_value(pathprefix, t0, t1, h, curve, live);
}

template<typename CurveSpec>
void plot_value(
		std::string const& pathprefix,
		CurveSpec const & curve,
		bool live = true){
	typename CurveSpec::TimePointType t0 = curve.leftDomain();
	typename CurveSpec::TimePointType t1 = curve.rightDomain();
	typename CurveSpec::RealType h = (t1 - t0); h /= 128.0; // arbitrarily...
	plot_value(pathprefix, t0, t1, h, curve, live);
}

/**
 * generates apropriate gnuplot files to draw a file under /path/to/a/file[.dat]
 * You give path to a dat file (with optional extension, if different from .dat)
 * it produces files: /path/to/a/file.gp
 * and /path/to/a/file.png (if not live)
 */
void plot_datfile(
		std::string const& pathprefix,
		std::vector<std::string> const& gpoptions={},
		bool live = true,
		std::string filesuffix = "");

/**
 * generates appropriate gnuplot files to draw with gnuplot
 * you give a path and prefix to a file:
 *
 * /a/dir/path/prefix
 *
 * and it creates a file prefix.gp and prefix.png (if live = false) in the directory /a/dir/path/
 *
 * in gpplots you must give full plot info, i.e. "'file.dat' u 1:2 with lines", etc.
 */
void xplot_many(
		std::string const& pathprefix,
		std::string plottype,
		std::vector<std::string> const& gpplots={},
		bool live = true,
		std::string filesuffix = "");

/** xplot_many() but with default mode 'plot' (2d) */
void plot_many(
		std::string const& pathprefix,
		std::vector<std::string> const& gpplots={},
		bool live = true,
		std::string filesuffix = "");

/** xplot_many() but with default mode 'splot' (3d) */
void splot_many(
		std::string const& pathprefix,
		std::vector<std::string> const& gpplots={},
		bool live = true,
		std::string filesuffix = "");

template<typename TimePointSpec, typename CurveSpec, typename StepSpec>
void plot_phasespace(
		std::string& pathprefix,
		TimePointSpec const &t0,
		TimePointSpec const &t1,
		StepSpec const & h,
		CurveSpec const & curve,
		bool live = true,
		std::vector<int> coords = std::vector<int>()){
	int dim = coords.size();
	if (dim == 0) {
		dim = curve.dimension() <= 3 ? curve.dimension() : 3;
		for (int i = 0; i < dim; ++i ) coords.push_back(i);
	}else if (dim > 3){
		dim = 3;
		coords.resize(3);
	}
	if (dim == 1){
		coords.resize(2, 0);
	}
	auto dirandname = split_pathprefix(pathprefix);
	auto dirpath = dirandname.first;
	auto fileprefix = dirandname.second;
	#ifdef DDES_ALLOW_SYSTEM
	if (dirpath != "") capd::ddeshelper::runSystemCommand(std::string("mkdir -p " + dirpath));
	#endif
	std::string datpath = pathprefix + "ddes-phasespace.dat";
	std::string datname = fileprefix + "ddes-phasespace.dat";
	std::string pngname = fileprefix + "ddes-phasespace.png";
	std::string gppath = pathprefix + "ddes-phasespace.gp";
	std::string gpname = fileprefix + "ddes-phasespace.gp";
	std::ofstream outf_sol(datpath);
	capd::ddeshelper::value_to_gnuplot(outf_sol, curve.pastTime(), curve.currentTime(), h, curve);
	outf_sol.close();
	std::ofstream outg(gppath);
	if (!live){
	outg << "set terminal png size 1600,1200" << std::endl;
	outg << "set output '" << pngname << "'" << std::endl;
	}

	std::ostringstream coordsstr;
//	coordsstr << "(0.5*($" << (2*coords[0]+2) << "+$" << (2*coords[0]+3) << "))";
//	for (int i = 1; i < coords.size(); ++i)
//		coordsstr << ":(0.5*($" << (2*coords[i]+2) << "+$" << (2*coords[i]+3) << "))";

	// na pozycji 1 i 2 jest czas oraz jego grubość (zazwyczaj 0)
	// na pozycji 3 i 4 jest wspolrzedna 0 i jej grubosc, etc.
	// TODO: (NOT URGENT) na pozycji 2*coords[i] + 4 jest grubosc
	// TODO: (NOT URGENT) przemyslec, czy ni zsumowac grubosci i np. traktowac jako grubosc punktu na wykresie? Powinno sie dac narysowac liniowo-punktowy wykres i powinno byc widac grubosc
	coordsstr << 2*coords[0]+3;
	for (int i = 1; i < coords.size(); ++i)
		coordsstr << ":" << 2*coords[i]+3;

	if (dim == 2){
		outg << "plot '" << datname << "' using " << coordsstr.str() << " with lines" << std::endl;
	}else {
		outg << "splot '" << datname << "' using " << coordsstr.str() << " with lines" << std::endl;
	}
	outg.close();
	#ifdef DDES_ALLOW_SYSTEM
	std::ostringstream cmd; cmd << "cd '" << dirpath << "' && gnuplot '" << gpname << "'";
	capd::ddeshelper::runSystemCommand(cmd.str());
	#endif

}

template<typename CurveSpec, typename StepSpec>
void plot_phasespace(
		std::string const& pathprefix,
		StepSpec const & h,
		CurveSpec const & curve,
		bool live = true){
	typename CurveSpec::TimePointType t0 = curve.leftDomain();
	typename CurveSpec::TimePointType t1 = curve.rightDomain();
	plot_phasespace(pathprefix, t0, t1, h, curve, live);
}

template<typename CurveSpec>
void plot_phasespace(
		std::string const& pathprefix,
		CurveSpec const & curve,
		bool live = true){
	typename CurveSpec::TimePointType t0 = curve.leftDomain();
	typename CurveSpec::TimePointType t1 = curve.rightDomain();
	typename CurveSpec::RealType h = (t1 - t0); h /= 128.0; // arbitrarily...
	plot_phasespace(pathprefix, t0, t1, h, curve, live);
}

/**
 * outputs to the stream the pairs (t, curve(t))
 * for n equidistant points in [t0, t1].
 * It requires curve to be defined on all [t0, t1].
 * For equidistant timepoints on the grid (e.g. DiscreteGrid)
 * we recommend using other function.
 */
template<typename TimePointSpec, typename CurveSpec>
void value_to_gnuplot(
		std::ostream& out,
		TimePointSpec const &t0,
		TimePointSpec const &t1,
		int n,
		CurveSpec const & curve){
	typedef CurveSpec CurveType;
	typedef typename CurveType::RealType RealType;
	n = n < 0 ? -n : n;
	RealType h = RealType(t1 - t0) / RealType(n);
	value_to_gnuplot(out, t0, t1, h, curve);
}

/**
 * outputs to the stream the pairs (curve(t-tau), curve(t))
 * for all t = n * dt in t0+tau, t1. It requires tau > 0, t1 >= t0,
 * and curve to be defined on all [t0, t1].
 * It works best with TimePoints that are supporting a discrete
 * equidistant grid (see for example DiscreteGrid type),
 * but in principle it should work with any representation of time.
 */
template<typename TimePointSpec, typename CurveSpec, typename StepSpec>
void phasespace_to_gnuplot(
		std::ostream& out,
		TimePointSpec const &t0,
		TimePointSpec const &t1,
		StepSpec const &dt,
		TimePointSpec const &tau,
		CurveSpec const & curve){
	to_dat(out, curve.eval(t0-tau)); out << " ";
	to_dat(out, curve.eval(t0)); out << "\n";
	for (TimePointSpec t = t0 + dt; t < t1; t += dt){
		to_dat(out, curve.eval(t-tau)); out << " ";
		to_dat(out, curve.eval(t)); out << "\n";
	}
	to_dat(out, curve.eval(t1-tau)); out << " ";
	to_dat(out, curve.eval(t1)); out << "\n";
}

/**
 * outputs detailed plots (large) of the representation.
 */
template<typename CurveSpec>
void representation_to_gnuplot(
		std::string const& pathprefix,
		CurveSpec const & curve){
	auto dirandname = split_pathprefix(pathprefix);
	auto dirpath = dirandname.first;
	auto fileprefix = dirandname.second;
	#ifdef DDES_ALLOW_SYSTEM
	if (dirpath != "") capd::ddeshelper::runSystemCommand(std::string("mkdir -p " + dirpath));
	#endif
	// TODO: (NOT URGENT) we assume that curve.begin() has the lowest order
	// TODO: (NOT URGENT) we also assume that all have the same order, so we plot Xi
	typedef typename CurveSpec::size_type size_type;
	typedef typename CurveSpec::VectorType VectorType;
	typedef typename CurveSpec::RealType RealType;
	size_type n = (*curve.begin())->order();
	std::vector<std::string> dat_filenames;
	for (size_type k = 0; k <= n; k++){
		std::ostringstream filename;
		filename << fileprefix << (fileprefix != "" ? "--" : "") << k << "-th--coeff.dat";
		dat_filenames.push_back(filename.str());
		std::ofstream out(dirpath + dat_filenames.back()); out.precision(15);
		for (auto jet = curve.begin(); jet != curve.end(); ++jet){
			VectorType v = (**jet)[k];
			to_dat(out, RealType((*jet)->t0())); out << " "; to_dat(out, v); out << "\n";
		}
		out.close();
	}
	std::ostringstream filename;
	filename << fileprefix << (fileprefix != "" ? "--" : "") << (n+1) << "-th--coeff.dat";
	dat_filenames.push_back(filename.str());
	std::ofstream out(dirpath + dat_filenames.back());  out.precision(15);
	for (auto jet = curve.begin(); jet != curve.end(); ++jet){
		VectorType v = (*jet)->get_Xi();
		to_dat(out, RealType((*jet)->t0())); out << " "; to_dat(out, v); out << "\n";
	}
	out.close();

	std::ostringstream gp_filename;
	std::ostringstream png_filename;
	gp_filename << fileprefix << (fileprefix != "" ? "--" : "") << "repr-plot.gp";
	png_filename<< fileprefix << (fileprefix != "" ? "--" : "") << "repr-plot.png";
	auto gp_filename_str = gp_filename.str();
	std::ofstream gp(dirpath + gp_filename_str);
	gp << "set terminal png size 1200," << (n+2)*600 << " enhanced" << std::endl;
	gp << "set output '" << png_filename.str() << "'" << std::endl;
//	gp << "OPTSHIFT=0" << std::endl;
	gp << "set multiplot layout " << n+2 << ",1" << std::endl;
	size_type k = 0;
	for (auto file = dat_filenames.begin(); file != dat_filenames.end(); ++file, ++k){
		gp << "plot '" << *file << "' using 1:($3-$4):($3+$4) with filledcu notitle, '" << *file << "' using 1:3 with lines title '" << k << "-th'" << std::endl;
	}
	gp.close();

	#ifdef DDES_ALLOW_SYSTEM
	std::string cmd = "";
	cmd += (dirpath != "" ? std::string("cd " + dirpath + " && ") : std::string("")) + "gnuplot " + gp_filename_str;
	capd::ddeshelper::runSystemCommand(cmd);
	#endif
}



} // namespace ddeshelper
} // namespace capd


#endif /* _CAPD_DDEHELPERDRAWING_H_ */
