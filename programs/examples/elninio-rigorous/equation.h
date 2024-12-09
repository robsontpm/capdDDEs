#ifndef ELNINIO_RIG_H_
#define ELNINIO_RIG_H_

#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>

template<class T, capd::filib::RoundingStrategy R, capd::filib::IntervalMode M> 
struct fadbad::Op<capd::filib::Interval<T, R, M>>
{
  typedef capd::filib::Interval<T, R, M> Base;

  static Base myInteger(const int i) { return Base(i); }
  static Base myZero() { return myInteger(0); }
  static Base myOne() { return myInteger(1);}
  static Base myTwo() { return myInteger(2); }
  static Base myPI() { return Base::pi(); }
  static Base myPos(const Base& x) { return +x; }
  static Base myNeg(const Base& x) { return -x; }
  template <typename U> static Base& myCadd(Base& x, const U& y) { return x+=y; }
  template <typename U> static Base& myCsub(Base& x, const U& y) { return x-=y; }
  template <typename U> static Base& myCmul(Base& x, const U& y) { return x*=y; }
  template <typename U> static Base& myCdiv(Base& x, const U& y) { return x/=y; }
  static Base myInv(const Base& x) { return myOne()/x; }
  static Base mySqr(const Base& x) { return sqr(x); }
  template <typename X, typename Y>
  static Base myPow(const X& x, const Y& y) { return power(x,y); }
  static Base mySqrt(const Base& x) { return sqrt(x); }
  static Base myLog(const Base& x) { return log(x); }
  static Base myExp(const Base& x) { return exp(x); }
  static Base mySin(const Base& x) { return sin(x); }
  static Base myCos(const Base& x) { return cos(x); }
  static Base myBasean(const Base& x) { return tan(x); }
  static Base myAsin(const Base& x) { return asin(x); }
  static Base myAcos(const Base& x) { return acos(x); }
  static Base myAtan(const Base& x) { return atan(x); }
  static bool myEq(const Base& x, const Base& y) { return x==y; }
  static bool myNe(const Base& x, const Base& y) { return x!=y; }
  static bool myLt(const Base& x, const Base& y) { return x<y; }
  static bool myLe(const Base& x, const Base& y) { return x<=y; }
  static bool myGt(const Base& x, const Base& y) { return x>y; }
  static bool myGe(const Base& x, const Base& y) { return x>=y; }
};


template<typename T> inline T get_pi(){ return T::pi(); }

template<> inline double get_pi(){ return 3.14159265358979323846; }

template<typename T> inline T tanh(T const& x){
	return (exp(2 * x) - 1) / (exp(2 * x) + 1);
}

template<typename ScalarSpec, typename ParamSpec = ScalarSpec>
class ElNinoRig{
public:
	typedef ScalarSpec ScalarType;
	typedef unsigned int size_type;
	typedef ScalarType RealType;
	typedef ParamSpec ParamType;
	typedef capd::vectalg::Vector<ParamSpec, 0> ParamsVectorType;
	typedef capd::vectalg::Vector<ScalarType, 0> VectorType;

	ElNinoRig(): alpha(1), beta(1), kappa(1) {}
	ElNinoRig(ParamType alpha, ParamType beta, ParamType kappa): alpha(alpha), beta(beta), kappa(kappa) {}
	ElNinoRig(capd::vectalg::Vector<ParamSpec, 0> const & params): alpha(params[0]), beta(params[1]), kappa(params[2]) {}

	static size_type imageDimension() { return 1; };

	static size_type dimension() { return 2; }

	static size_type getParamsCount() { return 3; }

	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		typedef typename OutVectorSpec::ScalarType OutScalarType;
		auto Tt = x[0];
		auto Ttau = x[1];
		fx[0] = -alpha * tanh(kappa * Ttau) + beta * cos(2 * get_pi<ScalarType>() * t);
	}

	static std::string show(){
		return "El Nino equation as autonomous system: $T'(t) = \\beta \\cos(2 \\pi s) - \\alpha \\tanh( \\kappa T(t-\\tau))$.";
	}

protected:
	ParamType alpha;
	ParamType beta;
	ParamType kappa;
};


#endif /* ELNINIO_RIG_H_ */

