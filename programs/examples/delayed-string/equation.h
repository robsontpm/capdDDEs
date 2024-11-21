#pragma once

#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>
#include <capd/ddeshelper/DDEHelperNonrigorous.hpp>

template<class ScalarSpec, class ParamSpec = ScalarSpec>
class DelayedStringEq {
public:
  using ScalarType = ScalarSpec;
  using size_type = unsigned int;
  using RealType = ScalarType;
  using ParamType = ParamSpec;
  using ParamsVectorType = capd::vectalg::Vector<ParamType, 0>;
  using VectorType = capd::vectalg::Vector<ScalarType, 0>;

  DelayedStringEq(const ParamsVectorType& params):
    m_alpha(params[0])
  {}

  static size_type imageDimension() {
    return 2;
  }

  static size_type dimension() {
    return 4;
  }

  static size_type getParamsCount() {
    return 1;
  }

  template<class RealSpec, class InVectorSpec, class OutVectorSpec>
  void operator()(const RealSpec& t, const InVectorSpec& x, OutVectorSpec& fx) const {
    /*const auto& T0 = x[0];*/
    /*const auto& T1 = x[1];*/
    /*const auto& Ttau0 = x[2];*/
    /*const auto& Ttau1 = x[3];*/
    fx[0] = x[1];
    fx[1] = -m_alpha * x[2];
  }

private:
  ParamType m_alpha;
};
