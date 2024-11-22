#pragma once

#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>
#include <capd/ddeshelper/DDEHelperNonrigorous.hpp>

/**
 * @brief A class representing an equation \f$ x''(t) = -\alpha x(t - \tau) \f$
 * (Hooke's law)
 *
 * @details
 * An example of a minimal class implementation that can be used with
 * @ref capd::ddeshelper::NonrigorousHelper. It specifies a function
 * \f$ f:\mathbb{R} \times \mathbb{R}^4 \to \mathbb{R}^2 \f$:
 * \f\[
 *   f(t, (u, v, w, z)) = (v, -\alpha w) 
 * \f\]
 * When viewed in the context of delay differential equation
 * \f$ (x_1'(t), x_2'(t)) = f(t, (x_1(t), x_2(t), x_1(t - \tau), x_2(t - \tau))) \f$
 * arguments \f$ u \f$, \f$ w \f$ represent the values of function
 * \f$ x_1 \f$ at times \f$ t \f$ and \f$ t - \tau \f$ respectively.
 * Similarily, arguments \f$ v \f$, \f$ z \f$ represent the values
 * of function \f$ x_2 \f$ at times \f$ t \f$ and \f $ t - \tau \f$
 *
 * @tparam ScalarSpec A type representing the scalars over which the vector
 * field is defined. (NOTE: currently assumed to be real)
 *
 * @tparam ParamSpec Type of the optional parameters in the function
 */
template<class ScalarSpec, class ParamSpec = ScalarSpec>
class DelayedSpring {
public:
  /** @brief Type of the scalars in the vector field */
  using ScalarType = ScalarSpec;

  using size_type = unsigned int;

  /** @brief Type representing real numbers. Used by the time variable */
  using RealType = ScalarType;

  /** @brief Type of parameters in the equation */
  using ParamType = ParamSpec;

  /** @brief Preffered way to store a collection of parameters */
  using ParamsVectorType = capd::vectalg::Vector<ParamType, 0>;

  /** @brief Type used to store return value of the function */
  using VectorType = capd::vectalg::Vector<ScalarType, 0>;

  /** @brief Constructor from parameter vector
   *  
   *  @details
   *  Creates a new object from parameters specified in the @ref params
   *  vector. The first @ref getParamsCount() elements in the vector represent
   *  parameters. @ref params can have more elements, but they should be ignored
   *
   *  @param params Vector of parameters. params.size() >= @ref getParamsCount()
   */
  DelayedSpring(const ParamsVectorType& params):
    m_alpha(params[0])
  {}

  /** @brief The dimension of the image of the function represented by
   *  the class.
   *
   *  @details
   *  Returns the expected number of elements in the output of the
   *  @ref operator(). It is equal to the dimension of equation
   *  we want to solve.
   *
   *  @code
   *  using Eq = DelayedSpring<double>;
   *  Eq eq(Eq::ParamsVectorType{ 1.0 });
   *  Eq::VectorType x{ 1.0, 2.0, 3.0, 4.0 };
   *  Eq::VectorType fx(2);
   *  eq(0.0, x, fx);
   *  assert(fx.dimension() == Eq::imageDimension());
   *  @endcode
   *
   *  @return The dimension of the image.
   */
  static size_type imageDimension() {
    return 2;
  }

  /** @brief The dimension of the domain of the function represented
   *  the class (without the time).
   *  
   *  @details
   *  Returns the expected length of the argument vector in @ref operator().
   *  It is equal to the number of delays multiplied by the dimension of the
   *  image.
   *  
   *  @code
   *  using Eq = DelayedSpring<double>;
   *  Eq eq(Eq::ParamsVectorType{ 1.0 });
   *  Eq::VectorType x{ 1.0, 2.0, 3.0, 4.0 };
   *  Eq::VectorType fx(2);
   *  eq(0.0, x, fx);
   *  assert(x.dimension() == Eq::dimension());
   *  @endcode
   *
   *  @return The dimension of the domain.
   */
  static size_type dimension() {
    return 4;
  }

  /** @brief The number of parameters used by the function represented
   *  by the class.
   *
   *  @return The nubmer of parameters in the function.
   */
  static size_type getParamsCount() {
    return 1;
  }

  /** @brief Calculates the value of the function represented by the class.
   *  
   *  @details
   *  Calculates the value of the function for time argument @ref t
   *  and argument @ref x and inserts the result into @ref fx.
   *  The structure of @ref x is the following:
   *
   *  When used in the equation \f$ x'(t) = f(t, x(t), x(t - \tau)) \f$
   *  the structure of @ref x is going to be:
   *
   *  x[0], x[1], ..., x[DelayedSpring::imageDimension() - 1] are the
   *  coordinates \f$ x \f$ at time \f$ t \f$.
   *
   *  x[DelayedSpring::imageDimension()], ..., x[DelayedSpring::dimension() - 1]
   *  are the coordinates of \f$ x \f$ at time \f$ t - \tau \f$. (Note, that
   *  the equation in the example has only one delay, so
   *  DelayedSpring::dimension() == 2 * DelayedSpring::imageDimension())
   *
   *  @param [in] t Time argument
   *
   *  @param [in] x Argument.
   *  Satisfies x.dimension() == DelayedSpring::dimension()
   *
   *  @param [out] fx The value of the function.
   *  Satisfies fx.dimension() == DelayedSpring::imageDimension()
   */
  template<class RealSpec, class InVectorSpec, class OutVectorSpec>
  void operator()(const RealSpec& t, const InVectorSpec& x, OutVectorSpec& fx) const {
    // x[0] = x_1(t)
    // x[1] = x_2(t)
    // x[2] = x_1(t - \tau)
    // x[3] = x_2(t - \tau)
    fx[0] = x[1];
    fx[1] = -m_alpha * x[2];
  }

private:
  ParamType m_alpha;
};
