#define BOOST_TEST_MODULE DDEJetSectionTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <capd/ddes/DDEJetSection.h>
// Include GenericJet.hpp because we use GenericJet with custom types (via MockCurve)
// and DDEJetSection.h only includes GenericJet.h
#include <capd/ddes/storage/GenericJet.hpp>
#include <capd/capdlib.h>
#include <vector>
#include <iostream>

using namespace capd;
using namespace capd::ddes;

// MockCurve declaration
template<typename ScalarT>
struct MockCurve;

// Specialization for double
template<>
struct MockCurve<double> {
    typedef double ScalarType;
    typedef capd::DMatrix MatrixType;
    typedef capd::DVector VectorType;

    typedef capd::ddes::DiscreteTimeGrid<ScalarType> GridType;
    typedef typename GridType::TimePointType TimePointType;
    typedef size_t size_type;

    typedef GenericJet<TimePointType, VectorType, VectorType, MatrixType, false> JetType;

    std::vector<JetType*> m_jets;
    VectorType m_val;

    MockCurve(size_type dim = 1) : m_val(dim) {}

    MockCurve(const MockCurve& other) : m_val(other.m_val) {
        for(auto p : other.m_jets) {
            m_jets.push_back(new JetType(*p));
        }
    }

    MockCurve& operator=(const MockCurve& other) {
        if(this != &other) {
            for(auto p : m_jets) delete p;
            m_jets.clear();
            m_val = other.m_val;
            for(auto p : other.m_jets) {
                m_jets.push_back(new JetType(*p));
            }
        }
        return *this;
    }

    ~MockCurve() {
        for(auto p : m_jets) delete p;
    }

    void addJet(const JetType& jet) {
        m_jets.push_back(new JetType(jet));
    }

    typedef typename std::vector<JetType*>::const_iterator const_iterator;
    typedef typename std::vector<JetType*>::iterator iterator;
    typedef typename std::vector<JetType*>::const_reverse_iterator const_reverse_iterator;
    typedef typename std::vector<JetType*>::reverse_iterator reverse_iterator;

    iterator begin() { return m_jets.begin(); }
    iterator end() { return m_jets.end(); }
    const_iterator begin() const { return m_jets.begin(); }
    const_iterator end() const { return m_jets.end(); }

    reverse_iterator rbegin() { return m_jets.rbegin(); }
    reverse_iterator rend() { return m_jets.rend(); }
    const_reverse_iterator rbegin() const { return m_jets.rbegin(); }
    const_reverse_iterator rend() const { return m_jets.rend(); }

    VectorType getValueAtCurrent() const { return m_val; }
    void setValueAtCurrent(const VectorType& v) { m_val = v; }

    operator VectorType() const { return m_val; }

    template<typename SectionType>
    ScalarType dot(const SectionType& section) const {
        VectorType s = (VectorType)section;
        VectorType v = m_val;
        if(s.dimension() != v.dimension())
             throw std::runtime_error("Dimension mismatch in dot");

        ScalarType res = 0;
        for(size_t i=0; i<v.dimension(); ++i) res += s[i]*v[i];
        return res;
    }

    VectorType get_x() const { return m_val; }

    void operator*=(ScalarType s) {
        m_val *= s;
        if (s == ScalarType(0)) {
             // mock zeroing
        }
    }
};

typedef double ScalarType;
typedef MockCurve<ScalarType> CurveType;
typedef DDEJetSection<CurveType> SectionType;
typedef SectionType::VectorType VectorType;
typedef SectionType::JetType JetType;

BOOST_AUTO_TEST_SUITE(DDEJetSectionTests)

BOOST_AUTO_TEST_CASE(DefaultConstructor) {
    SectionType section;
    BOOST_CHECK_EQUAL(section.dimension(), 0);
    BOOST_CHECK_EQUAL(section.storageDimension(), 0);
}

BOOST_AUTO_TEST_CASE(CoordinateConstructor) {
    size_t d = 3;
    size_t i = 1;
    ScalarType c = 5.0;
    SectionType section(d, i, c);

    BOOST_CHECK_EQUAL(section.dimension(), d);
    BOOST_CHECK_EQUAL(section.get_c(), c);

    VectorType s = section.get_s();
    BOOST_CHECK_EQUAL(s.dimension(), d);
    BOOST_CHECK_EQUAL(s[0], 0.0);
    BOOST_CHECK_EQUAL(s[1], 1.0);
    BOOST_CHECK_EQUAL(s[2], 0.0);

    VectorType orig = (VectorType)section;
    BOOST_CHECK_EQUAL(orig.dimension(), d);
    BOOST_CHECK_EQUAL(orig[1], 1.0);
}

BOOST_AUTO_TEST_CASE(VectorConstructor) {
    size_t d = 2;
    size_t p = 1; // 1 jet
    size_t n = 1; // order 1
    ScalarType c = 2.0;

    VectorType vec(6);
    for(size_t k=0; k<6; ++k) vec[k] = (double)(k+1);

    SectionType section(d, p, n, vec, c);

    BOOST_CHECK_EQUAL(section.get_c(), c);
    BOOST_CHECK_EQUAL(section.dimension(), d);

    VectorType s = section.get_s();
    BOOST_CHECK_EQUAL(s[0], 1.0);
    BOOST_CHECK_EQUAL(s[1], 2.0);

    int jetCount = 0;
    for(auto it = section.begin(); it != section.end(); ++it) {
        jetCount++;
        JetType& jet = *it;
        BOOST_CHECK_EQUAL(jet.order(), n);

        VectorType c0 = (VectorType)jet[0];
        BOOST_CHECK_EQUAL(c0[0], 3.0);
        BOOST_CHECK_EQUAL(c0[1], 4.0);

        VectorType c1 = (VectorType)jet[1];
        BOOST_CHECK_EQUAL(c1[0], 5.0);
        BOOST_CHECK_EQUAL(c1[1], 6.0);
    }
    BOOST_CHECK_EQUAL(jetCount, 1);
}

BOOST_AUTO_TEST_CASE(CurveConstructor) {
    CurveType curve(2);
    curve.m_val[0] = 1.0;
    curve.m_val[1] = 2.0;

    JetType jet(typename CurveType::TimePointType(), 2, 1);
    VectorType v(2); v[0]=3.0; v[1]=4.0;
    jet[0] = v;
    v[0]=5.0; v[1]=6.0;
    jet[1] = v;

    curve.addJet(jet);

    ScalarType c = 10.0;
    SectionType section(curve, c);

    BOOST_CHECK_EQUAL(section.get_c(), c);
    BOOST_CHECK_EQUAL(section.dimension(), 2);

    VectorType s = section.get_s();
    BOOST_CHECK_EQUAL(s[0], 1.0);
    BOOST_CHECK_EQUAL(s[1], 2.0);

    int jetCount = 0;
    for(auto it = section.begin(); it != section.end(); ++it) {
        jetCount++;
        JetType& j = *it;
        VectorType c0 = (VectorType)j[0];
        BOOST_CHECK_EQUAL(c0[0], 3.0);
        BOOST_CHECK_EQUAL(c0[1], 4.0);
    }
    BOOST_CHECK_EQUAL(jetCount, 1);
}

BOOST_AUTO_TEST_CASE(Evaluation) {
    size_t d = 2;
    VectorType vec(2); vec[0] = 1.0; vec[1] = 1.0;
    ScalarType c = 5.0;

    SectionType section(d, 0, 0, vec, c);

    CurveType curve(2);
    curve.m_val[0] = 2.0;
    curve.m_val[1] = 3.0;

    ScalarType val = section(curve);
    BOOST_CHECK_CLOSE(val, 0.0, 1e-10);

    curve.m_val[1] = 4.0;
    val = section(curve);
    BOOST_CHECK_CLOSE(val, 1.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(Gradient) {
    size_t d = 2;
    ScalarType c = 0.0;
    VectorType vec(6);
    for(size_t k=0; k<6; ++k) vec[k] = (double)(k+1);

    SectionType section(d, 1, 1, vec, c);
    CurveType curve(d);
    JetType jet(typename CurveType::TimePointType(), d, 1);
    curve.addJet(jet);

    VectorType grad = section.getGradient(curve);

    BOOST_CHECK_EQUAL(grad.dimension(), d);
    BOOST_CHECK_EQUAL(grad[0], 1.0);
    BOOST_CHECK_EQUAL(grad[1], 2.0);
}

BOOST_AUTO_TEST_CASE(SettersGetters) {
    SectionType section(2, 0);
    section.set_c(10.0);
    BOOST_CHECK_EQUAL(section.get_c(), 10.0);

    VectorType s(2); s[0]=3.0; s[1]=4.0;
    section.set_s(s);
    VectorType s_out = section.get_s();
    BOOST_CHECK_EQUAL(s_out[0], 3.0);
    BOOST_CHECK_EQUAL(s_out[1], 4.0);

    VectorType orig = (VectorType)section;
    BOOST_CHECK_EQUAL(orig[0], 3.0);
    BOOST_CHECK_EQUAL(orig[1], 4.0);
}

BOOST_AUTO_TEST_CASE(Extend) {
    size_t d = 2;
    SectionType section(d, 0, 0.0);

    JetType jet(typename CurveType::TimePointType(), d, 0);
    VectorType v(2); v[0]=9.0; v[1]=8.0;
    jet[0] = v;

    section.extend(jet);

    int jetCount = 0;
    for(auto it = section.begin(); it != section.end(); ++it) jetCount++;
    BOOST_CHECK_EQUAL(jetCount, 1);

    // Explicitly use global VectorType (since it is defined)
    VectorType origVec = (VectorType)section;
    BOOST_CHECK_EQUAL(origVec.dimension(), 4);
    BOOST_CHECK_EQUAL(origVec[0], 1.0);
    BOOST_CHECK_EQUAL(origVec[1], 0.0);
    BOOST_CHECK_EQUAL(origVec[2], 9.0);
    BOOST_CHECK_EQUAL(origVec[3], 8.0);
}

BOOST_AUTO_TEST_CASE(ExceptionTests) {
    BOOST_CHECK_THROW(SectionType(3, 3, 0.0), std::logic_error);

    VectorType vec(5);
    BOOST_CHECK_THROW(SectionType(2, 1, 1, vec, 0.0), std::logic_error);

    SectionType section(2, 0);
    VectorType s(3);
    BOOST_CHECK_THROW(section.set_s(s), std::logic_error);

    JetType jet(typename CurveType::TimePointType(), 3, 0);
    BOOST_CHECK_THROW(section.extend(jet), std::logic_error);
}

BOOST_AUTO_TEST_SUITE_END()
