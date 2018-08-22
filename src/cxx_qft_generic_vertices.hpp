// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef CXX_QFT_GENERIC_VERTICES_H
#define CXX_QFT_GENERIC_VERTICES_H

#include "error.hpp"
#include "numerics2.hpp"

#include <boost/mpl/erase.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/vector.hpp>

#include <complex>
#include <tuple>

namespace flexiblesusy {

namespace cxx_qft {

namespace detail {

template<class Field>
struct number_of_field_indices {
   static constexpr int value =
      std::tuple_size<typename field_indices<Field>::type>::value;
   using type = boost::mpl::int_<value>;
};

template<class Sequence>
struct total_number_of_field_indices {
   using type = typename boost::mpl::fold<
      Sequence,
      boost::mpl::int_<0>,
      boost::mpl::plus<
         boost::mpl::_1,
         number_of_field_indices<
            boost::mpl::_2
            >
         >
      >::type;
   static constexpr int value = type::value;
};

} // namespace detail

/**
 * @class ScalarVertex
 */
class ScalarVertex {
public:
   explicit ScalarVertex(std::complex<double> v)
      : val(v) {}

   std::complex<double> value() const { return val; }

   bool isZero() const {
      return (is_zero(val.real()) && is_zero(val.imag()));
   }

private:
   std::complex<double> val;
};

/**
 * @class ChiralVertex
 */
class ChiralVertex {
public:
   ChiralVertex(const std::complex<double>& left,
                const std::complex<double>& right)
      : value(left, right) {}

   std::complex<double> left() const { return value.first; }
   std::complex<double> right() const { return value.second; }

   bool isZero() const {
      return (is_zero(value.first.real()) && is_zero(value.first.imag()) &&
              is_zero(value.second.real()) && is_zero(value.second.imag()));
   }

private:
   std::pair<std::complex<double>, std::complex<double>> value;
};

/**
 * @class MomentumDifferenceVertex
 */
class MomentumDifferenceVertex {
public:
   MomentumDifferenceVertex(std::complex<double> v, int mi, int si )
      : val(v), minuendIndex(mi), subtrahendIndex(si) {}

   std::complex<double> value(int mi, int si) const {
      if( mi == minuendIndex && si == subtrahendIndex ) {
         return val;
      }

      if( mi == subtrahendIndex && si == minuendIndex ) {
         return -val;
      }

      throw std::invalid_argument(
         "MomentumDifferenceVertex: Wrong index combination" );
      return 0.0;
   }

   bool isZero() const {
      return (is_zero(val.real()) && is_zero(val.imag()));
   }

private:
  std::complex<double> val;
  int minuendIndex;
  int subtrahendIndex;
};

/**
 * @class InverseMetrixVertex
 */
class InverseMetricVertex {
public:
   explicit InverseMetricVertex(std::complex<double> v)
      : val(v) {}

   std::complex<double> value() const { return val; }

   bool isZero() const {
      return (is_zero(val.real()) && is_zero(val.imag()));
   }

private:
  std::complex<double> val;
};

template<class ...Fields>
struct VertexData;

template<class ...Fields>
class Vertex {
  using Data = VertexData<Fields...>;
public:
  using index_bounds = typename boost::mpl::fold<
    boost::mpl::vector<Fields...>,
    boost::mpl::pair<
      boost::mpl::vector<>,
      boost::mpl::vector<>
    >,
    boost::mpl::pair<
      boost::mpl::joint_view<
        boost::mpl::first<boost::mpl::_1>,
        boost::mpl::first<
          cxx_qft::index_bounds<boost::mpl::_2>
        >
      >,
      boost::mpl::joint_view<
        boost::mpl::second<boost::mpl::_1>,
        boost::mpl::second<
          cxx_qft::index_bounds<boost::mpl::_2>
        >
      >
    >
  >::type;
  using indices_type = std::array<
    int,
    detail::total_number_of_field_indices<
      boost::mpl::vector<Fields...>
    >::value
  >;
  using vertex_type = typename Data::vertex_type;

  template<int FieldIndex>
  static typename field_indices<
    typename boost::mpl::at_c<
      boost::mpl::vector<Fields...>,
      FieldIndex
    >::type
  >::type
  field_indices(const indices_type& indices)
  {
    using namespace boost::mpl;
    using fields = vector<Fields...>;

    using result_type = typename cxx_qft::field_indices<
      typename boost::mpl::at_c<fields, FieldIndex>::type
    >::type;

    using preceeding_fields = typename erase<
      fields,
      typename advance<
        typename begin<fields>::type,
        int_<FieldIndex>
      >::type,
      typename end<fields>::type
    >::type;

    constexpr int offset = detail::total_number_of_field_indices<
      preceeding_fields
    >::value;
    constexpr int length = std::tuple_size<result_type>::value;

    result_type result_indices;
    std::copy( indices.begin() + offset,
               indices.begin() + offset + length,
               result_indices.begin() );
    return result_indices;
  }

   template <class EvaluationContext>
   static vertex_type evaluate(const indices_type& indices,
                               const EvaluationContext& context);
};

} // namespace cxx_qft

} // namespace flexiblesusy

#endif