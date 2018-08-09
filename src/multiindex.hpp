#ifndef H_FS_MULTIINDEX
#define H_FS_MULTIINDEX

#include "find_if.hpp"

#include <array>
#include <iterator>

#include <boost/array.hpp>

#include <boost/mpl/at.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/range_c.hpp>

#include <boost/range/iterator_range.hpp>

#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/include/copy.hpp>
#include <boost/fusion/adapted/boost_array.hpp>

namespace detail
{
template<class End, class Location>
struct can_increment_multiindex
{
  template<class Indices>
  bool operator()( const Indices &indices )
  {
    return (indices[Location::value] !=
      boost::mpl::at<End, Location>::type::value - 1);
  }
};

template<class End>
struct increment_checker
{
  template<class Location>
  using type = can_increment_multiindex<End, Location>;
};

template<class Begin>
struct increment_multiindex
{
  template<class Location>
  struct type
  {
    template<class Indices>
    void operator()( Indices &indices )
    {
      boost::array<
        int,
        std::tuple_size<Indices>::value
      > boost_temp;
      boost::fusion::copy( Begin(), boost_temp );
      std::copy( boost_temp.begin(),
        boost_temp.begin() + Location::value,
        indices.begin() );

      ++(indices[Location::value]);
    }
  };
};

template<class Begin, class End, bool>
class multiindex_impl;

/* FIXME: This copying back and forth
          using boost::array<> between
          boost::mpl and std::array<>
          is stupid...
          Using boost::hana would be the
          way to go! */
template<class Begin, class End>
class multiindex_impl<Begin, End, false>
{
  static constexpr auto NumIndices = boost::mpl::size<Begin>::value;
  using data_type = std::array<int, NumIndices>;
    
  data_type data;

  multiindex_impl( data_type &&d )
  : data( std::move( d ) ) {}
public:
  static multiindex_impl begin()
  {
    boost::array<int, NumIndices> boost_temp;
    std::array<int, NumIndices> std_temp;

    boost::fusion::copy( Begin(), boost_temp );
    std::copy( boost_temp.begin(), boost_temp.end(), std_temp.begin() );

    return multiindex_impl{ std::move( std_temp ) };
  }

  static multiindex_impl end()
  {
    boost::array<int, NumIndices> boost_temp;
    std::array<int, NumIndices> std_temp;

    boost::fusion::copy( End(), boost_temp );
    std::copy( boost_temp.begin(), boost_temp.end(), std_temp.begin() );

    return multiindex_impl{ std::move( std_temp ) };
  }
  
  multiindex_impl &operator++()
  {
    using locations = boost::mpl::range_c<int, 0, NumIndices>;

    if( meta::find_if<
          locations,
          detail::increment_checker<End>::template type,
          detail::increment_multiindex<Begin>::template type
        >( data ) == false )
    {
      boost::array<int, NumIndices> boost_temp;
      boost::fusion::copy( End(), boost_temp );

      std::copy( boost_temp.begin(), boost_temp.end(), data.begin() );
    }
    
    return *this;
  }

  multiindex_impl operator++( int )
  {
     multiindex_impl copy( *this );
     operator++();
     return copy;
  }

  const data_type &operator*( void ) const
  { return data; }

  const data_type *operator->( void ) const
  { return &data; }

  bool operator==( const multiindex_impl &other ) const
  { return data == other.data; }

  bool operator!=( const multiindex_impl &other ) const
  { return !(*this == other); }
};

template<class Begin, class End>
class multiindex_impl<Begin, End, true>
{
  using data_type = std::array<int, 0>;
  bool is_incremented;

  data_type data;

  multiindex_impl( bool is_inc )
  : is_incremented( is_inc ) {}
public:
  static multiindex_impl begin()
  { return multiindex_impl{ false }; }

  static multiindex_impl end()
  { return multiindex_impl{ true }; }
  
  multiindex_impl &operator++()
  {
    is_incremented = true;
    return *this;
  }

  multiindex_impl operator++( int )
  {
     multiindex_impl copy( *this );
     operator++();
     return copy;
  }

  const data_type &operator*( void ) const
  { return data; }

  const data_type *operator->( void ) const
  { return &data; }

  bool operator==( const multiindex_impl &other ) const
  { return is_incremented == other.is_incremented; }

  bool operator!=( const multiindex_impl &other ) const
  { return !(*this == other); }
};
}

template<class Begin, class End>
class multiindex
: public detail::multiindex_impl<
    Begin,
    End,
    std::is_same<Begin, End>::type::value
>
{
  using impl = detail::multiindex_impl<
    Begin,
    End,
    std::is_same<Begin, End>::type::value
  >;

  multiindex( impl &&index )
  : impl( std::move( index ) ) {}
public:
  static multiindex begin()
  { return impl::begin(); }

  static multiindex end()
  { return impl::end(); }
};

namespace std
{
  template<class Begin, class End>
  struct iterator_traits<multiindex<Begin, End>>
  {
    using difference_type = std::ptrdiff_t;
    using value_type = typename std::decay<decltype(
      *std::declval<multiindex<Begin, End>>() )
    >::type;
    using pointer = const value_type *;
    using reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;
  };
}

#endif

