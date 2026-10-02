// -----------------------------------------------------------------------------------------------------------------------
// Andrew W. Rose, 2026
// Imperial College London Particles community
// -----------------------------------------------------------------------------------------------------------------------

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <regex>


template < typename Derived >
struct magic
{
  template< typename U > void Apply( const U& aFn )
  {
    static_cast<       Derived* >( this )->Implement( aFn );
  }

  template< typename U > void Apply( const U& aFn ) const
  {
    static_cast< const Derived* >( this )->Implement( aFn );
  }

  const std::vector< std::string >& MagicFields() const
  {
    static const std::vector<std::string> v = MagicDelimeter( static_cast< const Derived* >( this )->MagicString() );
    return v;
  }

  bool operator== ( const magic< Derived >& aB ) const
  {
    bool lBool;
    aB.Apply( [&]( auto&... params ){ lBool = MagicCompare( params... ); } ); // Variadic lambda
    return lBool;
  }

private:
  static std::vector< std::string > MagicDelimeter( const std::string& aStr )
  {
    const std::regex lRegex( "\\s*,\\s*" );
    return std::vector<std::string>( std::sregex_token_iterator( aStr.begin(), aStr.end(), lRegex , -1 ), std::sregex_token_iterator() );
  }

  template< typename... U >
  bool MagicCompare ( const U&... aB ) const
  {
    bool lBool;
    Apply( [&]( auto&&... params ){ lBool = ( (params==aB) && ... ); } ); // Variadic lambda invoking a C++17 fold-expression
    return lBool;
  }

};


template< typename T >
std::ostream& operator<< ( std::ostream& aStr , const magic< T >& aArg )
{
  aStr << "(";
  auto lIt = aArg.MagicFields().begin();
  aArg.Apply( [&]( auto&&... params ){ ( ( aStr << ' ' << *lIt++ << '=' << params ) , ... ); } ); // Variadic lambda invoking a C++17 fold-expression
  aStr << " )";
  return aStr;
}

// template< typename T >
// std::istream& operator>> ( std::istream& aStr , magic< T >& aArg )
// {
  // aArg.Apply( [ &aStr ]( auto&&... params ){ ( ( aStr >> params ) , ... ); } ); // Variadic lambda invoking a C++17 fold-expression
  // return aStr;
// }


#define MAGIC( ... ) \
  template< typename __________U__________ > void Implement( const __________U__________& aFn )       { aFn( __VA_ARGS__ ); } \
  template< typename __________U__________ > void Implement( const __________U__________& aFn ) const { aFn( __VA_ARGS__ ); } \
  inline std::string MagicString() const { return #__VA_ARGS__ ; }


#define MAGIC_CONSTRUCTOR( CLASS ) \
  template< typename... U > CLASS( const U&... aArgs ){ Apply( [&]( auto&&... params ){ ( (params=aArgs) , ... ); } ); }
