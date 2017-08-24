#ifndef ALEPH_UTILITIES_CONTAINER_OPERATORS_HH__
#define ALEPH_UTILITIES_CONTAINER_OPERATORS_HH__

/**
  @file  ContainerOperators.hh
  @brief Defines (arithmetic) operators to use on standard containers

  In some cases, it makes sense to treat an STL container as
  a mathematical vector. These objects then should support a
  set of arithmetic operations, such as addition. As the STL
  does not provide these things, we have to do them manually
  in this file.

  Note that the operations are kept in *aleph::utilities* in
  order not to pollute the global namespace.
*/

#include <algorithm>
#include <iterator>
#include <vector>

#include <cassert>

namespace aleph
{

namespace utilities
{

/**
  Element-wise addition of two vectors of the same size.

  @param lhs First vector
  @param rhs Second vector

  @returns Addition of the two vectors
*/

template <class T> std::vector<T> operator+( const std::vector<T>& lhs, const std::vector<T>& rhs )
{
  assert( lhs.size() == rhs.size() );

  std::vector<T> result;
  result.reserve( lhs.size() );

  std::transform( lhs.begin(), lhs.end(), rhs.begin(), std::back_inserter( result ), std::plus<T>() );
  return result;
}

/**
  Element-wise subtraction of two vectors of the same size.

  @param lhs First vector
  @param rhs Second vector

  @returns Addition of the two vectors
*/

template <class T> std::vector<T> operator-( const std::vector<T>& lhs, const std::vector<T>& rhs )
{
  assert( lhs.size() == rhs.size() );

  std::vector<T> result;
  result.reserve( lhs.size() );

  std::transform( lhs.begin(), lhs.end(), rhs.begin(), std::back_inserter( result ), std::minus<T>() );
  return result;
}

/**
  Element-wise multiplikation of a vector with a scalar. The scalar type S and 
  the vector-element type T must support multiplikation with result form type 
  T.
  
  @param scalar scalar to multiplikate to the vector elements
  @param vector vector
  
  @returns Vector of multiplikations of scalar and vector-elements
 */

template <typename T, typename S>
std::vector<T> operator* (const S& scalar, const std::vector<T>& vector)
{
  std::vector<T> result;
  for (auto x : vector)
  {
    result.push_back(scalar * x);
  }
  return result;
}

/**
  Element-wise addition from elements of the right hand side vector on the 
  elements of the left hand side vector. 
  
  @param lhs vector to add the elements of the other vector on
  @param vector vector with the added elements
  
  @returns reference to the accumulated vector lhs
 */

template <typename T>
std::vector<T>& operator+= (std::vector<T>& lhs, const std::vector<T>& rhs)
{
  assert(lhs.size() == rhs.size());
  size_t N = lhs.size();
  for (size_t i = 0; i < N; i++)
  {
    lhs[i] += rhs[i];
  }
  return lhs;
}

/**
  Casts an Object of type vector<vector<pair>> to an object of type pair
  <vector<vector>>. 
  
  @param input vector<vector<pair>> object to cast
  
  @returns pair<vector<vector>> object
 */

// TODO: probably unnecessary think about it.
template<typename T>
std::pair< std::vector<std::vector<T>>, std::vector<std::vector<T>> > cast_pair(std::vector<std::vector<std::pair<T,T> >> input)
{
  std::vector<std::vector<T>> X;
  std::vector<std::vector<T>> Y;
  for (size_t k = 0; k < input.size(); k++)
  {
    X.push_back(std::vector<T>() );
    Y.push_back(std::vector<T>() );
    
    for (auto critPoint : input[k])
    {
      X[k].push_back(critPoint.first);
      Y[k].push_back(critPoint.second);
    }
  }
  return std::make_pair(X,Y);
}  
  
} // namespace utilities

} // namespace aleph



#endif
