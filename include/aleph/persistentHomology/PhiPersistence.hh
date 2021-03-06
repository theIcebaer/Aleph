#ifndef ALEPH_PERSISTENT_HOMOLOGY_PHI_PERSISTENCE_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_PHI_PERSISTENCE_HH__

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Conversions.hh>
#include <aleph/topology/Intersections.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <initializer_list>
#include <map>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace aleph
{

/**
  Partitions a simplicial complex according to its $\phi$-persistence
  values. This follows the persistent intersection homology algorithm
  in:

    Persistent Intersection Homology\n
    Paul Bendich and John Harer\n

  The function uses a simplicial complex \f$K\f$ and a function \f$\phi\f$
  that determines whether a simplex is proper or not. The function is
  going to create a new simplicial complex. This complex contains all
  proper simplices (in their original order) followed by all improper
  ones.
*/

template <class Simplex, class Function> std::pair<topology::SimplicialComplex<Simplex>, std::size_t> partition( const topology::SimplicialComplex<Simplex>& K, Function phi )
{
  topology::SimplicialComplex<Simplex> L;

  for( auto&& simplex : K )
  {
    if( phi(simplex) )
      L.push_back( simplex );
  }

  auto s = L.size();

  for( auto&& simplex : K )
  {
    if( !phi(simplex) )
      L.push_back( simplex );
  }

  return std::make_pair( L, s );
}

/**
  @class Perversity
  @brief Perversity model in the sense of intersection homology

  Models a perversity in the sense of intersection homology. The class
  ensures that all values satisfy

  \f
    -1 \leq p_k \leq k-1
  \f
*/

class Perversity
{
public:

  /**
    Creates a new perversity from a range of values. The values must be
    at least implicitly convertible to integers.
  */

  template <class InputIterator> Perversity( InputIterator begin, InputIterator end )
    : _values( begin, end )
  {
    for( std::size_t k = 0; k < _values.size(); k++ )
    {
      if( _values[k] < -1 )
        _values[k] = -1;
      // There is an index shift going on here: since $k$ runs from $0$
      // to $d-1$, there is no need to shift the upper bound.
      else if( _values[k] > static_cast<int>( k ) )
        _values[k] = static_cast<int>( k );
    }
  }

  /** Creates a new perversity from an initializer list of values */
  Perversity( std::initializer_list<int> values )
    : Perversity( values.begin(), values.end() )
  {
  }

  /**
    Queries the perversity value in a given dimension $d$. Invalid
    dimension values only cause the function to return a zero.
  */

  int operator()( std::size_t d ) const noexcept
  {
    if( d < _values.size() + 1 )
      return _values[ static_cast<std::size_t>( d-1 ) ];
    else
      return 0;
  }

  using const_iterator = typename std::vector<int>::const_iterator;

  const_iterator begin() const noexcept { return _values.begin(); }
  const_iterator end()   const noexcept { return _values.end();   }

private:
  std::vector<int> _values;
};

std::ostream& operator<<( std::ostream& o, const Perversity p )
{
  o << "[";

  for( auto it = p.begin(); it != p.end(); ++it )
  {
    if( it != p.begin() )
      o << ",";

    o << *it;
  }

  o << "]";

  return o;
}

template <class Simplex> auto calculateIntersectionHomology( const aleph::topology::SimplicialComplex<Simplex>& K,
                                                             const std::vector< aleph::topology::SimplicialComplex<Simplex> >& X,
                                                             const Perversity& p ) -> std::vector< PersistenceDiagram<typename Simplex::DataType> >
{
  // 0. Check consistency of strata
  // 1. Create allowability function based on the dimensionality of the
  //    intersection of simplices with individual strata.
  // 2. Calculate $phi$-persistence
  // 3. Convert the result into a persistence diagram.

  // Check consistency of filtration -----------------------------------
  //
  // The maximum dimension of each complex in the filtration has to
  // match the dimension of the simplicial complex.

  {
    std::size_t minDimension = K.dimension();
    std::size_t maxDimension = 0;

    for( auto&& x : X )
    {
      if( !K.empty() )
      {
        minDimension = std::min( minDimension, x.dimension() );
        maxDimension = std::max( maxDimension, x.dimension() );
      }
    }

    if( maxDimension != K.dimension() )
      throw std::runtime_error( "Invalid filtration" );
  }

  // Check whether simplex is allowable --------------------------------

  std::map<Simplex, bool> phi;

  {
    auto d = K.dimension();

    for( auto&& s : K )
    {
      bool admissible = true;

      for( std::size_t k = 1; k <= d; k++ )
      {
        // The notation follows Bendich and Harer, so $i$ is actually
        // referring to a dimension instead of an index. Beware!
        auto i            = s.dimension();
        auto intersection = aleph::topology::lastLexicographicalIntersection( X.at( d - k ), s );
        auto dimension    = intersection.empty() ? -1 : static_cast<long>( intersection.dimension() );
        admissible        = admissible && intersection.empty() ? true : static_cast<long>( dimension ) <= ( long(i) - long(k) + long( p(k) ) );
      }

      phi[s] = admissible;
    }
  }

  // Partition according to allowable simplices ------------------------

  aleph::topology::SimplicialComplex<Simplex> L;
  std::size_t s = 0;

  std::tie( L, s ) =
    aleph::partition( K, [&phi] ( const Simplex& s )
                         {
                           return phi.at(s);
                         } );

  // Calculate persistent intersection homology ------------------------

  auto boundaryMatrix             = aleph::topology::makeBoundaryMatrix( L, s );
  using IndexType                 = typename decltype(boundaryMatrix)::Index;
  bool includeAllUnpairedCreators = true;
  auto pairing                    = aleph::calculatePersistencePairing( boundaryMatrix, includeAllUnpairedCreators, static_cast<IndexType>(s) );
  auto persistenceDiagrams        = aleph::makePersistenceDiagrams( pairing, L );

  return persistenceDiagrams;
}

} // namespace aleph

#endif
