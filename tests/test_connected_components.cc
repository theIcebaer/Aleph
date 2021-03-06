#include <aleph/config/Base.hh>

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/RipsSkeleton.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <tests/Base.hh>

#include <aleph/persistentHomology/Calculation.hh>
#include <aleph/persistentHomology/ConnectedComponents.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <vector>

using namespace aleph::containers;
using namespace aleph::geometry;
using namespace aleph::topology;
using namespace aleph;

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Point cloud loading" );

  using PointCloud = PointCloud<T>;
  using Distance   = aleph::distances::Euclidean<T>;

  PointCloud pointCloud = load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_colon_separated.txt" ) );

  ALEPH_ASSERT_THROW( pointCloud.size()      == 150 );
  ALEPH_ASSERT_THROW( pointCloud.dimension() ==   4);

  ALEPH_TEST_END();

  ALEPH_TEST_BEGIN( "Rips skeleton calculation" );

  using Wrapper      = BruteForce<PointCloud, Distance>;
  using RipsSkeleton = RipsSkeleton<Wrapper>;

  Wrapper wrapper( pointCloud );
  RipsSkeleton ripsSkeleton;

  auto K = ripsSkeleton( wrapper, 1.0 );


  using Simplex = typename decltype(K)::ValueType;

  ALEPH_ASSERT_THROW( K.empty() == false );
  ALEPH_ASSERT_THROW( std::count_if( K.begin(), K.end(), [] ( const Simplex& s ) { return s.dimension() == 1; } ) != 0 );

  ALEPH_TEST_END();

  ALEPH_TEST_BEGIN( "Zero-dimensional persistent homology calculation" );

  K.sort( filtrations::Data<Simplex>() );

  auto diagrams = calculatePersistenceDiagrams( K );
  auto diagram2 = std::get<0>( calculateZeroDimensionalPersistenceDiagram( K ) );

  ALEPH_ASSERT_THROW( diagrams.empty() == false );

  auto diagram1 = diagrams.front();

  ALEPH_ASSERT_THROW( diagram1.empty() == false );
  ALEPH_ASSERT_EQUAL( diagram1.size(), pointCloud.size() );

  diagram1.removeDiagonal();

  ALEPH_ASSERT_THROW( diagram2.empty() == false );
  ALEPH_ASSERT_EQUAL( diagram1.size(), diagram2.size() );

  using Point = typename decltype(diagram1)::Point;

  auto sortPoints = [] ( const Point& p, const Point& q )
  {
    return p.x() < q.x() || ( p.x() == q.x() && p.y() < q.y() );
  };

  std::sort( diagram1.begin(), diagram1.end(), sortPoints );
  std::sort( diagram2.begin(), diagram2.end(), sortPoints );

  ALEPH_ASSERT_THROW( diagram1 == diagram2 );

  ALEPH_TEST_END();
}

int main()
{
  test<float> ();
  test<double>();
}
