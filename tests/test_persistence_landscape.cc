#include <tests/Base.hh>
#include <aleph/persistenceDiagrams/PersistenceLandscape.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/containers/PointCloud.hh>
#include <aleph/geometry/SphereSampling.hh>
#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/distances/Euclidean.hh>

#include <cmath>
#include <deque>
#include <limits>
#include <random>

using DataType              = double;
using Distance              = aleph::distances::Euclidean<DataType>;
using PointCloud            = aleph::containers::PointCloud<DataType>;
using PersistenceDiagram    = aleph::PersistenceDiagram<DataType>;
using PersistenceLandscape  = aleph::PersistenceLandscape<>;
using Interval              = aleph::Interval<double>;



template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
    for (auto& el : vec)
    {
        os << el << ' ';
    }
    return os;
}

template <class T> 
aleph::PersistenceDiagram<T> createRandomPersistenceDiagram( unsigned n )
{
  std::random_device rd;
  std::default_random_engine rng( rd() );
  std::uniform_real_distribution<T> 
    distribution( T(0), 
                  T( std::nextafter( T(1), std::numeric_limits<T>::max() ) )
                );

  PersistenceDiagram D;

  for( unsigned i = 0; i < n; i++ )
  {
    auto x = distribution( rng );
    auto y = distribution( rng );

    if( x > y )
      std::swap( x,y );

    D.add( x,y );
  }

  return D;
}

/*
 * This function tests the basic functionalities of the PersistenceLandscape 
 * class.
 */

void __attribute__((optimize("O0")))testPersistenceLandscape()
{
  double inf = std::numeric_limits<double>::infinity();
  
  PersistenceLandscape l0;
  
  assert( l0[0](0) == 0 );
  assert( l0[0](-12) == 0 );
  assert( l0[0](15) == 0 );
  assert( l0[0](0.5) == 0 );
  assert( l0[0](1.21324) == 0 );
  
  std::deque<Interval> data_1 = 
    { Interval(1. , 5.), Interval(2. , 8.), 
      Interval(3. , 4.), Interval(5. , 9.),
      Interval(6. , 7.), Interval(9.,10.) 
    };
      
  std::deque<Interval> data_2 = 
    { Interval(1. , 6.), Interval(2. , 3.),
      Interval(3. , 7.), Interval(4. , 5.) 
    };
    
  std::deque<Interval> twisted_data_1 = 
    { Interval(2. , 8.), Interval(1. , 5.),
      Interval(3. , 4.), Interval(9.,10.) , 
      Interval(5. , 9.), Interval(6. , 7.)
    };

  std::deque<Interval> doubled_data_1 = 
    { Interval(2. , 8.), Interval(1. , 5.),
      Interval(3. , 4.), Interval(9.,10.) , 
      Interval(5. , 9.), Interval(6. , 7.),
      Interval(2. , 8.), Interval(1. , 5.),
      Interval(3. , 4.), Interval(9.,10.) 
    };    
    
  PersistenceLandscape l1(data_1);
  PersistenceLandscape l2(data_2);
  PersistenceLandscape twisted_l1(twisted_data_1);
  PersistenceLandscape doubled_l1(doubled_data_1);
  
  assert( l1.getX() == twisted_l1.getX() );
  assert( l1.getY() == twisted_l1.getY() );
  
  assert( l1.getX() == std::vector<std::vector<double>>(
    {{-inf, 1.0, 3.0, 3.5, 5.0, 6.5, 7.0, 9.0, 9.5, 10, inf}, 
     {-inf, 2.0, 3.5, 5.0, 6.5, 8.0, inf},
     {-inf, 3.0, 3.5, 4.0, 6.0, 6.5, 7.0, inf}} ) );
  assert( l1.getY() == std::vector<std::vector<double>>(
    {{ 0.0, 0.0, 2.0, 1.5, 3.0, 1.5, 2.0, 0.0, 0.5, 0, 0.0 }, 
     { 0.0, 0.0, 1.5, 0.0, 1.5, 0.0, 0.0 },
     { 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0 }} ) );
  
  assert( l2.getX() == std::vector<std::vector<double>>(
    {{-inf, 1.0, 3.5, 4.5, 5.0, 7.0, inf}, 
     {-inf, 2.0, 2.5, 3.0, 4.5, 6.0, inf},
     {-inf, 4.0, 4.5, 5.0, inf}} ) );
  assert( l2.getY() == std::vector<std::vector<double>>(
    {{ 0.0, 0.0, 2.5, 1.5, 2.0, 0.0, 0.0 }, 
     { 0.0, 0.0, 0.5, 0.0, 1.5, 0.0, 0.0 },
     { 0.0, 0.0, 0.5, 0.0, 0.0 }} ) );
  

  assert( l1(0,0.0) == 0 ); // test x-value smaller than first critical point
  assert( l1(0,1.0) == 0 ); // test x-value at first critical point
  assert( l1(0,2.0) == 1 ); // test x-value after first critical point
  assert( l1(0,1.001) >  0 ); // test x-value short after critical point
  assert( l1(0,4.0) == 2 ); // test x-value between critical points
  assert( l1(0,200) == 0 ); // test x-value larger than last critical point
  
  assert( l1(0,0.0) == l1[0](0) );
  assert( l1(0,1.0) == l1[0](1.0) );
  assert( l1(0,2.0) == l1[0](2.0) );  
  assert( l1(0,1.001) == l1[0](1.001) );
  assert( l1(0,4.0) == l1[0](4.0) );
  assert( l1(0,200) == l1[0](200) );
  
  PersistenceLandscape l3(l1); // test copy 

  l3 = l1 * 3.0;
  
  assert( l3.getX() == std::vector<std::vector<double>>(
    {{-inf, 1.0, 3.0, 3.5, 5.0, 6.5, 7.0, 9.0, 9.5, 10, inf}, 
     {-inf, 2.0, 3.5, 5.0, 6.5, 8.0, inf},
     {-inf, 3.0, 3.5, 4.0, 6.0, 6.5, 7.0, inf}} ) );
  assert( l3.getY() == std::vector<std::vector<double>>(
    {{ 0.0, 0.0, 6.0, 4.5, 9.0, 4.5, 6.0, 0.0, 1.5, 0, 0.0 }, 
     { 0.0, 0.0, 4.5, 0.0, 4.5, 0.0, 0.0 },
     { 0.0, 0.0, 1.5, 0.0, 0.0, 1.5, 0.0, 0.0 }} ) );
  
  l3 = 3.0 * l1;
  assert( l3.getX() == std::vector<std::vector<double>>(
    {{-inf, 1.0, 3.0, 3.5, 5.0, 6.5, 7.0, 9.0, 9.5, 10, inf}, 
     {-inf, 2.0, 3.5, 5.0, 6.5, 8.0, inf},
     {-inf, 3.0, 3.5, 4.0, 6.0, 6.5, 7.0, inf}} ) );
  
  assert( l3.getY() == std::vector<std::vector<double>>(
    {{ 0.0, 0.0, 6.0, 4.5, 9.0, 4.5, 6.0, 0.0, 1.5, 0, 0.0 }, 
     { 0.0, 0.0, 4.5, 0.0, 4.5, 0.0, 0.0 },
     { 0.0, 0.0, 1.5, 0.0, 0.0, 1.5, 0.0, 0.0 }} ) );
  
  l3 = l1 + l1;
  assert( l3.getX() == std::vector<std::vector<double>>(
    {{-inf, 1.0, 3.0, 3.5, 5.0, 6.5, 7.0, 9.0, 9.5, 10, inf}, 
     {-inf, 2.0, 3.5, 5.0, 6.5, 8.0, inf},
     {-inf, 3.0, 3.5, 4.0, 6.0, 6.5, 7.0, inf}} ) );
  assert( l3.getY() == std::vector<std::vector<double>>(
    {{ 0.0, 0.0, 4.0, 3.0, 6.0, 3.0, 4.0, 0.0, 1.0, 0, 0.0 }, 
     { 0.0, 0.0, 3.0, 0.0, 3.0, 0.0, 0.0 },
     { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0 }} ) );
  
  l3 = l1 + l2;
  assert( l3.getX() == std::vector<std::vector<double>>(
    {{-inf, 1.0, 3.0, 3.5, 4.5, 5.0, 6.5, 7.0, 9.0, 9.5, 10, inf}, 
     {-inf, 2.0, 2.5, 3.0, 3.5, 4.5, 5.0, 6.0, 6.5, 8.0, inf},
     {-inf, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 6.5, 7.0, inf}} ) );
  
  assert( l3.getY() == std::vector<std::vector<double>>(
    {{0, 0, 4, 4, 4, 5, 2, 2, 0, 0.5, 0, 0}, 
     {0, 0, 1, 1, 2, 2, 1, 1, 1.5, 0, 0},
     {0, 0, 0.5, 0, 0.5, 0, 0, 0.5, 0, 0}} ) );
  
  assert(l1 == l1);
  assert(l2 == l2);
  assert(l3 == l3);
  
  l2 += l1;
  assert(l2 == l3);
  
  l3 = l2 * 3.0;
  l2 *= 3.0;
  assert(l2 == l3);
  
  // Test output -------------------------------------------------------------
  l2 = PersistenceLandscape(data_2);
  l1.fileOutput("l1.dat");
  l2.fileOutput("l2.dat");
  
  // Test addition of landscapes with different number of layers.
  PersistenceLandscape l5(std::vector<std::vector<double>> {{-inf, 1.0, 3.0, 3.5, 4.5, 5.0, 6.5, 7.0, 9.0, 9.5, 10, inf}}, 
                          std::vector<std::vector<double>> {{0, 0, 4, 4, 4, 5, 2, 2, 0, 0.5, 0, 0}}
                         );
  l2 = l1 + l5;
  l2 = l5 + l1;
  assert( l1 + l5 == l5 + l1 );
  l2 = l1 + l0;
  assert ( l2 == l1 );
  assert ( l2 == (l2 + l0) );
  // Test Persistence Diagram IO ---------------------------------------------
  unsigned n = 10;
  PersistenceDiagram diag = createRandomPersistenceDiagram<double>(n);
  PersistenceLandscape l4(diag);
  
  // Test norm calculation
  std::deque<Interval> data_3 = { Interval(0.0, 1.0) };
  std::deque<Interval> data_4 = { Interval(1.0, 3.0) };
  
  l1 = PersistenceLandscape(data_3);
  l2 = PersistenceLandscape(data_4);
  auto res1 = l1.norm(1);
  auto res2 = l2.norm(1);
  assert( res1 == 0.25 );
  assert( res2 == 1.0 );

  l1 = PersistenceLandscape(data_1);
  l2 = PersistenceLandscape(data_2);
  assert( l1.norm(1) == 17.75 );
  //std::cout<< "l1 with p = 2" << l1.norm(2) << std::endl;
  //std::cout << "tu was " << std::endl;
  assert( l1.norm(2) < 5.37743 ); // dirty workaround for float comparison maybe do it better Later
  assert( l1.norm(2) > 5.37742 );
  //assert( l2.norm(1) = 
  //PersistenceDiagram<double> sphereDiag = createRandomSpherePersistenceDiagram(;
  // test Landscape Calculation based on a real data diagram from reddit graphs
  
  
  PersistenceDiagram D;
  
  D.add( 1, inf); 
  D.add( 1, inf); 
  D.add( 1, 58); 
  D.add( 1, 57); 
  D.add( 1, 58); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 58); 
  D.add( 1, 57); 
  D.add( 1, 12); 
  D.add( 1, 57); 
  D.add( 1, 11); 
  D.add( 1, 57); 
  D.add( 1, 58); 
  D.add( 1, 7 ); 
  D.add( 1, 7 ); 
  D.add( 1, 57); 
  D.add( 1, 57);
  D.add( 1, 10); 
  D.add( 1, 57); 
  D.add( 1, 9 ); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 10); 
  D.add( 1, 57); 
  D.add( 1, 8 ); 
  D.add( 1, 57);
  D.add( 1, 59); 
  D.add( 1, 58); 
  D.add( 1, 57); 
  D.add( 1, 57);
  D.add( 1, 7 ); 
  D.add( 1, 13); 
  D.add( 1, 7 ); 
  D.add( 1, 10); 
  D.add( 1, 57); 
  D.add( 1, 11); 
  D.add( 1, 5 ); 
  D.add( 1, 57); 
  D.add( 1, 7 ); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 58); 
  D.add( 1, 57); 
  D.add( 1, 8 ); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 10); 
  D.add( 1, 8 ); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 8 ); 
  D.add( 1, 5 ); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 57); 
  D.add( 1, 5 ); 
  D.add( 1, 8 ); 
  D.add( 1, 57); 

  
  PersistenceLandscape landscape_from_diag(D,0,115);
  landscape_from_diag.fileOutput("/tmp/test_diagram_with_equal_points.dat");
  PersistenceDiagram D1;
  
  D1.add(7 , inf); 
  D1.add(8 , inf); 
  D1.add(8 , inf); 
  D1.add(12, inf); 
  D1.add(14, inf); 
  D1.add(19, inf); 
  D1.add(75, inf); 
  D1.add(76, inf); 
  D1.add(78, inf); 
  D1.add(82, inf);
  
  PersistenceLandscape landscape_from_diag1(D1,0,115);

  // test absolute funktion --------------------------------------------------

  std::vector<std::vector<double>> X = {{-inf, 1.0, 2.0, 4.0, 5.0, 6.0, 8.0, 9.0, inf}, {-inf, 1.5, 2.5, 3.5, 4.5, inf}};
  std::vector<std::vector<double>> Y = {{ 0.0, 0.0, 1.0,-1.0, 0.0,-1.0, 1.0, 0.0, 0.0}, { 0.0, 0.0,-1.0,-1.0, 0.0, 0.0 }};

  PersistenceLandscape l6 = PersistenceLandscape(X,Y);
  PersistenceLandscape l7 = aleph::abs(l6);

  assert( l7.getX() == std::vector<std::vector<double>>({{-inf, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, inf}, {-inf, 1.5, 2.5, 3.5, 4.5, inf}}) );
  assert( l7.getY() == std::vector<std::vector<double>>({{ 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0}, { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0}}) );
  return;
}

int main ()
{
  testPersistenceLandscape();
}


/* 
  for debuggin outputs of landscapes the following code piece can be used:
  
  std::cout << "l3: " << std::endl;
  for (auto layer : l3.getX())
  {
    std::cout << "{ ";
    for (auto el : layer)
    {
      std::cout << el << " ";
    }
    std::cout << "}" << std::endl;
  }
*/
