#include <tests/Base.hh>
#include <aleph/persistenceDiagrams/PersistenceLandscape.hh>

#include <cmath>
#include <deque>
#include <limits>


void testPersistenceLandscape()
{
  
  double inf = std::numeric_limits<double>::infinity();
  //using namespace global;
  
  std::deque<pl::Interval<>> data_1 = {pl::Interval<>(1. , 5.), pl::Interval<>(2. , 8.), pl::Interval<>(3. , 4.), pl::Interval<>(5. , 9.), pl::Interval<>(6. , 7.), pl::Interval<>(9.,10.) };
  std::deque<pl::Interval<>> data_2 = {pl::Interval<>(1. , 6.), pl::Interval<>(2. , 3.), pl::Interval<>(3. , 7.), pl::Interval<>(4. , 5.) };

  pl::PersistenceLandscape<> l1(data_1);
  pl::PersistenceLandscape<> l2(data_2);
  
  //std::cout << "l1: " << std::endl;
  /*for (auto layer : l1.getX())
  {
    std::cout << "{ ";
    for (auto el : layer)
    {
      std::cout << el << " ";
    }
    std::cout << "}" << std::endl;
  }*/
  assert( l1.getX() == std::vector<std::vector<double>>({{-inf, 1.0, 3.0, 3.5, 5.0, 6.5, 7.0, 9.0, 9.5, 10, inf}, 
                                                         {-inf, 2.0, 3.5, 5.0, 6.5, 8.0, inf},
                                                         {-inf, 3.0, 3.5, 4.0, 6.0, 6.5, 7.0, inf}} ) );
  assert( l1.getY() == std::vector<std::vector<double>>({{ 0.0, 0.0, 2.0, 1.5, 3.0, 1.5, 2.0, 0.0, 0.5, 0, 0.0 }, 
                                                         { 0.0, 0.0, 1.5, 0.0, 1.5, 0.0, 0.0 },
                                                         { 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0 }} ) );
  

  assert( l1(0,0.0) == 0); // test x-value smaller than first critical point
  assert( l1(0,1.0) == 0); // test x-value at first critical point
  assert( l1(0,2.0) == 1); // test x-value after first critical point
  assert( l1(0,1.001) >  0); // test x-value short after critical point
  assert( l1(0,4.0) == 2); // test x-value between critical points
  assert( l1(0,200) == 0); // test x-value larger than last critical point
  
  assert( l1(0,0.0) == l1[0](0) );
  assert( l1(0,1.0) == l1[0](1.0));
  assert( l1(0,2.0) == l1[0](2.0));  
  assert( l1(0,1.001) == l1[0](1.001));
  assert( l1(0,4.0) == l1[0](4.0));
  assert( l1(0,200) == l1[0](200));
  
  // std::cout<<" check 1" << std::endl;
  
  pl::PersistenceLandscape<> l3(l1);
  // std::cout<<" check 4" << std::endl;

  l3 = l1 * 3.0;
  // std::cout << "l3: " << std::endl;
  /*for (auto layer : l3.getX())
  {
    std::cout << "{ ";
    for (auto el : layer)
    {
      std::cout << el << " ";
    }
    std::cout << "}" << std::endl;
  }*/
  assert( l3.getX() == std::vector<std::vector<double>>({{-inf, 1.0, 3.0, 3.5, 5.0, 6.5, 7.0, 9.0, 9.5, 10, inf}, 
                                                         {-inf, 2.0, 3.5, 5.0, 6.5, 8.0, inf},
                                                         {-inf, 3.0, 3.5, 4.0, 6.0, 6.5, 7.0, inf}} ) );
  assert( l3.getY() == std::vector<std::vector<double>>({{ 0.0, 0.0, 6.0, 4.5, 9.0, 4.5, 6.0, 0.0, 1.5, 0, 0.0 }, 
                                                         { 0.0, 0.0, 4.5, 0.0, 4.5, 0.0, 0.0 },
                                                         { 0.0, 0.0, 1.5, 0.0, 0.0, 1.5, 0.0, 0.0 }} ) );
  l3 = 3.0 * l1;
  
  
  // TODO: ergebnisse ausrechnen und korrekte asserts einfügen 
  l1 + l2;
  l1 * 3.0;
  l2 += l1;
  l2 *= 3.0;
  
  
  
  /*
  {-inf, 0} {1, 0} {3, 2} {3.5, 1.5} {5, 3} {6.5, 1.5} {7, 2} {9, 0} {9.5, 0.5} {10, 0} {inf, 0} 
  {-inf, 0} {2, 0} {3.5, 1.5} {5, 0} {6.5, 1.5} {8, 0} {inf, 0} 
  {-inf, 0} {3, 0} {3.5, 0.5} {6, 0} {4, 0} {6.5, 0.5} {7, 0} {inf, 0}
  
  {-inf, 0} {1, 0} {3.5, 2.5} {4.5, 1.5} {5, 2} {7, 0} {inf, 0} 
  {-inf, 0} {2, 0} {2.5, 0.5} {3, 0} {4.5, 1.5} {6, 0} {inf, 0} 
  {-inf, 0} {4, 0} {4.5, 0.5} {5, 0} {inf, 0}
  */
  return;
}

int main ()
{
  testPersistenceLandscape();
}
