#ifndef PERSISTENCELANDSCAPE_H_
#define PERSISTENCELANDSCAPE_H_

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <limits>
#include <cmath>
#include <vector>
#include <deque>

#include <stdexcept>

#include <boost/icl/continuous_interval.hpp>

#include <aleph/utilities/ContainerOperators.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>


// TODO: skalierungsinvarianz einbauen

namespace aleph
{
// typedefs
template<typename T = double, typename R = double>
using CritPoint = std::pair<T,R>;

template<typename T = double>
using Interval = std::pair<T,T>;

template<typename T>
T upper(Interval<T> interval)
{
  return interval.second;
}

template<typename T>
T lower(Interval<T> interval)
{
  return interval.first;
}
// Comparing structs
bool equal (double a, double b, double precision = 1e-04)
{
  return (std::abs( a - b ) < precision);
}
  
struct intervalCompare
{
  inline bool operator() (const Interval<>& a, const Interval<>& b)
  {
    // its important to evaluate this first, because of the special nature of float equality
    // but this solution is not very elegant, because we have to evaluate the equality term 
    // for all cases, while in most of them it will be unneccessary
    // TODO: Ask for a better solution
    if ( ((lower(a) == lower(b))) && (upper(a) > upper(b)) ) 
    {
      return true;
    }
    else if (lower(a) < lower(b))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
};


template <typename T = double ,typename R = double ,typename S = double>
class PersistenceLandscape 
{
  using ArgumentType = T; // Real numbers
  using ReturnType = R; // Real numbers with infinity
  using ScalarType = S; // Real numbers
  
private:
  //Interval<T> intervals;
  //std::vector<Functor> functions;
  std::vector<std::vector<T> > X; // Critical x-values
  std::vector<std::vector<R> > Y; // Critical y-values
  T intervalMax;
  T intervalMin;
  
  class LandscapeLayer
  {
  private:
    PersistenceLandscape& outer;
    size_t k;
    //std::vector<T>& X_k;
    //std::vector<R>& Y_k;
  public:
    LandscapeLayer(PersistenceLandscape& outer, size_t k);
    T X() const;
    R Y() const;
    R operator() (T t);
  };
  
  //std::vector<LandscapeLayer> layerV; // stores founctors for the single layers
  
  const static T Infinity;// = std::numeric_limits<T>::infinity();
  const static T minusInfinity;// = -std::numeric_limits<T>::infinity();
  
  static PersistenceLandscape<T,R,S> linComb (std::vector<PersistenceLandscape<T,R,S>> landscapeV_, std::vector<S> scalarV);
  static std::pair<std::vector<T>, std::vector<R> > linComb_singleLayer ( const std::vector<PersistenceLandscape<T,R,S>>& landscapeV, const std::vector<S>& a,const size_t k);
  
public:
  // constructors:
  PersistenceLandscape ();
  PersistenceLandscape (std::deque<Interval<T>> A); // Constructor by birth-death pairs  
  PersistenceLandscape (const aleph::PersistenceDiagram<T>& diag, T intervalMin, T intervalMax);
  PersistenceLandscape (const aleph::PersistenceDiagram<T>& diag);
  PersistenceLandscape (std::vector<std::vector<T>>, std::vector<std::vector<R>>); // Constructor by already claculated critical points
  PersistenceLandscape (const PersistenceLandscape& other);
  
  PersistenceLandscape& operator= (const PersistenceLandscape& other);
  ReturnType operator() (size_t k, ArgumentType t) const;
  LandscapeLayer operator[] (size_t k);
  
  bool operator== (const PersistenceLandscape& other) const;
  
  // calculation operators
  PersistenceLandscape& operator+= (const PersistenceLandscape& other);
  PersistenceLandscape operator+ (const PersistenceLandscape& other) const;
  PersistenceLandscape& operator*= (const ScalarType& scalar);
  PersistenceLandscape operator* (const ScalarType& scalar) const;
  
  
  // non operator functionalities
  R norm(double p) const; // schoener als eigene Klasse mit template p und argument lambda?
  ReturnType integral (size_t k, double p) const;
  //ReturnType integral (size_t p) const;
  ReturnType integral (T lower, T upper) const;
  std::vector<std::vector<T> > getX() const; // maybe it should return a reference?
  std::vector<std::vector<T> > getY() const; // s.o.
  size_t layer() const;
  size_t size() const;
  //void restrictIntervals
  void fileOutput(std::string filename) const;
  T getIntervalMin() const;
  T getIntervalMax() const;
};


template <typename T, typename R, typename S> 
constexpr T PersistenceLandscape<T,R,S>::Infinity = std::numeric_limits<T>::infinity();

template <typename T, typename R, typename S> 
constexpr T PersistenceLandscape<T,R,S>::minusInfinity = -std::numeric_limits<T>::infinity();

template <typename T, typename R, typename S> 
PersistenceLandscape<T,R,S> PersistenceLandscape<T,R,S>::linComb (std::vector<PersistenceLandscape<T,R,S>> landscapeV_, std::vector<S> scalarV)
{
    // dirty hack for implementing addition of two landscapes of different size
    // TODO implement iterators and + as well as += better
    size_t layerN = 0;
    std::vector<std::vector<T>> X_tail;
    std::vector<std::vector<T>> Y_tail;
    
    if ( landscapeV_.size() == 2 )
    {
      int diff = static_cast<int>(landscapeV_[0].layer()) - static_cast<int>(landscapeV_[1].layer());
      if (diff > 0) 
      {
        layerN = landscapeV_[1].layer();
        auto tmpX = landscapeV_[0].getX();
        auto tmpY = landscapeV_[0].getY();
        X_tail = std::vector<std::vector<T>>( (tmpX.end() - diff) , tmpX.end() );
        Y_tail = std::vector<std::vector<T>>( (tmpY.end() - diff) , tmpY.end() );
      }
      else if (diff < 0)
      {
        layerN = landscapeV_[0].layer();
        auto tmpX = landscapeV_[1].getX();
        auto tmpY = landscapeV_[1].getY();
        X_tail = std::vector<std::vector<T>>( (tmpX.end() - std::abs(diff)) , tmpX.end() );
        Y_tail = std::vector<std::vector<T>>( (tmpY.end() - std::abs(diff)) , tmpY.end() );
      }
      else 
      {
        layerN = landscapeV_[0].layer();
      }
    }
    else
    {
      // bad but i have to assume here that this function only gets called by operators who add two vectors
      // this should be implemented in a better way at a later time.
      layerN = landscapeV_[0].layer();
      if (landscapeV_.size() > 2) assert(false);
    }
    
    //bad coding style
    //auto minLandscape = std::min_element( landscapeV_.begin(), landscapeV.end(), 
    //                                  [] (const PersistenceLandscape& a, const PersistenceLandscape& a) 
    //                                   { return a.layer() < b.layer(); }
    //                  );
    //size_t layerN = minLandscape.layer();
    std::vector<std::vector<T>> X_ = {};
    std::vector<std::vector<R>> Y_ = {};
    // std::cout<< "check 2" << std::endl;
    for(size_t k = 0; k < layerN; k++)
    {
      //std::vector<std::pair< std::vector<T>, std::vector<R> >> inputV;
        //assert (landscapeV_[].layer() == layerN);
        //inputV.push_back(std::pair< std::vector<T>, std::vector<R> >(landscape.getX()[k], landscape.getY()[k]) );
        auto result = linComb_singleLayer(landscapeV_, scalarV, k);
        X_.push_back(result.first);
        Y_.push_back(result.second);
    }
    X_.insert( X_.end(), X_tail.begin(), X_tail.end() );
    Y_.insert( Y_.end(), Y_tail.begin(), Y_tail.end() );
    return PersistenceLandscape<T,R,S>(X_, Y_);
}

template <typename T, typename R, typename S>
std::pair<std::vector<T>, std::vector<R> > PersistenceLandscape<T,R,S>::linComb_singleLayer ( const std::vector<PersistenceLandscape<T,R,S>>& landscapeV, const std::vector<S>& a,const size_t k )
{
    //size_t layerN = landscapeV[0].layer();
    size_t N = landscapeV.size();
    // merge the vectors evtl. with std::accumulate TODO: capsule this in a structure
    std::vector<T> X;
    for (auto landscape : landscapeV)
    {
      //assert (landscape.layer() == layerN);

      std::vector<T> tmp = {};
      std::vector<T> X_i = landscape.getX()[k];
      std::merge(X.begin(), X.end(), X_i.begin(), X_i.end(), std::back_inserter(tmp));
      X = tmp;
    //  for (auto x : X)
    //  { std::cout << x << " ";}
    //  std::cout << std::endl;
    }
    // erase duplicates
    auto last = std::unique(X.begin(), X.end() /*, [](const double& a, const double& b) { return a == b; }*/ );
    X.erase(last, X.end());
    // calculate stuff
  //  std::cout << "erased: " << std::endl;
  //  for (auto x : X)
  //    { std::cout << x << " ";}
  //    std::cout << std::endl;
  //  std::cout << std::endl;
    std::vector<std::vector<R> > Y_;
    for (size_t j = 0; j < N; j++)
    {
      Y_.push_back(std::vector<R>() );
      for (auto x_i : X )
      {
        Y_[j].push_back( landscapeV[j](k, x_i) ); 
      }
    }
    // encapsule this in tools. accumulate or something
    std::vector<R> Y(X.size(),0); // = std::accumulate(Y_.begin(), Y_.end(), std::vector<double>(X.size(), 0), [a] (std::vector<double> Y1, std::vector<double> Y2) { return Y1 + a[j]  * Y2
    for (size_t j = 0; j < N; j++)
    { 
      { 
        using aleph::utilities::operator*;
        using aleph::utilities::operator+=;
        Y += (a[j] * Y_[j]);
      }
    }
    
    return std::make_pair(X,Y);
}  

template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S>::PersistenceLandscape()
{
  this->X = {{minusInfinity, Infinity}};
  this->Y = {{0, 0}};
}


template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S>::PersistenceLandscape(std::vector<std::vector<T>> X_, std::vector<std::vector<R>> Y_) : X(X_), Y(Y_)
{}

template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S>::PersistenceLandscape(const PersistenceLandscape& other)
{
  this->X = other.getX();
  this->Y = other.getY();
}


template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S>::PersistenceLandscape (std::deque<Interval<T>> A) 
{
  //R MAX = 2000;
  size_t k = 0;
  std::vector<std::vector<CritPoint<T,R>> > L;

  std::sort(A.begin(), A.end(), intervalCompare());
  //std::unique(A.begin(), A.end()); 
  while ( !(A.empty()) )
  {
    // for (auto it : A ) { std::cout << it << " "; } std::cout << std::endl;
    L.push_back(std::vector<CritPoint<T,R>>());
    assert(L.size() == k+1);
    // std::cout<< "check 1" << std::endl;
    Interval<T> current = A.front();
    A.pop_front();
    T b = lower(current); 
    T d = upper(current);
    auto p = A.begin();
    if( b == minusInfinity && d == Infinity )
    {
      throw;
      L[k].push_back(CritPoint<T,R>(minusInfinity, Infinity));
      L[k].push_back(CritPoint<T,R>(Infinity, Infinity));
    }
    else
    {
      if (d == Infinity)
      {
        throw;
        L[k].push_back(CritPoint<T,R>(minusInfinity,0));
        L[k].push_back(CritPoint<T,R>(b,0));
        L[k].push_back(CritPoint<T,R>(Infinity,Infinity));
      }
      else
      {
        if(b == minusInfinity)
        {
          throw;
          L[k].push_back(CritPoint<T,R>(minusInfinity,Infinity));
        }
        else
        {
          L[k].push_back(CritPoint<T,R>(minusInfinity,0.) );
          L[k].push_back(CritPoint<T,R>(b,0.) );
          L[k].push_back(CritPoint<T,R>( (b + d) / 2. , (d - b) / 2. ));
        }
      }
    }
    
    //std::cout << *p << std::endl;
    //std::cout << "k is: " << k << " L[k]: " << std::endl;
    //std::cout << L.size() << std::endl;
    //std::cout <<"check 1.1"<< std::endl;
    //std::cout<< "check 1.2" << std::endl;
    //std::cout<< "check 1.3" << std::endl;
    
    
    while ( L[k].back() != (CritPoint<T,R>(Infinity,0.)) && L[k].back() != CritPoint<T,R>(Infinity,Infinity) ) // anders herum??
    {
      //std::cout<< "check 2" << std::endl;
      //for (auto it : A ) { std::cout << it << " "; } std::cout << std::endl;

      p = std::find_if(p, A.end(), [d] (auto interval)->bool {return upper(interval) > d;} ); //wsl hier segfault
      //std::cout << *p << std::endl;
      if ( p == A.end() )
      {
        L[k].push_back(CritPoint<T,R>(d,0.));
        L[k].push_back(CritPoint<T,R>(Infinity,0.));
      }
      else
      {
        T b_ = lower(*p); //std::cout<< "b_: " << b_<< " "; std::cout<< "b: " << b<< std::endl;
        T d_ = upper(*p); //std::cout<< "d_: " << d_<< " "; std::cout<< "d: " << d<< std::endl;
        auto h = p;
        //std::cout<< "h: " << *h << std::endl;
        p = A.erase(h);// evtl: erase p 
        //for (auto it : A ) { std::cout << it << " "; } std::cout << std::endl;
        //std::cout << *p << std::endl;
        if (b_ == d)
        {
          //std::cout<< " b_ equals d" << std::endl;
          L[k].push_back(CritPoint<T,R>(b_,0.) );
          if( (L[k].end()-2)->first == L[k].back().first)
          {
            throw std::logic_error("Double X value added");
          }
        }
        else if (b_ > d)
        {
          //std::cout<< " b_ > d" << std::endl;
          L[k].push_back(CritPoint<T,R>(d ,0.) );
          L[k].push_back(CritPoint<T,R>(b_,0.) );
          if( ((L[k].end()-2)->first) == (L[k].back().first))
          {
            throw std::logic_error("Double X value added");
          }
        }
        else
        {

          //std::cout<< " b_ < d" << std::endl; std::cout<< "b_: " << b_<< " "; std::cout<< "d: " << d<< std::endl;
          L[k].push_back(CritPoint<T,R>( (b_+d) / 2. , (d-b_) / 2. ));
          if( (L[k].end()-2)->first == L[k].back().first)
          {
            throw std::logic_error("Double X value added");
          }
          
          //std::cout<< " p: " << *p<< std::endl;
          //std::cout<< " b_ < d" << std::endl; std::cout<< "b_: " << b_<< " "; std::cout<< "d: " << d<< std::endl;
          //for (auto it : A ) { std::cout << it << " "; } std::cout <<"check" << std::endl;
          //std::cout << "p: " <<*p << std::endl;
          p = A.insert( std::upper_bound(p, A.end(), Interval<T>(b_, d), intervalCompare() ),
                        Interval<T>(b_, d) );
          //for (auto it : A ) { std::cout << it << " "; } std::cout <<"check" << std::endl;
          //std::cout<< " b_ < d" << std::endl; std::cout<< "b_: " << b_<< " "; std::cout<< "d: " << d<< std::endl;
          p++;
        }
        if( d_ == Infinity )
        {
          throw;
          L[k].push_back(CritPoint<T,R>(Infinity,Infinity));
        }
        else
        {
          L[k].push_back(CritPoint<T,R>( (b_+d_) / 2. , (d_-b_) / 2. ));
          if( (L[k].end()-2)->first == L[k].back().first)
          {
            throw std::logic_error("Double X value added");
          }
          b = b_;
          d = d_;
        }
      }
    }
    ++k;
  }
  // TODO: think about this, im not sure if this behaves correctly
  //std::unique(L.begin(), L.end());
  //L.erase(last);
  auto result = aleph::utilities::cast_pair(L);

  
  this->X = result.first;
  this->Y = result.second;
  
  // assert that nothing went wrong
  
  //std::for_each(X.begin(), X.end(), [] (auto l) { return std::for_each(l.begin(),l.end(), [] (auto b) {return b > 0;}); });
  
  for (size_t i = 0; i < X.size(); i++)
  {
    for (size_t j = (X[i].size()-1); j > 0; j--)
    {
      if ((X[i][j] == X[i][j-1]) )//&& (Y[i][j] == Y[i][j-1]))
      {
        throw std::logic_error( "double X value at (" 
                                + std::to_string(X[i][j])
                                + ","
                                + std::to_string(Y[i][j])
                                + ") with index "
                                + std::to_string(j)
                                + "in layer "
                                + std::to_string(i)
                              );
       // X[i].erase(X[i].begin() + j);
       // Y[i].erase(Y[i].begin() + j);
      }
      if (Y[i][j] < 0)
      {
        throw std::logic_error( "invalid Y value at ("
                                + std::to_string(X[i][j])
                                + ","
                                + std::to_string(Y[i][j])
                                + ") with index "
                                + std::to_string(j)
                                + "in layer "
                                + std::to_string(i)
                              );
      }
    }
  }
  
  //std::unique(X[k].begin(),X[k].end());
  //std::unique(Y[k].begin(),Y[k].end());
}

template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S>::PersistenceLandscape( const aleph::PersistenceDiagram<T>& diag) : PersistenceLandscape<T,R,S>(diag, std::numeric_limits<T>::min(), std::numeric_limits<T>::max()) 
{}


template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S>::PersistenceLandscape( const aleph::PersistenceDiagram<T>& diag, T intervalMin_, T intervalMax_)
{
  this->intervalMin = intervalMin_;
  this->intervalMax = intervalMax_;
  
  std::deque<Interval<T>> queue = {};
  for( auto point : diag )
  {
    auto lowerBound = std::sqrt(2.0) * point.x();
    auto upperBound = std::sqrt(2.0) * point.y();
    if (lowerBound < intervalMin) {
      lowerBound = intervalMin;
    }
    if (upperBound > intervalMax)
    {
      upperBound = intervalMax;
    }
    queue.push_back( Interval<T>(lowerBound, upperBound) );
  }
  *this = PersistenceLandscape(queue);
}

template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S>& PersistenceLandscape<T,R,S>::operator= (const PersistenceLandscape& other)
{
  this->X = other.getX();
  this->Y = other.getY();
  this->intervalMax = other.getIntervalMax();
  this->intervalMin = other.getIntervalMin();
  return *this;
}

template <typename T, typename R, typename S>
R PersistenceLandscape<T,R,S>::operator() (size_t k, T t) const
{
  auto xIter = std::lower_bound(X[k].begin(), X[k].end(), t);
  if (*xIter == Infinity) {return 0;}
  else if( *xIter <= X[k][1]) { return 0;}
  size_t index = xIter - X[k].begin();
  // linear interpolation
  R m = (Y[k][index] - Y[k][index - 1]) / (X[k][index] - X[k][index - 1]); // TODO: catch 0 division
  return (Y[k][index] - (m * (X[k][index] - t) ));
}

template <typename T, typename R, typename S>
typename PersistenceLandscape<T,R,S>::LandscapeLayer PersistenceLandscape<T,R,S>::operator[] (size_t k ) // maybe return reference here?
{
  return LandscapeLayer(*this, k); // is it better or worse to construct the matching Layer object here?
}

template <typename T, typename R, typename S>
bool PersistenceLandscape<T,R,S>::operator== (const PersistenceLandscape& other) const
{
  return (this->X == other.getX() && this->Y == other.getY() );
}

// unterstützt derzeit nur Landscapes auf den selben typen. TODO: unterschiedliche Typen ermöglichen.
template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S> PersistenceLandscape<T,R,S>::operator+ (const PersistenceLandscape& other) const
{
  std::vector<PersistenceLandscape> inputV = {*this, other};
  std::vector<S> scalarV = {1.0, 1.0};
  return PersistenceLandscape<T,R,S>::linComb(inputV, scalarV);
}

template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S>& PersistenceLandscape<T,R,S>::operator+= (const PersistenceLandscape& other)
{
  std::vector<PersistenceLandscape> inputV = {*this, other};
  std::vector<S> scalarV = {1.0, 1.0};
  *this = PersistenceLandscape<T,R,S>::linComb(inputV, scalarV);
  return *this;
}

template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S> PersistenceLandscape<T,R,S>::operator* (const S& scalar) const
{
  std::vector<PersistenceLandscape> inputV = {*this};
  std::vector<S> scalarV = {scalar};
  return PersistenceLandscape<T,R,S>::linComb(inputV, scalarV);
}

template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S>& PersistenceLandscape<T,R,S>::operator*= (const S& scalar)
{
  std::vector<PersistenceLandscape> inputV = {*this};
  std::vector<S> scalarV = {scalar};
  *this = PersistenceLandscape<T,R,S>::linComb(inputV, scalarV);
  return *this;
}

template <typename T, typename R, typename S>
R PersistenceLandscape<T,R,S>::integral (size_t k,double p) const
{
  //std::unique(X[k].begin(),X[k].end());
  //std::unique(Y[k].begin(),Y[k].end());
  assert (X[k].size() == Y[k].size());
  //std::cout << "layer k is: " << k << std::endl;
  //using aleph::utilities::operator*;
  //using aleph::utilities::operator+=;
  R surface(0);
  for (size_t i = 2; i < X[k].size()-1; i++)
  { //maybe look for problems with numerical stability??
  //  if (X[k][i] == X[k][i-1]) // extremly dirty fix why does this actually happen?????????? TODO
  //  {
  //    continue;
  //  }
    // steigung:
    R a = (Y[k][i]-Y[k][i-1]) / (X[k][i]-X[k][i-1]); 
    //std::cout<< "a: " << a << std::endl;
    // y-achsenabschnitt:
    R b = Y[k][i] - (a * X[k][i]);
    //std::cout << "b: " << b << std::endl;
    if ( a != 0) //besser Y[k][i] == Y[k][i-1] ?
    {
      auto res = ( 1/(a * (p + 1.0)) * ( pow( (a*X[k][i] + b) , p + 1 ) - pow( (a*X[k][i-1] + b), p + 1) ) ); 
      //std::cout<< "surface: " << res << std::endl;
      surface += res;
    }
    else 
    {
      auto res = (X[k][i] - X[k][i-1]) * pow( Y[k][i], p ); 
      //std::cout<< "surface: " << res << std::endl;
      surface += res;
    } // explanation insert
  }
  if (!std::isfinite(surface))
  {
    throw std::domain_error("integration result is infinity");
  }
  else if (surface < 0)
  {
    throw std::domain_error("integration result is negative");
  }
  //std::cout<< "surface for k = " << k << ": " << surface << std::endl;
  return surface;
}

template <typename T, typename R, typename S>
R PersistenceLandscape<T,R,S>::norm (double p) const
{
  //std::cout << "p: " << p << std::endl;
  if ( p <= 0 ) 
  {
    throw std::invalid_argument("p is not positive"); 
  }
  R surface = R(0);
  for (size_t k = 0; k < X.size(); k++)
  {
    surface += integral(k, p); // k layernumber, p norm index     
  }
  auto result = pow( surface, (1/p) );
  if (!std::isfinite(result))
  {
    throw std::domain_error("norm result is infinity");
  }
  return result;
}

template <typename T, typename R, typename S>
R PersistenceLandscape<T,R,S>::integral (T lower, T upper) const
{
  throw;
}
template <typename T, typename R, typename S>
std::vector<std::vector<T> > PersistenceLandscape<T,R,S>::getX() const
{
  return this->X;
}

template <typename T, typename R, typename S>
std::vector<std::vector<T> > PersistenceLandscape<T,R,S>::getY() const
{
  return this->Y;
}

template <typename T, typename R, typename S>
T PersistenceLandscape<T,R,S>::getIntervalMin() const
{
  return this->intervalMin;
}

template <typename T, typename R, typename S>
T PersistenceLandscape<T,R,S>::getIntervalMax() const
{
  return this->intervalMax;
}

template <typename T, typename R, typename S>
size_t PersistenceLandscape<T,R,S>::size() const
{
 // size_t result = 0;
 // for (size_t l = 0; l < this->layer(); l++)
 // {
 //   result += getX()
 // }
  return std::accumulate(X.begin(), X.end(), size_t(0), [] (auto a, auto b) { return a + b.size(); });
}

template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S> operator* (const S& scalar, const PersistenceLandscape<T,R,S>& landscape)
{
  return landscape * scalar;
}



template <typename T, typename R, typename S>
size_t PersistenceLandscape<T,R,S>::layer() const
{
  //std::cout << "check 3" << std::endl;
  assert(this->X.size() == this->Y.size());
  return this->X.size();
}

template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S>::LandscapeLayer::LandscapeLayer(PersistenceLandscape& outer_, size_t k_) : outer(outer_), k(k_)
{}

template <typename T, typename R, typename S>
R /*__attribute__((optimize("O0")))*/ PersistenceLandscape<T,R,S>::LandscapeLayer::operator() ( T t)
{
  // falls k > X.size() return trivial layer 0
  auto xIter = std::lower_bound(outer.X[k].begin(), outer.X[k].end(), t);
  if (*xIter == Infinity) {return 0;} // t > last crit point
  else if( *xIter == outer.X[k][0]) { return 0; }
  else if( *xIter == outer.X[k][1]) { return 0; } 
  size_t index = xIter - outer.X[k].begin();
  // linear interpolation
  double m = (outer.Y[k][index] - outer.Y[k][index - 1]) / (outer.X[k][index] - outer.X[k][index - 1]); 
  // TODO: catch 0 division, doesnt happen does it?
  return (static_cast<R>(outer.Y[k][index] - (m * (outer.X[k][index] - t) ))); 
}

template <typename T, typename R, typename S>
void __attribute__((optimize("O0"))) PersistenceLandscape<T,R,S>::fileOutput (std::string filename) const
{
  std::ofstream outputFile;
  outputFile.open(filename);
  outputFile << *this;
}

template <typename T, typename R, typename S>
std::ofstream& __attribute__((optimize("O0"))) operator<< ( std::ofstream& os, const PersistenceLandscape<T, R, S>& landscape)
{
  auto X = landscape.getX();
  auto Y = landscape.getY();
  
  os << "# persistence landscape" << std::endl;
  os << "# X Y" << std::endl;
  for ( size_t k = 0; k < X.size(); k++ )
  {
    os << "\"layer " << k << ":\"" << std::endl;
    for ( size_t i = 0; i < X[k].size(); i++ )
    {
      os << X[k][i] << " " << Y[k][i] << std::endl;
    }
    os << std::endl;
    os << std::endl;
  }
  
  return os;
}

} // namespace aleph

#endif // PERSISTENCE_LANDSCAPE_IMP_H_
     
