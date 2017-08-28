#ifndef PERSISTNCELANDSCAPE_H_
#define PERSISTNCELANDSCAPE_H_

#include <iostream>
#include <algorithm>
#include <functional>
#include <limits>
#include <cmath>
#include <vector>
#include <deque>
#include <boost/icl/continuous_interval.hpp>
#include <aleph/utilities/ContainerOperators.hh>


// TODO: skalierungsinvarianz einbauen

namespace pl
{
// typedefs
template<typename T = double, typename R = double>
using CritPoint = std::pair<T,R>;

template<typename T = double>
using Interval = boost::icl::continuous_interval<T>; 

// Comparing structs
struct doubleEquality 
{
  double precision = 1e-08;
  inline bool operator() (double a, double b)
  {
    return (std::abs(a - b ) < precision);
  }
};
  
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
  
  class LandscapeLayer
  {
  private:
    PersistenceLandscape& outer;
    size_t k;
    //std::vector<T>& X_k;
    //std::vector<R>& Y_k;
  public:
    LandscapeLayer(PersistenceLandscape& outer, size_t k);
    R operator() (T t);
  };
  
  std::vector<LandscapeLayer> layerV; // stores founctors for the single layers
  
  const static T Infinity;// = std::numeric_limits<T>::infinity();
  const static T minusInfinity;// = -std::numeric_limits<T>::infinity();
  
  static PersistenceLandscape<T,R,S> linComb (std::vector<PersistenceLandscape<T,R,S>> landscapeV_, std::vector<S> scalarV)
  {
    size_t layerN = landscapeV_[0].layer();
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
    return PersistenceLandscape<T,R,S>(X_, Y_);
  }
  static std::pair<std::vector<T>, std::vector<R> > linComb_singleLayer ( const std::vector<PersistenceLandscape<T,R,S>>& landscapeV, const std::vector<S>& a,const size_t k)
  {
    size_t layerN = landscapeV[0].layer();
    size_t N = landscapeV.size();
    // merge the vectors evtl. with std::accumulate TODO: capsule this in a structure
    std::vector<T> X;
    for (auto landscape : landscapeV)
    {
      assert (landscape.layer() == layerN);

      std::vector<T> tmp = {};
      std::vector<T> X_i = landscape.getX()[k];
      std::merge(X.begin(), X.end(), X_i.begin(), X_i.end(), std::back_inserter(tmp) );
      X = tmp;
    }
    // erase duplicates
    std::unique(X.begin(), X.end(), [](const double& a, const double& b) { return a == b; } );
    // calculate stuff
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
    std::vector<double> Y(X.size(),0); // = std::accumulate(Y_.begin(), Y_.end(), std::vector<double>(X.size(), 0), [a] (std::vector<double> Y1, std::vector<double> Y2) { return Y1 + a[j]  * Y2
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
public:
  // constructors:
  PersistenceLandscape (std::deque<Interval<T>> A); // Constructor by birth-death pairs
  PersistenceLandscape (std::vector<std::vector<T>>, std::vector<std::vector<R>>); // Constructor by already claculated critical points
  PersistenceLandscape (const PersistenceLandscape& other);
  
  PersistenceLandscape& operator= (const PersistenceLandscape& other);
  ReturnType operator() (size_t k, ArgumentType t) const;
  LandscapeLayer operator[] (size_t k);
  
  // calculation operators
  PersistenceLandscape& operator+= (const PersistenceLandscape& other);
  PersistenceLandscape operator+ (const PersistenceLandscape& other) const;
  PersistenceLandscape& operator*= (const ScalarType& scalar);
  PersistenceLandscape operator* (const ScalarType& scalar) const;
  
  // PersistenceLandscape& operator/= (const ScalarType& scalar); // erst wenn es einen anwendungsfall gibt
  // PersistenceLandscape operator/ (const ScalarType& scalar) const;
  
  // non operator functionalities
  ScalarType norm(size_t p) const; // schoener als eigene Klasse mit template p und argument lambda?
  ReturnType integral (size_t k, size_t p) const;
  ReturnType integral (size_t p) const;
  ReturnType integral (T lower, T upper) const;
  std::vector<std::vector<T> > getX() const; // maybe it should return a reference?
  std::vector<std::vector<T> > getY() const; // s.o.
  size_t layer() const;
};


template <typename T, typename R, typename S> 
constexpr T PersistenceLandscape<T,R,S>::Infinity = std::numeric_limits<T>::infinity();

template <typename T, typename R, typename S> 
constexpr T PersistenceLandscape<T,R,S>::minusInfinity = -std::numeric_limits<T>::infinity();

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
PersistenceLandscape<T,R,S>::PersistenceLandscape(std::deque<Interval<T>> A) 
{
  size_t k = 0;
  std::vector<std::vector<CritPoint<T,R>> > L;

  std::sort(A.begin(), A.end(), intervalCompare());

  while ( !(A.empty()) )
  {
    // for (auto it : A ) { std::cout << it << " "; } std::cout << std::endl;
    L.push_back(std::vector<CritPoint<T,R>>());
    assert(L.size() == k+1);
    // std::cout<< "check 1" << std::endl;
    Interval<T> current = A.front();
    A.pop_front();
    double b = lower(current); // aufpassen ohne referenz funktioniert es evtl nicht richtig
    double d = upper(current);
    auto p = A.begin();
    
    //std::cout << *p << std::endl;
    //std::cout << "k is: " << k << " L[k]: " << std::endl;
    //std::cout << L.size() << std::endl;
    //std::cout <<"check 1.1"<< std::endl;
    L[k].push_back(CritPoint<T,R>(minusInfinity,0.) );
    //std::cout<< "check 1.2" << std::endl;
    L[k].push_back(CritPoint<T,R>(b,0.) );
    //std::cout<< "check 1.3" << std::endl;
    
    L[k].push_back(CritPoint<T,R>( (b + d) / 2. , (d - b) / 2. ));
    
    while ( L[k].back() != CritPoint<T,R>(Infinity,0.) ) // anders herum??
    {
      //std::cout<< "check 2" << std::endl;
      //for (auto it : A ) { std::cout << it << " "; } std::cout << std::endl;

      p = std::find_if(p, A.end(), [d] (auto interval)->bool {return interval.upper() > d;} ); //wsl hier segfault
      //std::cout << *p << std::endl;
      if ( p == A.end() )
      {
        L[k].push_back(CritPoint<T,R>(d,0.));
        L[k].push_back(CritPoint<T,R>(Infinity,0.));
      }
      else
      {
        double b_ = p->lower(); //std::cout<< "b_: " << b_<< " "; std::cout<< "b: " << b<< std::endl;
        double d_ = p->upper(); //std::cout<< "d_: " << d_<< " "; std::cout<< "d: " << d<< std::endl;
        auto h = p;
        //std::cout<< "h: " << *h << std::endl;
        p = A.erase(h);// evtl: erase p 
        //for (auto it : A ) { std::cout << it << " "; } std::cout << std::endl;
        //std::cout << *p << std::endl;
        if (b_ == d)
        {
          //std::cout<< " b_ equals d" << std::endl;
          L[k].push_back(CritPoint<T,R>(b_,0.) );
        }
        else if (b_ > d)
        {
          //std::cout<< " b_ > d" << std::endl;
          L[k].push_back(CritPoint<T,R>(d ,0.) );
          L[k].push_back(CritPoint<T,R>(b_,0.) );
        }
        else
        {
          //std::cout<< " b_ < d" << std::endl; std::cout<< "b_: " << b_<< " "; std::cout<< "d: " << d<< std::endl;
          L[k].push_back(CritPoint<T,R>( (b_+d) / 2. , (d-b_) / 2. ));
          //std::cout<< " p: " << *p<< std::endl;
          //std::cout<< " b_ < d" << std::endl; std::cout<< "b_: " << b_<< " "; std::cout<< "d: " << d<< std::endl;
          //for (auto it : A ) { std::cout << it << " "; } std::cout <<"check" << std::endl;
          //std::cout << "p: " <<*p << std::endl;
          A.insert(p, Interval<T>(b_, d) );
          //for (auto it : A ) { std::cout << it << " "; } std::cout <<"check" << std::endl;
          //std::cout<< " b_ < d" << std::endl; std::cout<< "b_: " << b_<< " "; std::cout<< "d: " << d<< std::endl;
          p++;
        }
        L[k].push_back(CritPoint<T,R>( (b_+d_) / 2. , (d_-b_) / 2. ));
        b = b_; // kann nicht wirklich stimmen auf referenzen aufpassen
        d = d_; // s.o.
      }
    }
    ++k;
  }
  auto result = aleph::utilities::cast_pair(L);
  this->X = result.first;
  this->Y = result.second;
}

template <typename T, typename R, typename S>
PersistenceLandscape<T,R,S>& PersistenceLandscape<T,R,S>::operator= (const PersistenceLandscape& other)
{
  this->X = other.getX();
  this->Y = other.getY();
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
  double m = (Y[k][index] - Y[k][index - 1]) / (X[k][index] - X[k][index - 1]); // TODO: catch 0 division
  return (Y[k][index] - (m * (X[k][index] - t) ));
}

template <typename T, typename R, typename S>
typename PersistenceLandscape<T,R,S>::LandscapeLayer PersistenceLandscape<T,R,S>::operator[] (size_t k ) // maybe return reference here?
{
  return LandscapeLayer(*this, k); // is it better or worse to construct the matching Layer object here?
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
R PersistenceLandscape<T,R,S>::integral (size_t k,size_t p) const
{
  R surface;
  for (size_t i = 2; i < X.size(); i++)
  {
    // steigung
    R a = (Y[k][i]-Y[k][i-1]) / (X[k][i]-X[k][i]); 
    // y-achsenabschnitt
    R b = Y[i] - (a * X[i]);
    if ( a != 0) //besser Y[k][i] == Y[k][i-1] ?
    {
      surface += 1/(a * p + 1.0) * ( pow( (a*X[i] + b) , p + 1 ) - pow( (a*X[i-1] + b), p + 1) ); 
      //TODO: klären ob es nicht ein integral für positive Komponenten und eines für negative geben müsste.
    }
    else 
    {
      surface += (X[i] - X[i-1]) * pow( Y[i], p ); 
    } // explanation insert
  }
  return surface;
}

template <typename T, typename R, typename S>
R PersistenceLandscape<T,R,S>::integral (size_t p) const
{
  R surface = R(0); // does this work??
  for (size_t k = 0; k < X.size(); k++)
  {
    surface += integral(k, p);
  }
  return surface;
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
R __attribute__((optimize("O0"))) PersistenceLandscape<T,R,S>::LandscapeLayer::operator() ( T t)
{
  auto xIter = std::lower_bound(outer.X[k].begin(), outer.X[k].end(), t);
  if (*xIter == Infinity) {return 0;} // t > last crit point
  else if( *xIter == outer.X[k][0]) { return 0; }
  else if( *xIter == outer.X[k][1]) { return 0;} 
  size_t index = xIter - outer.X[k].begin();
  // linear interpolation
  double m = (outer.Y[k][index] - outer.Y[k][index - 1]) / (outer.X[k][index] - outer.X[k][index - 1]); 
  // TODO: catch 0 division, doesnt happen does it?
  return (outer.Y[k][index] - (m * (outer.X[k][index] - t) )); 
}


/* 
 * TODO: 
 * landscape layer 
 *  - operator+ (landscape layer)
 *  - operator- (landscape layer)
 *  - operator+ (landscape)
 * persistence landscape
 *  - operator- (landscape)
 *  - operator+ (landscape layer)
 */

//template <typename T, typename R, typename S>

}

#endif // PERSISTNCE_LANDSCAPE_IMP_H_
