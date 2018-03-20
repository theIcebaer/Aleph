/* ---------------------------------------------------------------------------
 * 
 * This exapmple demonstrates the creation of a persistence landscapes from
 * several geometrical objects e.g. sphere, torus, ...
 * 
 * For this purpose we will use the functions described in the example
 * "create_persistence_diagrams.cc" to build random persistence diagrams from
 * geometrical samples and calculate the persistence landscapes from those
 * diagrams.
 * 
 * ---------------------------------------------------------------------------
 */

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/SphereSampling.hh>
#include <aleph/geometry/TorusSampling.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/PersistenceLandscape.hh>
#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <stdexcept>

#include <cmath>

#include <getopt.h>

// We first have to specify the data type of the persistence diagram,
// i.e. the type that is used by its individual points.
using DataType              = double;
using Distance              = aleph::distances::Euclidean<DataType>;
using PointCloud            = aleph::containers::PointCloud<DataType>;
using PersistenceDiagram    = aleph::PersistenceDiagram<DataType>;
using PersistenceLandscape  = aleph::PersistenceLandscape<DataType>;

auto infty = std::numeric_limits<DataType>::infinity();

namespace ch = std::chrono;

DataType getMaximumDeath (PersistenceDiagram diag)
{
  auto max = DataType(0);
  for (auto point : diag)
  {
    if (point.y() == infty)
    {
      continue;
    }
    else if (point.y() > max)
    {
      max = point.y();
    }
  }
  return max;
}


int main ( int argc, char** argv )
{
  std::vector<std::string> filenames;
  if (argc > 2)
  {
    filenames.assign(argv + 1, argv + argc);
  }
  else
  {
    throw std::invalid_argument("Wrong number of Arguments");
  }

  unsigned n = argc - 1;

  //DataType sphere_r = 0.63;
  
  //DataType torus_r = 0.5;
  //DataType torus_R = 0.25;

  std::cerr << "* calculating persistence Landscapes" << "\n";

  std::vector<PersistenceLandscape> landscapeV;
  
  std::ofstream normout;
  normout.open("/home/jens/Uni/data_topology/Project/syntheticResults/landscapeNorms/norm.dat");
  
  for( unsigned i = 0; i < n; i++ )
  {
    PersistenceDiagram diagram = aleph::io::load<DataType>(filenames[i]);

    auto max = getMaximumDeath(diagram);

    std::cout<< " - loading persistence diagram " << i << " with  betti number: " << diagram.betti() << std::endl;

    PersistenceLandscape landsc = PersistenceLandscape(diagram, 0, 4 * max);
    landscapeV.push_back(landsc);
    
    std::string filename = "/home/jens/Uni/data_topology/Project/syntheticResults/landscapesFromfixedSample/landscape_"
                         + std::to_string(i)
                         + ".dat";
    //filename << "landscape_" << i << ".dat";
    landsc.fileOutput(filename);
    
    //  + std::to_string(i) + ".dat");
    normout << PersistenceLandscape(diagram, 0, 4 * max).norm(2) << std::endl;
  }
  //PersistenceLandscape mean = std::accumulate(landscapeV.begin(), landscapeV.end(), PersistenceLandscape());
  auto mean = std::accumulate(landscapeV.begin(), landscapeV.end(), PersistenceLandscape());
  mean *= (1/static_cast<double>(landscapeV.size()));

  mean.fileOutput("/home/jens/Uni/data_topology/Project/syntheticResults/landscapesFromfixedSample/mean.dat");

  /*
  std::cerr << "* calculating mean landscape" << "\n";

  // calculate mean landscape from given sample
  auto mean = std::accumulate(landscapeV.begin(), landscapeV.end(), PersistenceLandscape());
  mean *= (1/static_cast<double>(landscapeV.size()));

  std::string object;
  if (sampleFromSphere)
  {
    object = "sphere";
  }
  else if(sampleFromTorus)
  {
    object = "torus";
  }
  auto filename = "mean_"
                       + object
                       + "_r="
                       + std::to_string(r)
                       + "_R="
                       + std::to_string(R)
                       + "_landscapes="
                       + std::to_string(n)
                       + "_points="
                       + std::to_string(m);
                       + ".dat"
  mean.fileOutput(filename);
  */

  // calculate landscape distance matrix
  std::vector<std::vector<DataType>> distanceMatrix (n,std::vector<DataType>(n,0));

  // File output for distance matrix
  std::string filename = "landscapeDistanceMatrix.dat";

  std::ofstream matrixFile;
  matrixFile.open(filename);

  for(size_t i = 0; i < n; i++)
  {
    for(size_t j = 0; j < n; j++)
    {
      auto difference = landscapeV[i] - landscapeV[j];
      distanceMatrix[i][j] = difference.norm(2);
      matrixFile << distanceMatrix[i][j] << " ";
      
    }
    matrixFile << "\n";
  }
  
  return 0;
}
