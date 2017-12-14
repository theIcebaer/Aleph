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

#include <aleph/persistentHomology/Calculation.hh>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <algorithm>
#include <chrono>

#include <cmath>

#include <getopt.h>

// We first have to specify the data type of the persistence diagram,
// i.e. the type that is used by its individual points.
using DataType              = double;
using Distance              = aleph::distances::Euclidean<DataType>;
using PointCloud            = aleph::containers::PointCloud<DataType>;
using PersistenceDiagram    = aleph::PersistenceDiagram<DataType>;
using PersistenceLandscape  = aleph::PersistenceLandscape<DataType>;

namespace ch = std::chrono;

auto infty = std::numeric_limits<DataType>::infinity();


/*helper function for getting the maximum paired death point of a persistence diagram
 */
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


/**
  Auxiliary function for creating a random persistence diagram based on
  random samples from a torus. The torus is described by two radii, one
  for the outer part, the other for the inner part.
  Due to the sampling technique used, the specified number is merely an
  upper bound for the number of points that are to be sampled. Moreover
  the function will automatically handle Vietoris--Rips expansion.
*/

PersistenceDiagram createRandomTorusPersistenceDiagram( DataType R, DataType r, unsigned n )
{
  auto pointCloud = aleph::geometry::makeTorus(
    aleph::geometry::torusRejectionSampling( R, r, n ),
    R, r
  );
  std::ofstream pointfileTorus;
  pointfileTorus.open("/home/jens/Uni/data_topology/Project/syntheticResults/samples/torus/sampleTorus.txt");
  pointfileTorus << pointCloud;
  pointfileTorus.close();
  
  aleph::geometry::BruteForce<PointCloud, Distance> bruteForceWrapper( pointCloud );

  auto K
    = aleph::geometry::buildVietorisRipsComplex( bruteForceWrapper,2.5  * r, 2 );

  auto diagrams
    = aleph::calculatePersistenceDiagrams( K );

  diagrams.at(1).removeDiagonal();

  // We are only interested in the one-dimensional persistent homology
  // of the samples.
  return diagrams.at(1);
}

/**
  Auxiliary function for creating a random persistence diagram based on
  random samples from a sphere with a given radius.
  Note that this function automatically handles Vietoris--Rips
  expansion.
*/

PersistenceDiagram createRandomSpherePersistenceDiagram( DataType r, unsigned n )
{
  auto pointCloud = aleph::geometry::makeSphere(
    aleph::geometry::sphereSampling<DataType>( n ),
    r
  );
  std::ofstream pointfileSphere;
  pointfileSphere.open("/home/jens/Uni/data_topology/Project/syntheticResults/samples/sphere/sampleSphere.txt");
  pointfileSphere << pointCloud;
  pointfileSphere.close();

  aleph::geometry::BruteForce<PointCloud, Distance> bruteForceWrapper( pointCloud );

  auto K
    = aleph::geometry::buildVietorisRipsComplex( bruteForceWrapper, r, 2 );

  auto diagrams
    = aleph::calculatePersistenceDiagrams( K );

  diagrams.at(1).removeDiagonal();

  // We are only interested in the one-dimensional persistent homology
  // of the samples.
  return diagrams.at(1);
}



PersistenceLandscape createRandomSpherePersistenceLandscape ( DataType r, unsigned m)
{
  PersistenceDiagram diagram = createRandomSpherePersistenceDiagram( r, m);
  PersistenceLandscape landscape(diagram);
  
  return landscape;
}

PersistenceLandscape createRandomTorusPersistenceLandscape ( DataType R, DataType r, unsigned m)
{
  PersistenceDiagram diagram = createRandomTorusPersistenceDiagram(r, R,m);
  PersistenceLandscape landscape(diagram);
  
  return landscape;
}

int main ( int argc, char** argv )
{
  
    
  static option commandLineOptions[] = {
      { "m"     , required_argument, nullptr, 'm' },
      { "n"     , required_argument, nullptr, 'n' },
      { "R"     , required_argument, nullptr, 'R' },
      { "r"     , required_argument, nullptr, 'r' },
      { "sphere", no_argument      , nullptr, 's' },
      { "torus" , no_argument      , nullptr, 't' },
      { nullptr , 0                , nullptr,  0  }
  };
  
  
  unsigned m           = 50;
  unsigned n           = 50;
  DataType R           = DataType(0.25);
  DataType r           = DataType(0.50);

  bool sampleFromSphere = false;
  bool sampleFromTorus  = false;

  int option = 0;
  while( ( option = getopt_long( argc, argv, "m:n:R:r:t", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 'm':
      m = static_cast<unsigned>( std::stoul(optarg) );
      break;
    case 'n':
      n = static_cast<unsigned>( std::stoul(optarg) );
      break;
    case 'R':
      R = static_cast<DataType>( std::stod(optarg) );
      break;
    case 'r':
      r = static_cast<DataType>( std::stod(optarg) );
      break;
    case 's':
      sampleFromSphere = true;
      sampleFromTorus  = false;
      break;
    case 't':
      sampleFromSphere = false;
      sampleFromTorus  = true;
      break;
    default:
      break;
    }
  }

  std::cerr << "* Sampling " << n << " persistence diagrams\n";
  if( sampleFromSphere )
    std::cerr << "* Sampling " << m << " points from a sphere with r=" << r << "\n";
  else if( sampleFromTorus )
    std::cerr << "* Sampling at most " << m << " points from a torus with R=" << R << " and r=" << r << "\n";
  else
    std::cerr << "* Generating " << m << " random points per diagram\n";

  std::vector<PersistenceLandscape> landscapeV;
  PersistenceDiagram diagram;
  if (sampleFromSphere):
  {
    std::string normPathSphere = "/home/jens/Uni/data_topology/Project/syntheticResults/landscapeNorms/norm_sphere_r="
                               + std::to_string(r)
                               + ".txt";
    std::ofstream normFileSphere(normPathSphere);
  }
  else if (sampleFromTorus)
  {
    std::string normPathTorus = "/home/jens/Uni/data_topology/Project/syntheticResults/landscapeNorms/norm_torus_r="
                              + std::to_string(r)
                              + "_R="
                              + std::to_string(R)
                              + ".txt";
    std::ofstream normFileTorus(normPathTorus);
  }

                          
  
  for( unsigned i = 0; i < n; i++ )
  {
    if (sampleFromSphere)
    {
      // ch::high_resolution_clock::time_point t1 = ch::high_resolution_clock::now();
      diagram = createRandomSpherePersistenceDiagram(r, m);      
      std::cout<< "betti number in step " << i << ": " << diagram.betti() << std::endl;
      
      std::string diagramFileName = "/home/jens/Uni/data_topology/Project/syntheticResults/diagrams/d_"
                                  + std::to_string(diagram.dimension())
                                  + ".txt";
      std::ofstream diagramFile;
      diagramFile.open(diagramFileName);
      diagramFile << diagram;
      
      //std::ofstream normOut;
      //normout.open("/home/jens/Uni/data_topology/Project/syntheticResults/landscapeNorms/norm_"  + std::to_string( PersistenceLandscape(diagram, 0, r).norm(2));
      //diagram.removeUnpaired();
      // ch::high_resolution_clock::time_point t2 = ch::high_resolution_clock::now();
      // auto createDiagDuration = ch::duration_cast<ch::milliseconds>(t2 - t1).count();
      // std::cout << "diagram creation time in step " << i << ": " << createDiagDuration << " ms" << std::endl;
      // ch::high_resolution_clock::time_point t3 = ch::high_resolution_clock::now();
      auto max = getMaximumDeath(diagram);
 
      landscapeV.push_back( PersistenceLandscape(diagram, 0, 4 * max) );
      normFileSphere << PersistenceLandscape(diagram, 0, 10 * max).norm(2) << std::endl;
// ch::high_resolution_clock::time_point t4 = ch::high_resolution_clock::now();
      // auto createLandscapeDuration = ch::duration_cast<ch::microseconds>(t4 - t3).count();
      // std::cout << "landscape creation time in step " << i << ": " << createLandscapeDuration << " mys" << std::endl;
    }
    else if (sampleFromTorus)
    {
      diagram = createRandomTorusPersistenceDiagram(R, r, m);
      std::cout<< "betti number in step " << i << ": " << diagram.betti() << std::endl;
      //diagram.removeUnpaired();
      auto max = getMaximumDeath(diagram);

      landscapeV.push_back(PersistenceLandscape(diagram, 0, 10 * max));
      normFileTorus << PersistenceLandscape(diagram, 0, 10 * max).norm(2) << std::endl;
    }
    // std::stringstream filename;
    // filename << "landscape_" << i << ".dat";
    // PersistenceLandscape(diagram).fileOutput(filename.str());
  }
  
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
  auto filename = "/home/jens/Uni/data_topology/Project/syntheticResults/
                       + object 
                       + "/mean_"
                       + object
                       + "_r="
                       + std::to_string(r)
                       + "_R="
                       + std::to_string(R)
                       + "_landscapes="
                       + std::to_string(n)
                       + "_points="
                       + std::to_string(m)
                       + ".dat";
  mean.fileOutput(filename);
  //auto norm = mean.norm(2);


  // calculate landscape distance matrix
  std::vector<std::vector<DataType>> distanceMatrix (n, std::vector<DataType>(n,0));

  // File output for distance matrix
  std::string matrixFilename = "/home/jens/Uni/data_topology/Project/syntheticResults/distanceMatrix.txt";
  std::ofstream matrixFile;
  matrixFile.open(matrixFilename);

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
