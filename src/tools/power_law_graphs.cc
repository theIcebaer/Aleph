#include <aleph/geometry/RipsExpander.hh>

#include <aleph/topology/io/EdgeLists.hh>
#include <aleph/topology/io/GML.hh>
#include <aleph/topology/io/Pajek.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>
#include <aleph/topology/filtrations/Degree.hh>

#include <aleph/utilities/Filesystem.hh>
#include <aleph/utilities/Format.hh>

#include <aleph/persistenceDiagrams/PersistenceLandscape.hh>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <typeinfo>

//#include <boost/filesystem.hpp>
//#include <boost/range/iterator_range.hpp>


int main( /*int argc, char** argv*/ )
{
  using DataType              = double;
  using VertexType            = std::size_t;
  using Simplex               = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex     = aleph::topology::SimplicialComplex<Simplex>;
  using PersistenceLandscape  = aleph::PersistenceLandscape<DataType>;

  //namespace fs = boost::filesystem;
  
  //if( argc <= 1 )
  //  return -1;

  //std::string directory = argv[1];
  

  std::vector<SimplicialComplex> simplicialComplexes;

  aleph::topology::io::GMLReader reader;

/*if(!fs::is_directory(directory))
  {
    std::cerr << directory << "is not a directory. \n";
  }

  for(auto& file : boost::make_iterator_range(fs::directory_iterator(directory)))
  {
    std::stringstream filename;
    filename << file;
    std::cerr << "* Reading '" << file << "'...";

    simplicialComplexes.push_back( SimplicialComplex() );
    reader( filename.str(), simplicialComplexes.back() );

    std::cerr << "finished\n"
              << "* Read " << simplicialComplexes.size() << " simplicial complexes\n";
  }
*/
  for(unsigned i = 0; i <= 99; i++)
  {
    std::string filename = "/home/jens/Uni/data_topology/Graph_Data/Synthetic-PowerLaw/pdg_"
                         + std::to_string(i)
                         + ".gml";

    std::cerr << "* Reading '" << filename << "'...";

    simplicialComplexes.push_back( SimplicialComplex() );
    reader( filename, simplicialComplexes.back() );

    std::cerr << "finished\n"
              << "* Read " << simplicialComplexes.size() << " simplicial complexes\n";
  }

  // Expand simplicial complexes ---------------------------------------

  aleph::geometry::RipsExpander<SimplicialComplex> expander;

  // TODO: make expansion configurable; does it make sense to expand the
  // complexes at all?

  // Calculate degrees -------------------------------------------------

  DataType maxDegree = 0;
  
//  std::vector<SimplicialComplex> testComplexes;
//  std::copy_if(simplicialComplexes.begin(), simplicialComplexes.end(), std::back_inserter(testComplexes), [](SimplicialComplex c){ return c.size() < 300; });
//  simplicialComplexes = testComplexes;
  
  std::cerr << "* Calculating degree-based filtration...";

  for( auto&& K : simplicialComplexes )
  {
    std::vector<DataType> degrees;
    aleph::topology::filtrations::degrees( K, std::back_inserter( degrees ) );

    // TODO: check degree for isolated vertices?

    if( !degrees.empty() )
    {
      maxDegree
        = std::max( maxDegree,
                    *std::max_element( degrees.begin(), degrees.end() ) );
    }

    // Wähle filtrationsfunktion ??? -----------------------------------
    K = expander.assignMaximumData( K, degrees.begin(), degrees.end() );
    //K = expander.assignData( K, degrees.begin(), degrees.end(), DataType(0), [] ( DataType a, DataType b ) { return a+b; } );
    K.sort( aleph::topology::filtrations::Data<Simplex>() );
  }

  std::cerr << "finished\n"
            << "* Identified maximum degree as D=" << maxDegree << "\n";

  // Store graphs ------------------------------------------------------
  
  /*
  {
    aleph::topology::io::GMLWriter writer;

    for( std::size_t i = 0; i < simplicialComplexes.size(); i++ )
    {
      auto filename = "/tmp/"
                      + aleph::utilities::format( i, simplicialComplexes.size() )
                      + labels[i]
                      + ".gml";

      std::cerr << "* Storing graph in '" << filename << "'...";

      writer( filename, simplicialComplexes[i] );

      std::cerr << "finished\n";
    }
  }
  */

  // Calculate persistent homology -------------------------------------

  {
    std::size_t index = 0;

    std::vector<std::vector<PersistenceLandscape>> landscapes_a(2);
    std::vector<std::vector<PersistenceLandscape>> landscapes_b(2);
    
    

    for( auto&& K : simplicialComplexes )
    {
      // auto i = &K - &simplicialComplexes[0];
      
      bool dualize                    = true;
      bool includeAllUnpairedCreators = true;

      auto diagrams
        = aleph::calculatePersistenceDiagrams( K,
                                               dualize,
                                               includeAllUnpairedCreators );
      
      for( auto&& diagram : diagrams )
      {
        diagram.removeDiagonal();
        auto dim = diagram.dimension();
        
        auto output = "/tmp/"
                      + aleph::utilities::format( index, simplicialComplexes.size() )
                      + "_d"
                      + std::to_string( diagram.dimension() )
                      + ".txt";

        std::ofstream out( output );

        for( auto&& point : diagram )
        {
          if( point.isUnpaired() )
            out << point.x() << "\t" << 2 * maxDegree << "\n";
          else
            out << point.x() << "\t" << point.y() << "\n";
        }
        
        auto landscapeFile = "/tmp/"
                           + aleph::utilities::format( index, simplicialComplexes.size() )
                           + "_l"
                           + std::to_string( diagram.dimension() )
                           + ".dat";
        PersistenceLandscape landscape(diagram, 0, 2 * maxDegree);

        landscapes_a[dim].push_back(landscape);
        landscape.fileOutput(landscapeFile);
      }
      ++index;
    }      
    
      // calculate landscape means
      
      std::cout << "landscape in dimension 0, has size: " << landscapes_a[0].size() << std::endl;
      std::cout << "landscape in dimension 1, has size: " << landscapes_a[1].size() << std::endl;
      //std::cout << "landscape with label  1, in dimension 0, has size: " << landscapes_b[0].size() << std::endl;
      //std::cout << "landscape with label  1, in dimension 1, has size: " << landscapes_b[1].size() << std::endl;
      
      std::sort(landscapes_a[0].begin(), landscapes_a[0].end(), [] (auto a, auto b) { return a.size() < b.size(); });
      std::sort(landscapes_a[1].begin(), landscapes_a[1].end(), [] (auto a, auto b) { return a.size() < b.size(); });

      
      std::transform(landscapes_a[0].begin(), landscapes_a[0].end(), landscapes_a[0].begin(), [] (auto a) { return (a += PersistenceLandscape()) *= 1/a.norm(2); });
      std::transform(landscapes_a[1].begin(), landscapes_a[1].end(), landscapes_a[1].begin(), [] (auto a) { return (a += PersistenceLandscape()) *= 1/a.norm(2); });

      
      auto mean_a_0 = std::accumulate(landscapes_a[0].begin(), landscapes_a[0].end(), PersistenceLandscape());
      mean_a_0 *= (1/static_cast<double>(landscapes_a[0].size()));

      std::cout << "calculated mean from label -1 data in dimension 0" << std::endl;
      
      auto mean_a_1 = std::accumulate(landscapes_a[1].begin(), landscapes_a[1].end(), PersistenceLandscape());
      mean_a_1 *= (1/static_cast<double>(landscapes_a[1].size()));

      std::cout << "calculated mean from label -1 data in dimension 1" << std::endl;
      
      // write means to /tmp/
      
      auto mean_a_0_File = "/tmp/mean_a_0.dat";
      auto mean_a_1_File = "/tmp/mean_a_1.dat";

      mean_a_0.fileOutput(mean_a_0_File);
      std::cout << "wrote mean a to" << mean_a_0_File << std::endl;
      mean_a_1.fileOutput(mean_a_1_File);
      std::cout << "wrote mean a to" << mean_a_1_File << std::endl;
  }

  // Store labels ------------------------------------------------------

  /*
  {
    std::vector<std::string> labels;
    reader.graphLabels( std::back_inserter( labels ) );

    std::ofstream out( "/tmp/"
                      + aleph::utilities::basename( filename )
                      + ".txt" );

    for( auto&& label : labels )
      out << label << "\n";
  }
  */
}
