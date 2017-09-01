#include <aleph/geometry/RipsExpander.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>
#include <aleph/topology/filtrations/Degree.hh>

#include <aleph/topology/io/SparseAdjacencyMatrix.hh>

#include <aleph/utilities/Format.hh>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main( int argc, char** argv )
{
  using DataType          = unsigned;
  using VertexType        = std::size_t;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  if( argc <= 1 )
    return -1;

  std::string filename = argv[1];

  std::vector<SimplicialComplex> simplicialComplexes;
  std::vector<std::string> labels;

  aleph::topology::io::SparseAdjacencyMatrixReader reader;
  reader.setReadGraphLabels();

  std::cerr << "* Reading '" << filename << "'...";

  reader( filename, simplicialComplexes );

  std::cerr << "finished\n"
            << "* Read " << simplicialComplexes.size() << " simplicial complexes\n";

  // Expand simplicial complexes ---------------------------------------

  aleph::geometry::RipsExpander<SimplicialComplex> expander;

  // TODO: make expansion configurable; does it make sense to expand the
  // complexes at all?

  // Calculate degrees -------------------------------------------------

  std::cerr << "* Calculating degree-based filtration...";

  for( auto&& K : simplicialComplexes )
  {
    std::vector<DataType> degrees;
    aleph::topology::filtrations::degrees( K, std::back_inserter( degrees ) );

    K = expander.assignMaximumData( K, degrees.begin(), degrees.end() );

    K.sort( aleph::topology::filtrations::Data<Simplex>() );

    std::cout << K << "\n";
  }

  std::cerr << "finished\n";

  // Calculate persistent homology -------------------------------------

  {
    std::size_t index = 0;

    for( auto&& K : simplicialComplexes )
    {
      bool dualize                    = true;
      bool includeAllUnpairedCreators = true;

      auto diagrams
        = aleph::calculatePersistenceDiagrams( K,
                                               dualize,
                                               includeAllUnpairedCreators );

      for( auto&& diagram : diagrams )
      {
        diagram.removeDiagonal();

        auto output = "/tmp/"
                      + aleph::utilities::format( index, simplicialComplexes.size() )
                      + "_d"
                      + std::to_string( diagram.dimension() )
                      + ".txt";

        std::ofstream out( output );
        out << diagram << "\n";
      }

      ++index;
    }
  }
}