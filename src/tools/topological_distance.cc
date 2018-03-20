/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Given a set of persistence diagrams, it calculates various topological
  distances and returns a distance matrix.

  This tool can be helpful in different application scenarios:

  1. You want to determine the dissimilarity between two
     high-dimensional shapes, filtered by their distance
     function.

  2. You want to measure how a data descriptor, e.g. any
     density estimator, is changing over embeddings of a
     high-dimensional data set.

  3. You want to determine if certain samples of a space
     have the same characteristics than the original.

  The tool attempts to be smart and groups different inputs according to
  their common prefix. Currently, it only understands _d and _k as valid
  suffixes. Hence, the following input files are considered to belong to
  the same data set:

  - Test_d01
  - Test_d05
  - Test_d07

  Likewise:

  - Test_k1
  - Test_k7
  - Test_k9

  Please keep this in mind when using the tool.
*/

#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <string>
#include <vector>
#include <chrono>
#include <ctime>
#include <iomanip>

#include <cmath>

// TODO: Replace this as soon as possible with a more modern option
// parser interface.
#include <getopt.h>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/PersistenceIndicatorFunction.hh>
#include <aleph/persistenceDiagrams/PersistenceLandscape.hh>

#include <aleph/persistenceDiagrams/distances/Hausdorff.hh>
#include <aleph/persistenceDiagrams/distances/Wasserstein.hh>

#include <aleph/persistenceDiagrams/io/JSON.hh>
#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <aleph/utilities/Filesystem.hh>

using DataType                     = double;
using PersistenceDiagram           = aleph::PersistenceDiagram<DataType>;
using PersistenceIndicatorFunction = aleph::math::StepFunction<DataType>;
using PersistenceLandscape         = aleph::PersistenceLandscape<DataType>;

/*
  Auxiliary structure for describing a data set. I need this in order to
  figure out the corresponding dimension of the persistence diagram.
*/

struct DataSet
{
  std::string name;
  std::string filename;
  unsigned dimension;

  PersistenceDiagram persistenceDiagram;
  PersistenceIndicatorFunction persistenceIndicatorFunction;
};

/* Usage information */
void usage()
{
  std::cerr << "Usage: topological_distance [--power=POWER] [--kernel] [--exp] [--sigma]\n"
            << "                            [--hausdorff|indicator|wasserstein]\n"
            << "                            [--clean] FILES\n"
            << "\n"
            << "Calculates distances between a set of persistence diagrams, stored\n"
            << "in FILES. By default, this tool calculates Hausdorff distances for\n"
            << "all diagrams. This can be modified.\n"
            << "\n"
            << "If no other value is given, all distances are weighted using $p=2$\n"
            << "during the construction of a pairwise distance matrix. Furthermore\n"
            << "this tool can calculate kernels for use in kernel-based methods in\n"
            << "machine learning.\n"
            << "\n"
            << "The distance matrix is written to STDOUT. Rows and columns will be\n"
            << "separated by whitespace.\n"
            << "\n"
            << "This tool tries to be smart and is able to detect whether a set of\n"
            << "persistence diagrams belongs to the same group. This works only if\n"
            << "each file contains a suffix with digits that is preceded by either\n"
            << "a 'd' (for dimension) or a 'k' (for clique dimension).\n"
            << "\n"
            << "Flags:\n"
            << "  -c: clean persistence diagrams (remove unpaired points)\n"
            << "  -e: use exponential weighting for kernel calculation\n"
            << "  -h: calculate Hausdorff distances\n"
            << "  -i: calculate persistence indicator function distances\n"
            << "  -k: calculate kernel values instead of distances\n"
            << "  -n: normalize the persistence indicator function\n"
            << "  -s: use sigma as a scale parameter for the kernel\n"
            << "  -w: calculate Wasserstein distances\n"
            << "\n";
}

/*
  Stores a matrix in an output stream. The matrix is formatted such that
  individual values are separated by spaces and each row ends with '\n'.

  This format can be easily parsed by auxiliary programs such as gnuplot
  or R.
*/

void storeMatrix( const std::vector< std::vector<double> >& M, std::ostream& out )
{
  if( M.empty() )
    return;

  auto rows = M.size();
  auto cols = M.front().size();

  for( decltype(rows) row = 0; row < rows; row++ )
  {
    for( decltype(cols) col = 0; col < cols; col++ )
    {
      if( col != 0 )
        out << " ";

      out << M[row][col];
    }

    out << "\n";
  }
}

/*
  Calculates the topological distance between two data sets using persistence
  indicator functions. This requires enumerating all dimensions and finding a
  corresponding persistence indicator function. If no suitable function could
  be found, the calculation defaults to calculating the norm.
*/

double distancePIF( const std::vector<DataSet>& dataSet1,
                    const std::vector<DataSet>& dataSet2,
                    unsigned minDimension,
                    unsigned maxDimension,
                    double power,
                    bool normalize )
{
  auto getPersistenceIndicatorFunction = [] ( const std::vector<DataSet>& dataSet, unsigned dimension )
  {
    auto it = std::find_if( dataSet.begin(), dataSet.end(),
                            [&dimension] ( const DataSet& dataSet )
                            {
                              return dataSet.dimension == dimension;
                            } );

    if( it != dataSet.end() )
      return PersistenceIndicatorFunction( it->persistenceIndicatorFunction );
    else
      return PersistenceIndicatorFunction();
  };

  double d = 0.0;

  for( unsigned dimension = minDimension; dimension <= maxDimension; dimension++ )
  {
    auto f = getPersistenceIndicatorFunction( dataSet1, dimension );
    auto g = getPersistenceIndicatorFunction( dataSet2, dimension );

    if( normalize )
    {
      f = aleph::math::normalize( f );
      g = aleph::math::normalize( g );
    }

    g = -g;
    if( power == 1.0 )
      d = d + (f+g).abs().integral();
    else
      d = d + (f+g).abs().pow( power ).integral();
  }

  return d;
}

/*
  Calculates the topological distance between two data sets, using
  a standard distance between two persistence diagrams, for example
  the Hausdorff, Wasserstein, or bottleneck distance.

  By default, the Wasserstein distance is calculated.
*/

template <class Functor>
double persistenceDiagramDistance( const std::vector<DataSet>& dataSet1,
                                   const std::vector<DataSet>& dataSet2,
                                   unsigned minDimension,
                                   unsigned maxDimension,
                                   double power,
                                   Functor functor = [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2, double power )
                                   {
                                     return aleph::distances::wassersteinDistance( D1, D2, power );
                                   } )
{
  auto getPersistenceDiagram = [] ( const std::vector<DataSet>& dataSet, unsigned dimension )
  {
    auto it = std::find_if( dataSet.begin(), dataSet.end(),
                            [&dimension] ( const DataSet& dataSet )
                            {
                              return dataSet.dimension == dimension;
                            } );

    if( it != dataSet.end() )
      return PersistenceDiagram( it->persistenceDiagram );
    else
      return PersistenceDiagram();
  };

  double d = 0.0;

  for( unsigned dimension = minDimension; dimension <= maxDimension; dimension++ )
  {
    auto D1 = getPersistenceDiagram( dataSet1, dimension );
    auto D2 = getPersistenceDiagram( dataSet2, dimension );

    d += functor( D1, D2, power );
  }

  //d = std::pow( d, 1.0 / power );
  return d;
}


int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "power"      , required_argument, nullptr, 'p' },
    { "sigma"      , required_argument, nullptr, 's' },
    { "clean"      , no_argument      , nullptr, 'c' },
    { "exp"        , no_argument      , nullptr, 'e' },
    { "hausdorff"  , no_argument      , nullptr, 'h' },
    { "indicator"  , no_argument      , nullptr, 'i' },
    { "normalize"  , no_argument      , nullptr, 'n' },
    { "kernel"     , no_argument      , nullptr, 'k' },
    { "wasserstein", no_argument      , nullptr, 'w' },
    { "landscape"  , no_argument      , nullptr, 'l' },
    { nullptr      , 0                , nullptr,  0  }
  };

  double power                      = 2.0;
  double sigma                      = 1.0;
  bool cleanPersistenceDiagrams     = false;
  bool useExponentialFunction       = false;
  bool useIndicatorFunctionDistance = false;
  bool normalize                    = false;
  bool calculateKernel              = false;
  bool useWassersteinDistance       = false;
  bool useLandscapeDistance         = false;

  int option = 0;
  while( ( option = getopt_long( argc, argv, "p:s:cehinkw", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 'p':
      power = std::stod( optarg );
      break;
    case 's':
      sigma = std::stod( optarg );
      break;
    case 'c':
      cleanPersistenceDiagrams = true;
      break;
    case 'e':
      useExponentialFunction = true;
      break;
    case 'h':
      useWassersteinDistance       = false;
      useIndicatorFunctionDistance = false;
      useLandscapeDistance         = false;
      break;
    case 'i':
      useIndicatorFunctionDistance = true;
      useWassersteinDistance       = false;
      useLandscapeDistance         = false;
      break;
    case 'k':
      calculateKernel = true;
      break;
    case 'n':
      normalize = true;
      break;
    case 'w':
      useIndicatorFunctionDistance = false;
      useWassersteinDistance       = true;
      useLandscapeDistance         = false;
      break;
    case 'l':
      useWassersteinDistance       = false;
      useIndicatorFunctionDistance = false;
      useLandscapeDistance         = true;
    default:
      break;
    }
  }

  if( ( argc - optind ) <= 1 )
  {
    usage();
    return -1;
  }

  std::vector< std::vector<DataSet> > dataSets;

  // Get filenames & prefixes ------------------------------------------

  unsigned minDimension = std::numeric_limits<unsigned>::max();
  unsigned maxDimension = 0;

  {
    std::vector<std::string> filenames;
    filenames.reserve( argc - 1 );

    for( int i = optind; i < argc; i++ )
      filenames.push_back( argv[i] );

    // This should never happen...
    if( filenames.empty() )
      return -1;

    // If the first filename is a text file, I am assuming that the rest
    // of them also are. The program will then read all diagrams, try to
    // match them to a dimension, and store them.
    if( aleph::utilities::extension( filenames.front() ) == ".txt" )
    {
      // TODO: Make this regular expression more, well, 'expressive' and
      // support more methods of specifying a dimension.
      std::regex reDataSetPrefix( "(.*)_[dk]([[:digit:]]+)\\.txt" );
      std::smatch matches;

      unsigned index = 0;

      // Maps filenames to indices. I need this to ensure that the internal
      // ordering of files coincides with the shell's ordering.
      std::map<std::string, unsigned> filenameMap;

      for( auto&& filename : filenames )
      {
        auto name = filename;

        // Check whether the name contains a recognizable prefix and
        // suffix. If not, use the complete filename to identify the
        // data set.
        if( std::regex_match( filename, matches, reDataSetPrefix ) ){
          name = matches[1];
          //std::cerr << "**** YaY ****\n";
        }
        if( filenameMap.find( name ) == filenameMap.end() )
          filenameMap[ name ] = index++;
      }

      dataSets.resize( filenameMap.size() );

      for( auto&& filename : filenames )
      {
        auto name      = filename;
        auto dimension = 0u;

        // Check if a recognizable prefix and suffix exist so that we may
        // grab information about the data set and its dimension. If not,
        // use the complete filename to identify the data set.
        if( std::regex_match( filename, matches, reDataSetPrefix ) )
        {
          name      = matches[1];
          dimension = unsigned( std::stoul( matches[2] ) );
        }

        dataSets.at( filenameMap[name] ).push_back( { name, filename, dimension, {}, {} } );

        minDimension = std::min( minDimension, dimension );
        maxDimension = std::max( maxDimension, dimension );
      }

      // Load persistence diagrams & calculate indicator functions -----

      for( auto&& sets : dataSets )
      {
        for( auto&& dataSet : sets )
        {
          std::cerr << "* Processing '" << dataSet.filename << "'...";

          dataSet.persistenceDiagram = aleph::io::load<DataType>( dataSet.filename );

          if( cleanPersistenceDiagrams )
          {
            dataSet.persistenceDiagram.removeDiagonal();
            dataSet.persistenceDiagram.removeUnpaired();
          }

          // FIXME: This is only required in order to ensure that the
          // persistence indicator function has a finite integral; it
          // can be solved more elegantly by using a special value to
          // indicate infinite intervals.
          auto pd = dataSet.persistenceDiagram;
          pd.removeUnpaired();

          dataSet.persistenceIndicatorFunction
             = aleph::persistenceIndicatorFunction( pd );

          std::cerr << "finished\n";
        }
      }
    }
    else if( aleph::utilities::extension( filenames.front() ) == ".json" )
    {
      dataSets.reserve( filenames.size() );

      for( auto&& filename : filenames )
      {
        auto persistenceDiagrams = aleph::io::readJSON<DataType>( filename );

        std::vector<DataSet> dataSet;
        dataSet.reserve( persistenceDiagrams.size() );

        for( auto&& diagram : persistenceDiagrams )
        {
          auto dimension = static_cast<unsigned>( diagram.dimension() );
          minDimension   = std::min( minDimension, dimension );
          maxDimension   = std::max( maxDimension, dimension );

          auto name  = aleph::utilities::stem( filename );
          name      += "_";
          name      += "d" + std::to_string( diagram.dimension() );

          // FIXME: This is only required in order to ensure that the
          // persistence indicator function has a finite integral; it
          // can be solved more elegantly by using a special value to
          // indicate infinite intervals.
          auto pd = diagram;
          pd.removeUnpaired();

          if( cleanPersistenceDiagrams )
          {
            diagram.removeDiagonal();
            diagram.removeUnpaired();
          }

          dataSet.push_back( { name,
                               filename,
                               dimension,
                               diagram,
                               aleph::persistenceIndicatorFunction( pd ) } );
        }

        dataSets.push_back( dataSet );
      }
    }
  }

  // Setup distance functor --------------------------------------------

  std::function< double( const PersistenceDiagram&, const PersistenceDiagram&, double ) > functor
    = !useIndicatorFunctionDistance ? useWassersteinDistance ? [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2, double p )
                                                               {
                                                                 return aleph::distances::wassersteinDistance( D1, D2, p );
                                                               }
                                                             : useLandscapeDistance ? [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2, double p )
                                                                                      {
                                                                                        //std::cerr << "***YAY-Landscape***\n";
                                                                                        auto L1 = PersistenceLandscape(D1);
                                                                                        auto L2 = PersistenceLandscape(D2);
                                                                                        return (L2-L1).norm(2);
                                                                                      }
                                                                                    : [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2, double p )
                                                                                      {
                                                                                        return std::pow( aleph::distances::hausdorffDistance( D1, D2 ), p );
                                                                                      }
                                    : [] ( const PersistenceDiagram&, const PersistenceDiagram&, double )
                                    {
                                      return 0.0;
                                    };

  // Calculate all distances -------------------------------------------

  {
    auto name = useIndicatorFunctionDistance ? "persistence indicator function"
                                             : useWassersteinDistance ? "Wasserstein" 
                                                                      : useLandscapeDistance ? "Landscape"
                                                                                             : "Hausdorff";

    auto type = calculateKernel ? "kernel values" : "distances";

    std::cerr << "* Calculating pairwise " << type << " with " << name << " distance\n";
    std::cerr << "* Calculating pairwise " << type << " with p=" << power << "...";
  }

  std::vector< std::vector<double> > distances;
  distances.resize( dataSets.size(), std::vector<double>( dataSets.size() ) );
  
  // Measure Performance -----------------------------------------------
  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();
  
  #pragma omp parallel for collapse(2)
  for( std::size_t row = 0; row < dataSets.size(); row++ )
  {
    for( std::size_t col = 0; col < dataSets.size(); col++ )
    {
      if( row <= col )
        continue;

      double d = 0.0;

      if( useIndicatorFunctionDistance )
        d = distancePIF( dataSets.at(row), dataSets.at(col), minDimension, maxDimension, power, normalize );
      else
        d = persistenceDiagramDistance( dataSets.at(row), dataSets.at(col), minDimension, maxDimension, power, functor );

      if( calculateKernel )
      {
        d = -d;
        if( useExponentialFunction )
          d = std::exp( sigma * d );
      }

      distances[row][col] = d;
      distances[col][row] = d;
    }
  }
  
  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();
  
  std::ofstream performanceFile;
  performanceFile.open("performance.txt");
  performanceFile << std::fixed << std::setprecision(6) << "CPU time used: "
                  << (c_end-c_start) / (CLOCKS_PER_SEC * 60.0) << " m\n"
                  << "Wall clock time passed: "
                  << std::to_string(std::chrono::duration<double, std::ratio<60>>(t_end-t_start).count())
                  << " m\n";

  std::cerr << "finished\n";

  std::cerr << "Storing matrix...";

  storeMatrix( distances, std::cout );

  std::cerr << "finished\n";

  std::cerr << "Data sets were processed in the following order:\n";
  for( auto&& dataSet : dataSets )
  {
    if( !dataSet.empty() )
      std::cerr << "  - " << dataSet.front().name << "\n";
  }
}
