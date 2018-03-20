
# FIXME: Ugly hack to make generic lambdas work
SET_SOURCE_FILES_PROPERTIES(
  create_persistence_landscapes.cc,
  persistence_landscape_from_input.cc,
  distance-matrix_from_diagrams.cc,
  mean-landscape_from_diagrams.cc
  
  PROPERTIES COMPILE_FLAGS "-std=c++14"
  COMPILE_FLAGS "-fopenmp"
  COMPILE_FLAGS "-O0"
)

ADD_EXECUTABLE( 1D                               1D.cc                               )
ADD_EXECUTABLE( create_persistence_diagrams      create_persistence_diagrams.cc      )
ADD_EXECUTABLE( create_random_graph              create_random_graph.cc              )
ADD_EXECUTABLE( network_analysis                 network_analysis.cc                 )
ADD_EXECUTABLE( network_analysis_landscape       network_analysis_landscape.cc       )
ADD_EXECUTABLE( persistent_intersection_homology persistent_intersection_homology.cc )
ADD_EXECUTABLE( ply                              ply.cc                              )
ADD_EXECUTABLE( vtk                              vtk.cc                              )
ADD_EXECUTABLE( vietoris_rips                    vietoris_rips.cc                    )
ADD_EXECUTABLE( vietoris_rips_eccentricity       vietoris_rips_eccentricity.cc       )
ADD_EXECUTABLE( create_persistence_landscapes    create_persistence_landscapes.cc    )
ADD_EXECUTABLE( persistence_landscape_from_input persistence_landscape_from_input.cc )
ADD_EXECUTABLE( distance-matrix_from_diagrams    distance-matrix_from_diagrams.cc    )
ADD_EXECUTABLE( mean-landscape_from_diagrams     mean-landscape_from_diagrams.cc     )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-O3" )
