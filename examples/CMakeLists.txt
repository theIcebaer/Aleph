
# FIXME: Ugly hack to make generic lambdas work
SET_SOURCE_FILES_PROPERTIES(
  create_persistence_landscapes.cc,
  persistence_landscape_from_input.cc
  PROPERTIES COMPILE_FLAGS "-std=c++14"
)

ADD_EXECUTABLE( 1D                               1D.cc                               )
ADD_EXECUTABLE( create_persistence_diagrams      create_persistence_diagrams.cc      )
ADD_EXECUTABLE( create_random_graph              create_random_graph.cc              )
ADD_EXECUTABLE( network_analysis                 network_analysis.cc                 )
ADD_EXECUTABLE( persistent_intersection_homology persistent_intersection_homology.cc )
ADD_EXECUTABLE( ply                              ply.cc                              )
ADD_EXECUTABLE( vtk                              vtk.cc                              )
ADD_EXECUTABLE( vietoris_rips                    vietoris_rips.cc                    )
ADD_EXECUTABLE( vietoris_rips_eccentricity       vietoris_rips_eccentricity.cc       )
ADD_EXECUTABLE( create_persistence_landscapes    create_persistence_landscapes.cc    )
ADD_EXECUTABLE( persistence_landscape_from_input persistence_landscape_from_input.cc    )

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-O3" )
