SET(SPHERICALMESH_INCLUDES
../include/mirtk/M2SConnectivity.h
../include/mirtk/M2SDiffuser.h
../include/mirtk/M2SEdge.h
../include/mirtk/M2SGeodesicInterpolator.h
../include/mirtk/M2SMeshMerger.h
../include/mirtk/M2SOptimizer.h
../include/mirtk/M2SParameters.h
../include/mirtk/M2SPostProcessor.h
../include/mirtk/MeshToSphere.h
../include/mirtk/M2SRemesher.h
)

SET(SPHERICALMESH_SRCS
M2SConnectivity.cc
M2SDiffuser.cc
M2SGeodesicInterpolator.cc
M2SMeshMerger.cc
M2SOptimizer.cc
M2SParameters.cc
M2SPostProcessor.cc
MeshToSphere.cc
M2SRemesher.cc
)

include_directories(../ThirdParty/FastMarching/include)
ADD_LIBRARY(sphericalMesh ${SPHERICALMESH_INCLUDES} ${SPHERICALMESH_SRCS})
