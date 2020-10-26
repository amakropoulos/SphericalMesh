SET(HEADERS
  ../include/GW_Config.h
  ../include/GW_Face.h
  ../include/GW_Mesh.h
  ../include/GW_FaceIterator.h
  ../include/GW_Vertex.h
  ../include/GW_VertexIterator.h
  ../include/GW_SmartCounter.h
  ../include/GW_GeodesicFace.h
  ../include/GW_GeodesicMesh.h
  ../include/GW_GeodesicPath.h
  ../include/GW_GeodesicPoint.h
  ../include/GW_GeodesicVertex.h
  ../include/GW_TriangularInterpolation_Linear.h
  ../include/GW_TriangularInterpolation_Quadratic.h
  ../include/GW_TriangularInterpolation_Cubic.h
  ../include/vtkPolyDataGeodesicDistance.h
  ../include/vtkFastMarchingGeodesicDistance.h
  ../include/vtkFastMarchingGeodesicPath.h
  #vtkPolygonalSurfaceContourLineInterpolator2.h
 )
 
 set(SOURCES
  GW_Config.cpp
  GW_Face.cpp
  GW_FaceIterator.cpp
  GW_Mesh.cpp
  GW_SmartCounter.cpp
  GW_Vertex.cpp
  GW_VertexIterator.cpp
  GW_GeodesicFace.cpp
  GW_GeodesicMesh.cpp
  GW_GeodesicPath.cpp
  GW_GeodesicPoint.cpp
  GW_GeodesicVertex.cpp
  GW_TriangularInterpolation_Linear.cpp
  GW_TriangularInterpolation_Quadratic.cpp
  GW_TriangularInterpolation_Cubic.cpp
  vtkPolyDataGeodesicDistance.cxx
  vtkFastMarchingGeodesicDistance.cxx
  vtkFastMarchingGeodesicPath.cxx
  #vtkPolygonalSurfaceContourLineInterpolator2.cxx
)

set(DEPENDS
  ${VTK_LIBRARIES}
)

include_directories(../include)
ADD_LIBRARY(fastMarchingGeodesic ${HEADERS} ${SOURCES})
target_link_libraries(fastMarchingGeodesic ${VTK_LIBRARIES})
