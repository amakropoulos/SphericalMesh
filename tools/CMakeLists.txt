add_mesh_executable(volume-to-mesh)

add_mesh_executable(mesh-to-sphere)
target_link_libraries(mesh-to-sphere fastMarchingGeodesic)