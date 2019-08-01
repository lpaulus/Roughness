cmake_minimum_required (VERSION 2.8.11)
project (Roughness)

find_package(CGAL REQUIRED)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(Polyhedron)

add_executable(Roughness main.cpp extract_Vpropres.cpp)
add_executable(Color color.cpp)

if(CGAL_FOUND)
  #link it with the required CGAL libraries
  target_link_libraries(Roughness CGAL::CGAL)
else()
  message("ERROR: this program requires CGAL and will not be compiled.")
endif()
