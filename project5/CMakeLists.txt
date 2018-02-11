cmake_minimum_required(VERSION 3.10)

PROJECT(proj5)
SET(VTK_DIR /Users/shengjiezhai/Develop/VTK/VTKBuild)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj5 MACOSX_BUNDLE proj5)

if(VTK_LIBRARIES)
  target_link_libraries(proj5 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj5 vtkHybrid)
endif()


