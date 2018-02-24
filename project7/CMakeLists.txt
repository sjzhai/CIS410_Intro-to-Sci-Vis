cmake_minimum_required(VERSION 3.10)

PROJECT(proj6B)
SET(VTK_DIR /Users/shengjiezhai/Develop/VTK/VTKBuild)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj6B MACOSX_BUNDLE proj6B)

if(VTK_LIBRARIES)
  target_link_libraries(proj6B ${VTK_LIBRARIES})
else()
  target_link_libraries(proj6B vtkHybrid)
endif()


