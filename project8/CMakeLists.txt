cmake_minimum_required(VERSION 3.10)

PROJECT(proj8)
SET(VTK_DIR /Users/shengjiezhai/Develop/VTK/VTKBuild)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj8 MACOSX_BUNDLE proj8)

if(VTK_LIBRARIES)
  target_link_libraries(proj8 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj8 vtkHybrid)
endif()


