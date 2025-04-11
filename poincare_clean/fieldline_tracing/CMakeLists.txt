cmake_minimum_required(VERSION 3.10)
project(fieldline_tracing)

find_package(NetCDF REQUIRED)
find_package(VTKm REQUIRED)
#set(Eigen3_DIR "/Users/dpn/sw/eigen/install/lib")  # Adjust to your installation path
#find_package(Eigen3 REQUIRED NO_DEFAULT_PATH)
#include_directories(/Users/dpn/sw/eigen/install/include/eigen3)

#set(TinySplinecxx_DIR "/Users/dpn/sw/tinyspline/include/lib64/cmake/tinysplinecxx")
set(tinysplinecxx_DIR "/Users/dpn/sw/tinyspline/install/lib64/cmake/tinysplinecxx")
set(TINY_SPLINE_INCLUDE_DIR "/Users/dpn/sw/tinyspline/install/include")
set(TINY_SPLINE_LIB "/Users/dpn/sw/tinyspline/install/lib64/libtinysplinecxx.a")
find_package(tinysplinecxx REQUIRED NO_DEFAULT_PATH)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${TINY_SPLINE_INCLUDE_DIR} /opt/homebrew/include)
link_directories(/opt/homebrew/lib)

add_subdirectory(alglib)


add_executable(fieldline_tracing fieldline_tracing.cxx SplineInterpolation.cxx NetCDFLoader.cxx RK4.cxx)
add_executable(fieldline_tracing_vtkm fieldline_tracing_vtkm.cxx SplineInterpolation.cxx NetCDFLoader.cxx RK4.cxx)

#target_link_libraries(fieldline_tracing NetCDF::NetCDF)
target_link_libraries(fieldline_tracing netcdf alglib ${TINY_SPLINE_LIB})
target_link_libraries(fieldline_tracing_vtkm netcdf alglib ${TINY_SPLINE_LIB} vtkm::cont vtkm::io vtkm::filter_resampling)
