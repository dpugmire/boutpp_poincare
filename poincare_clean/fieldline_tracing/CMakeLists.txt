cmake_minimum_required(VERSION 3.10)
project(fieldline_tracing)

#find_package(NetCDF REQUIRED)

find_package(Viskores REQUIRED QUIET)
find_package(mpi REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(fieldline_tracing_vtkm fieldline_tracing_vtkm.cxx  NetCDFLoader.cxx  utils.cxx)


# Allow user-specified paths to override
set(NetCDF_INCLUDE_DIR "" CACHE PATH "NetCDF include directory")
set(NetCDF_LIBRARIES "" CACHE FILEPATH "NetCDF library file")

# Add include path and link libraries manually
target_include_directories(fieldline_tracing_vtkm PRIVATE ${NetCDF_INCLUDE_DIR})
target_link_libraries(fieldline_tracing_vtkm ${NetCDF_LIBRARIES})
target_link_libraries(fieldline_tracing_vtkm viskores::cont viskores::io viskores::filter_resampling)


add_executable(test_spline test_spline.cxx Spline.cxx)
target_link_libraries(test_spline viskores::cont viskores::io viskores::filter_resampling)

add_executable(test_spline3d test_spline3d.cxx Spline.cxx)
target_link_libraries(test_spline3d viskores::cont viskores::io viskores::filter_resampling)
