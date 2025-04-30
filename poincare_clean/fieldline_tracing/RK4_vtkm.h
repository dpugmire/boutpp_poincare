#pragma once

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSet.h>

#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <iostream>

void
writeArray1DToFile(std::vector<double>& array, const std::string& fname);

void
writeArray2DToFile(std::vector<std::vector<double>>& array, const std::string& fname);

void
writeArray3DToFile(const std::vector<std::vector<std::vector<double>>>& array, const std::string& fname);

vtkm::Vec3f
RK4_FLT1_vtkm(const vtkm::Vec3f& pStart,
            const vtkm::cont::DataSet& grid2D,
            const vtkm::cont::DataSet& grid2D_cfr,
            const vtkm::cont::DataSet& grid2D_xz,
            const vtkm::cont::DataSet& grid3D,
            const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& xarray,
            const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& zarray,
            int region,
            int direction, int nypf1, int nypf2,
            std::ofstream& rk4Out,
            int iline,
            int it,
            bool dumpFiles);
