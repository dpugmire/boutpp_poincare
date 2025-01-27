#pragma once

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

std::pair<double, double> RK4_FLT1(
    double xStart, double yStart, double zStart,
    const std::vector<std::vector<std::vector<double>>>& dxdy,
    const std::vector<std::vector<std::vector<double>>>& dzdy,
    std::vector<double>& xarray,
    std::vector<double>& zarray,
    int region,
    const std::vector<std::vector<double>>& dxdy_pm1,
    const std::vector<std::vector<double>>& dzdy_pm1,
    int direction, int nypf1, int nypf2,
    std::ofstream& rk4Out,
    int iline,
    int it,
    bool dumpFiles);
