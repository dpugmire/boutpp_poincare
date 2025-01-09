#pragma once

#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <iostream>

std::pair<float, float> RK4_FLT1(
    float xStart, float yStart, float zStart,
    const std::vector<std::vector<std::vector<float>>>& dxdy,
    const std::vector<std::vector<std::vector<float>>>& dzdy,
    std::vector<float>& xarray,
    std::vector<float>& zarray,
    int region,
    const std::vector<std::vector<float>>& dxdy_pm1,
    const std::vector<std::vector<float>>& dzdy_pm1,
    int direction, int nypf1, int nypf2);