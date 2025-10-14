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
