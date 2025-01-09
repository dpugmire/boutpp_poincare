#include "RK4.h"
#include "SplineInterpolation.h"
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iostream>

extern float float_mod(float val, float mod_base);

float bilinear_interp2(const std::vector<float>& xarray,
              const std::vector<float>& zarray, 
              const std::vector<std::vector<float>>& data, 
              float x, 
              float z)
{
    // Validate inputs
    if (xarray.empty() || zarray.empty() || data.empty() || data[0].empty())
        throw std::invalid_argument("Input arrays must not be empty.");

    if (data.size() != xarray.size() || data[0].size() != zarray.size())
        throw std::invalid_argument("Data dimensions must match xarray and zarray sizes.");

    // Find indices for x and z
    auto x_it = std::lower_bound(xarray.begin(), xarray.end(), x);
    auto z_it = std::lower_bound(zarray.begin(), zarray.end(), z);

    if (x_it == xarray.end() || z_it == zarray.end())
        throw std::out_of_range("Interpolation point is out of bounds.");

    int i1 = std::max(0, static_cast<int>(x_it - xarray.begin() - 1));
    int j1 = std::max(0, static_cast<int>(z_it - zarray.begin() - 1));

    int i2 = std::min(i1 + 1, static_cast<int>(xarray.size() - 1));
    int j2 = std::min(j1 + 1, static_cast<int>(zarray.size() - 1));

    // Get the corners of the grid
    float x1 = xarray[i1], x2 = xarray[i2];
    float z1 = zarray[j1], z2 = zarray[j2];
    float q11 = data[i1][j1], q21 = data[i2][j1];
    float q12 = data[i1][j2], q22 = data[i2][j2];

    // Perform bilinear interpolation
    float interp = (q11 * (x2 - x) * (z2 - z) +
                    q21 * (x - x1) * (z2 - z) +
                    q12 * (x2 - x) * (z - z1) +
                    q22 * (x - x1) * (z - z1)) /
                   ((x2 - x1) * (z2 - z1));

    return interp;
}

std::vector<float>
flatten(std::vector<std::vector<float>>& input)
{
    std::vector<float> flattened;
    for (const auto& row : input)
    {
        flattened.insert(flattened.end(), row.begin(), row.end());
    }
    return flattened;
}
void
writeArray1DToFile(std::vector<float>& array, const std::string& fname, std::vector<std::size_t> dims)
{
    auto fname2 = fname + ".c.txt";
    std::ofstream out(fname2, std::ofstream::out);
    out<<"(";
    for (const auto& d : dims)
        out<<d<<",";
    out<<")"<<std::endl;

    out<<std::scientific<<std::setprecision(10);
    int cnt = 0;
    for (const auto& v : array)
    {
        out<<v<<std::endl;
        cnt++;
        if (cnt > 10000) break;
    }

    out.close();
}

void
writeArray1DToFile(std::vector<float>& array, const std::string& fname)
{
    std::vector<std::size_t> dims = {array.size()};
    writeArray1DToFile(array, fname, dims);
}

void
writeArray2DToFile(std::vector<std::vector<float>>& array, const std::string& fname)
{
    auto flatArray = flatten(array);
    std::vector<std::size_t> dims = {array.size(), array[0].size()};
    writeArray1DToFile(flatArray, fname, dims);
}

static std::vector<std::vector<float>>
Slice(const std::vector<std::vector<std::vector<float>>>& data, size_t idx)
{
    // Check bounds
    if (idx >= data[0].size())
    {
        throw std::out_of_range("ystart is out of bounds");
    }

    size_t dim1 = data.size();    // First dimension
    size_t dim3 = data[0][0].size(); // Third dimension

    // Create a 2D vector to hold the slice
    std::vector<std::vector<float>> dataSlice(dim1, std::vector<float>(dim3));

    // Extract the slice
    for (size_t i = 0; i < dim1; ++i)
    {
        for (size_t k = 0; k < dim3; ++k)
        {
            dataSlice[i][k] = data[i][idx][k];
        }
    }

    return dataSlice;
}

std::vector<std::vector<float>>
pad2DEdge(const std::vector<std::vector<float>>& input)
{
    if (input.empty() || input[0].empty())
    {
        throw std::invalid_argument("Input 2D vector must not be empty.");
    }

    // Create a new 2D vector with an additional column
    size_t rows = input.size();
    size_t cols = input[0].size();
    std::vector<std::vector<float>> padded(rows, std::vector<float>(cols + 1, 0.0f));

    // Copy original data and pad the last column
    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
            padded[i][j] = input[i][j];
        // Pad the last column with the last value of the current row
        padded[i][cols] = input[i][cols - 1];
    }

    return padded;
}

std::vector<std::vector<float>>
avgArrays(const std::vector<std::vector<float>>& x,
          const std::vector<std::vector<float>>& y)
{
    std::vector<std::vector<float>> result(x.size(), std::vector<float>(x[0].size(), 0.0f));
    for (size_t i = 0; i < x.size(); ++i)
        for (size_t j = 0; j < x[i].size(); ++j)
            result[i][j] = 0.5f * (x[i][j] + y[i][j]);

    return result;
}

std::pair<float, float> RK4_FLT1(
    float xStart, float yStart, float zStart,
    const std::vector<std::vector<std::vector<float>>>& dxdy,
    const std::vector<std::vector<std::vector<float>>>& dzdy,
     std::vector<float>& xarray,  std::vector<float>& zarray,
    int region,
    const std::vector<std::vector<float>>& dxdy_pm1,
    const std::vector<std::vector<float>>& dzdy_pm1,
    int direction, int nypf1, int nypf2)
{
    float hh = 0.5f;
    float h6 = 1.0f / 6.0f;

    bool dumpFiles = true;

    if (direction != 1 && direction != -1)
        throw std::invalid_argument("Direction parameter must be 1 or -1.");

    std::cout<<"dxdy.shape= "<<dxdy.size()<<" "<<dxdy[0].size()<<" "<<dxdy[0][0].size()<<std::endl;
    std::vector<std::vector<float>> dxdyp = Slice(dxdy, yStart);
    std::vector<std::vector<float>> dzdyp = Slice(dzdy, yStart);
    std::vector<std::vector<float>> dxdyn, dzdyn, dxdyh, dzdyh;
    std::cout<<"dxdyp.shape= "<<dxdyp.size()<<" "<<dxdyp[0].size()<<std::endl;

    if (direction == 1)
    {
        dxdyp = Slice(dxdy, yStart);
        dzdyp = Slice(dzdy, yStart);
        if (dumpFiles)
        {
            writeArray2DToFile(dxdyp, "dxdyp_1");
            writeArray2DToFile(dzdyp, "dzdyp_1");
        }
        if (region == 0 && yStart == nypf2)
        {
            dxdyn = dxdy_pm1;
            dzdyn = dzdy_pm1;
            dxdyh = avgArrays(dxdyp, dxdyn);
            dzdyh = avgArrays(dzdyp, dzdyn);
        }
        else
        {
            dxdyn = Slice(dxdy, yStart+1);
            dzdyn = Slice(dzdy, yStart+1);
            dxdyh = avgArrays(Slice(dxdy, yStart), Slice(dxdy, yStart+1));
            dzdyh = avgArrays(Slice(dzdy, yStart), Slice(dzdy, yStart+1));
            auto tmp0 = Slice(dxdy, yStart);
            auto tmp1 = Slice(dxdy, yStart+1);
            writeArray2DToFile(tmp0, "tmp0");
            writeArray2DToFile(tmp1, "tmp1");
            writeArray2DToFile(dxdyh, "sum0");
            writeArray2DToFile(dzdyh, "sum1");
        }
    }
    else
    {
        throw std::invalid_argument("Backwards RK4 not supported yet.");
    }

    if (dumpFiles)
    {
        writeArray1DToFile(xarray, "xarray");
        writeArray1DToFile(zarray, "zarray");
        writeArray2DToFile(dxdyp, "dxdyp");
        writeArray2DToFile(dxdyn, "dxdyn");
        writeArray2DToFile(dxdyh, "dxdyh");
        writeArray2DToFile(dzdyp, "dzdyp");
        writeArray2DToFile(dzdyn, "dzdyn");
        writeArray2DToFile(dzdyh, "dzdyh");
    }

    std::cout<<"Add the padding in Z!!! "<<std::endl;
    dxdyp = pad2DEdge(dxdyp);
    dxdyn = pad2DEdge(dxdyn);
    dxdyh = pad2DEdge(dxdyh);
    dzdyp = pad2DEdge(dzdyp);
    dzdyn = pad2DEdge(dzdyn);
    dzdyh = pad2DEdge(dzdyh);
    writeArray2DToFile(dxdyp, "dxdyp_");
    writeArray2DToFile(dxdyn, "dxdyn_");
    writeArray2DToFile(dxdyh, "dxdyh_");
    writeArray2DToFile(dzdyp, "dzdyp_");
    writeArray2DToFile(dzdyn, "dzdyn_");
    writeArray2DToFile(dzdyh, "dzdyh_");

    // Interpolation using the SplineInterpolation class
    //SplineInterpolation splineDxdy(xarray, zarray, dxdyp); //[zStart]);
    //SplineInterpolation splineDzdy(xarray, zarray, dzdyp); //[zStart]);

    //std::cout<<"RK4 step1: "<<xStart<<" "<<zStart<<" xx: "<<xarray.front()<<" "<<xarray.back()<<std::endl;
    printf("RK4 step1: %12.10e %12.10e\n", xStart, zStart);
    // RK4 Step 1
    float tmp1 = bilinear_interp2(xarray, zarray, dxdyp, xStart, zStart);
    float tmp2 = bilinear_interp2(xarray, zarray, dzdyp, xStart, zStart);

    float dxdy1 = bilinear_interp2(xarray, zarray, dxdyp, xStart, zStart);
    float dzdy1 = bilinear_interp2(xarray, zarray, dzdyp, xStart, zStart); //splineDzdy.evaluate(xStart);
    float x1 = xStart + direction * hh * dxdy1;
    float z1 = zStart + direction * hh * dzdy1;
    printf("  RES: %12.10e %12.10e\n", dxdy1, dzdy1);
    printf("  bi-RES: %12.10e %12.10e\n", tmp1, tmp2);
    //std::cout<<"   RES= "<<dxdy1<<" "<<dzdy1<<std::endl;

    // RK4 Step 2
    std::cout<<"RK4 step2: "<<x1<<" "<<z1<<" xz: "<<xarray[0]<<" "<<zarray[0]<<std::endl;

    float dxdy2 = bilinear_interp2(xarray, zarray, dxdyh, x1, std::fmod(z1, M_2_PI));
    float dzdy2 = bilinear_interp2(xarray, zarray, dzdyh, x1, std::fmod(z1, M_2_PI));
    float x2 = xStart + direction * hh * dxdy2;
    float z2 = zStart + direction * hh * dzdy2;
    std::cout<<"   RES= "<<dxdy2<<" "<<dzdy2<<std::endl;

    // RK4 Step 3
    std::cout<<"RK4 step3: "<<x2<<" "<<z2<<" xz: "<<xarray[0]<<" "<<zarray[0]<<std::endl;
    float dxdy3 = bilinear_interp2(xarray, zarray, dxdyh, x2, z2);
    float dzdy3 = bilinear_interp2(xarray, zarray, dzdyh, x2, z2);
    float x3 = xStart + direction * dxdy3;
    float z3 = zStart + direction * dzdy3;
    std::cout<<"   RES= "<<dxdy3<<" "<<dzdy3<<std::endl;

    // RK4 Step 4
    std::cout<<"RK4 step4: "<<x3<<" "<<z3<<" xz: "<<xarray[0]<<" "<<zarray[0]<<std::endl;
    float dxdy4 = bilinear_interp2(xarray, zarray, dxdyn, x3, z3);
    float dzdy4 = bilinear_interp2(xarray, zarray, dzdyn, x3, z3);
    std::cout<<"   RES= "<<dxdy4<<" "<<dzdy4<<std::endl;

    // Compute final x and z
    float xEnd = xStart + direction * h6 * (dxdy1 + 2.0f * dxdy2 + 2.0f * dxdy3 + dxdy4);
    float zEnd = zStart + direction * h6 * (dzdy1 + 2.0f * dzdy2 + 2.0f * dzdy3 + dzdy4);

    printf("   STEP= %12.10e %12.10e\n", xEnd, zEnd);

    return {xEnd, zEnd};
}
