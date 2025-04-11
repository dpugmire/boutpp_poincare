#include "RK4.h"
#include "SplineInterpolation.h"
#include "Evaluation.h"
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iostream>

extern double double_mod(double val, double mod_base);

double bilinear_interp2(const std::vector<double>& xarray,
              const std::vector<double>& zarray,
              const std::vector<std::vector<double>>& data,
              double x,
              double z)
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
    {
        std::cout<<"Interpolation point is out of bounds: "<<x<<" "<<xarray[0]<<" "<<xarray[xarray.size()-1]<<" "<<z<<" "<<zarray[0]<<" "<<zarray[zarray.size()-1]<<std::endl;
        throw std::out_of_range("Interpolation point is out of bounds.");

        if (x_it == xarray.end()) return xarray[xarray.size()-1];
        if (z_it == zarray.end()) return zarray[zarray.size()-1];
    }

    int i1 = std::max(0, static_cast<int>(x_it - xarray.begin() - 1));
    int j1 = std::max(0, static_cast<int>(z_it - zarray.begin() - 1));

    int i2 = std::min(i1 + 1, static_cast<int>(xarray.size() - 1));
    int j2 = std::min(j1 + 1, static_cast<int>(zarray.size() - 1));

    // Get the corners of the grid
    double x1 = xarray[i1], x2 = xarray[i2];
    double z1 = zarray[j1], z2 = zarray[j2];
    double q11 = data[i1][j1], q21 = data[i2][j1];
    double q12 = data[i1][j2], q22 = data[i2][j2];

    // Calculate the denominator and check for division by zero
    double denom = (x2 - x1) * (z2 - z1);
    if (denom == 0.0)
        throw std::runtime_error("Division by zero in bilinear interpolation: check grid spacing.");

    // Perform bilinear interpolation
    double interp = (q11 * (x2 - x) * (z2 - z) +
                    q21 * (x - x1) * (z2 - z) +
                    q12 * (x2 - x) * (z - z1) +
                    q22 * (x - x1) * (z - z1)) / denom;

    return interp;
}

std::vector<double>
flatten(std::vector<std::vector<double>>& input)
{
    std::vector<double> flattened;
    for (const auto& row : input)
    {
        flattened.insert(flattened.end(), row.begin(), row.end());
    }
    return flattened;
}

std::vector<double>
flatten(const std::vector<std::vector<std::vector<double>>>& input)
{
    std::vector<double> flattened;
    auto n0 = input.size();
    auto n1 = input[0].size();
    auto n2 = input[0][0].size();

    flattened.resize(n0*n1*n2);
    size_t cnt = 0;
    for (auto i = 0; i < n0; i++)
        for (auto j = 0; j < n1; j++)
            for (auto k = 0; k < n2; k++)
                flattened[cnt++] = input[i][j][k];

    return flattened;
}

void
writeArray1DToFile(std::vector<double>& array, const std::string& fname)
{
    auto fname2 = "/Users/dpn/" + fname + ".c.txt";
    std::ofstream out(fname2, std::ofstream::out);
    auto nx = array.size();
    out<<"("<<nx<<")"<<std::endl;

    out<<std::scientific<<std::setprecision(10);
    int cnt = 0;
    for (size_t i = 0; i < nx; i++)
    {
        if (cnt > 5000) break;
        out<<i<<", "<<array[i]<<std::endl;
        cnt++;
    }

    out.close();
}

void
writeArray2DToFile(std::vector<std::vector<double>>& array, const std::string& fname)
{
    auto fname2 = "/Users/dpn/" + fname + ".c.txt";
    std::ofstream out(fname2, std::ofstream::out);
    auto nx = array.size();
    auto ny = array[0].size();
    out<<"("<<nx<<", "<<ny<<")"<<std::endl;

    out<<std::scientific<<std::setprecision(10);
    size_t x0 = 0, x1 = nx, y0 = 0, y1 = ny;
    int cnt = 0;
    int maxCnt = 5000;
    maxCnt = -1;

    //x1 /= 2;
    //y1 /= 2;

    for (size_t i = x0; i < x1; i++)
        for (size_t j = y0; j < y1; j++)
        {
            if (maxCnt > 0 && cnt > maxCnt) break;
            out<<i<<", "<<j<<", "<<array[i][j]<<std::endl;
            cnt++;
        }
    out.close();
}

void
writeArray3DToFile(const std::vector<std::vector<std::vector<double>>>& array, const std::string& fname)
{
    auto val0 = array[123][79][101];
    auto val1 = array[192][47][200];

    auto fname2 = "/Users/dpn/" + fname + ".c.txt";
    std::ofstream out(fname2, std::ofstream::out);
    auto nx = array.size();
    auto ny = array[0].size();
    auto nz = array[0][0].size();
    out<<"("<<nx<<", "<<ny<<", "<<nz<<")"<<std::endl;

    out<<std::scientific<<std::setprecision(10);
    size_t x0 = 0, x1 = nx, y0 = 0, y1 = ny, z0 = 0, z1 = nz;

    x1 /= 2;
    y1 /= 2;
    z1 /= 2;

    int cnt = 0;
    for (size_t i = x0; i < x1; i++)
        for (size_t j = y0; j < y1; j++)
            for (size_t k = z0; k < z1; k++)
            {
                if (cnt > 5000) break;
                out<<i<<", "<<j<<", "<<k<<", "<<array[i][j][k]<<std::endl;
                cnt++;
            }
    out.close();
}

static std::vector<std::vector<double>>
Slice(const std::vector<std::vector<std::vector<double>>>& data, size_t idx)
{
    // Check bounds
    if (idx >= data[0].size())
    {
        throw std::out_of_range("ystart is out of bounds");
    }

    size_t dim1 = data.size();    // First dimension
    size_t dim3 = data[0][0].size(); // Third dimension

    // Create a 2D vector to hold the slice
    std::vector<std::vector<double>> dataSlice(dim1, std::vector<double>(dim3));

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

std::vector<std::vector<double>>
pad2DEdge(const std::vector<std::vector<double>>& input)
{
    if (input.empty() || input[0].empty())
    {
        throw std::invalid_argument("Input 2D vector must not be empty.");
    }

    // Create a new 2D vector with an additional column
    size_t rows = input.size();
    size_t cols = input[0].size();
    std::vector<std::vector<double>> padded(rows, std::vector<double>(cols + 1, 0.0f));

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

std::vector<std::vector<double>>
avgArrays(const std::vector<std::vector<double>>& x,
          const std::vector<std::vector<double>>& y)
{
    std::vector<std::vector<double>> result(x.size(), std::vector<double>(x[0].size(), 0.0f));
    for (size_t i = 0; i < x.size(); ++i)
        for (size_t j = 0; j < x[i].size(); ++j)
            result[i][j] = (x[i][j] + y[i][j]) / 2.0;

    return result;
}

static int count = 0;
double maxErrx1, maxErrz1, maxErrx2, maxErrz2;
double maxX1, maxZ1, maxX2, maxZ2;

void splineTest(double xStart, double zStart,
                std::vector<std::vector<double>>& dxdyp,
                std::vector<std::vector<double>>& dzdyp,
                std::vector<double>& xarray,  std::vector<double>& zarray)
{
    double linear_dxdy = bilinear_interp2(xarray, zarray, dxdyp, xStart, zStart);
    double linear_dzdy = bilinear_interp2(xarray, zarray, dzdyp, xStart, zStart);

    double spline1_dxdy = alglib_spline(xarray, zarray, dxdyp, xStart, zStart);
    double spline1_dzdy = alglib_spline(xarray, zarray, dzdyp, xStart, zStart);
    SplineInterpolation spline1dx(xarray, zarray, dxdyp), spline1dz(xarray, zarray, dzdyp);
    double spline2_dxdy = spline1dx.evaluate(xStart, zStart);
    double spline2_dzdy = spline1dz.evaluate(xStart, zStart);

    //double spline3_dxydy = interp2Spline(xarray, zarray, dxdyp, xStart, zStart);
    //double spline3_dzydy = interp2Spline(xarray, zarray, dzdyp, xStart, zStart);

    double errx1 = std::abs(linear_dxdy-spline1_dxdy), errz1 = std::abs(linear_dzdy-spline1_dzdy);
    double errx2 = std::abs(linear_dxdy-spline2_dxdy), errz2 = std::abs(linear_dzdy-spline2_dzdy);
    if (count == 0)
    {
        maxErrx1 = errx1;
        maxErrz1 = errz1;
        maxErrx2 = errx2;
        maxErrz2 = errz2;
        maxX1 = maxX2 = xStart;
        maxZ1 = maxZ2 = zStart;
    }
    if (errx1 > maxErrx1) {maxErrx1 = errx1; maxX1 = xStart;}
    if (errz1 > maxErrz1) {maxErrz1 = errz1; maxZ1 = zStart;}
    if (errx2 > maxErrx2) {maxErrx2 = errx2; maxX2 = xStart;}
    if (errz2 > maxErrz2) {maxErrz2 = errz2; maxZ2 = zStart;}

    std::cout<<"*************************************"<<std::endl;
    std::cout<<"pt: "<<xStart<<" "<<zStart<<std::endl;
    std::cout<<" Linear: "<<linear_dxdy<<" "<<linear_dzdy<<std::endl;
    std::cout<<"   Spline1: ("<<spline1_dxdy<<" "<<spline1_dzdy<<") "<<std::abs(linear_dxdy-spline1_dxdy)<<" "<<std::abs(linear_dxdy-spline1_dzdy)<<std::endl;
    std::cout<<"   Spline2: ("<<spline2_dxdy<<" "<<spline2_dzdy<<") "<<std::abs(linear_dxdy-spline2_dxdy)<<" "<<std::abs(linear_dxdy-spline2_dzdy)<<std::endl;

    std::cout<<"***** MAX: "<<maxErrx1<<" "<<maxErrz1<<" "<<maxErrx2<<" "<<maxErrz2<<std::endl;
    std::cout<<"        max1 at: "<<maxX1<<" "<<maxZ1<<std::endl;
    std::cout<<"        max2 at: "<<maxX2<<" "<<maxZ2<<std::endl<<std::endl;
    count++;
}

std::pair<double, double> RK4_FLT1(
    double xStart, double yStart, double zStart,
    const std::vector<std::vector<std::vector<double>>>& dxdy,
    const std::vector<std::vector<std::vector<double>>>& dzdy,
     std::vector<double>& xarray,  std::vector<double>& zarray,
    int region,
    const std::vector<std::vector<double>>& dxdy_pm1,
    const std::vector<std::vector<double>>& dzdy_pm1,
    int direction, int nypf1, int nypf2,
    std::ofstream& rk4Out, int iline, int it, bool dumpFiles)
{
    const double twoPi = 2.0 * M_PI;
    double h = 1.0;
    double hh = h / 2.0;
    double h6 = h / 6.0;

    double val0 = dxdy[123][79][101];
    double val1 = dxdy[192][47][200];

    if (dumpFiles)
    {
        writeArray3DToFile(dxdy, "dxdy");
        writeArray3DToFile(dzdy, "dzdy");
    }

    if (direction != 1 && direction != -1)
        throw std::invalid_argument("Direction parameter must be 1 or -1.");

    std::vector<std::vector<double>> dxdyp = Slice(dxdy, yStart);
    double val00 = dxdyp[123][101];
    double val01 = dxdyp[192][200];
    std::vector<std::vector<double>> dzdyp = Slice(dzdy, yStart);
    std::vector<std::vector<double>> dxdyn, dzdyn, dxdyh, dzdyh;

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
        }
    }
    else
    {
        throw std::invalid_argument("Backwards RK4 not supported yet.");
    }

    std::vector<int> _idx0 = { 124-1, 102-1}, _idx1 = {193-1, 201-1};


    double _val00, _val01, _val02, _val10, _val11, _val12;
    _val00 = dxdyn[_idx0[0]][_idx0[1]];
    _val01 = dxdyp[_idx0[0]][_idx0[1]];
    _val02 = dxdyh[_idx0[0]][_idx0[1]];
    _val10 = dxdyn[_idx1[0]][_idx1[1]];
    _val11 = dxdyp[_idx1[0]][_idx1[1]];
    _val12 = dxdyh[_idx1[0]][_idx1[1]];

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

    dxdyp = pad2DEdge(dxdyp);
    dxdyn = pad2DEdge(dxdyn);
    dxdyh = pad2DEdge(dxdyh);
    dzdyp = pad2DEdge(dzdyp);
    dzdyn = pad2DEdge(dzdyn);
    dzdyh = pad2DEdge(dzdyh);

    _val00 = dxdyn[_idx0[0]][_idx0[1]];
    _val01 = dxdyp[_idx0[0]][_idx0[1]];
    _val02 = dxdyh[_idx0[0]][_idx0[1]];
    _val10 = dxdyn[_idx1[0]][_idx1[1]];
    _val11 = dxdyp[_idx1[0]][_idx1[1]];
    _val12 = dxdyh[_idx1[0]][_idx1[1]];

    if (dumpFiles)
    {
        writeArray2DToFile(dxdyp, "dxdyp_");
        writeArray2DToFile(dxdyn, "dxdyn_");
        writeArray2DToFile(dxdyh, "dxdyh_");
        writeArray2DToFile(dzdyp, "dzdyp_");
        writeArray2DToFile(dzdyn, "dzdyn_");
        writeArray2DToFile(dzdyh, "dzdyh_");
    }

    bool doTest = false;
    if (doTest)
    {
        splineTest(xStart, zStart, dxdyp, dzdyp, xarray, zarray);
    }

    bool useSplineInterp = true;

    // Interpolation using the SplineInterpolation class
    //SplineInterpolation splineDxdy(xarray, zarray, dxdyp); //[zStart]);
    //SplineInterpolation splineDzdy(xarray, zarray, dzdyp); //[zStart]);

    //std::cout<<"RK4 step1: "<<xStart<<" "<<zStart<<" xx: "<<xarray.front()<<" "<<xarray.back()<<std::endl;
    //printf("RK4 step1: %12.10e %12.10e\n", xStart, zStart);
    // RK4 Step 1
    double dxdy1, dzdy1;
    if (useSplineInterp)
    {
        SplineInterpolation spline1dx(xarray, zarray, dxdyp), spline1dz(xarray, zarray, dzdyp);
        dxdy1 = spline1dx.evaluate(xStart, zStart);
        dzdy1 = spline1dz.evaluate(xStart, zStart);
    }
    else
    {
        dxdy1 = bilinear_interp2(xarray, zarray, dxdyp, xStart, zStart);
        dzdy1 = bilinear_interp2(xarray, zarray, dzdyp, xStart, zStart);

#if 0
        SplineInterpolation spline1dx(xarray, zarray, dxdyp);
        SplineInterpolation spline1dz(xarray, zarray, dzdyp);
        auto sdzdy1 = spline1dz.evaluate(xStart, zStart);
        auto sdxdy1 = spline1dx.evaluate(xStart, zStart);
        auto dx = std::fabs(dxdy1 - sdxdy1), dz = std::fabs(dzdy1 - sdzdy1);
        //auto dx2 = std::fabs(dxdy1 - tmp1), dz2 = std::fabs(dzdy1 - tmp2);
        if (dx > 1e-9 || dz > 1e-9)
            std::cout<<" Error_sp: "<<dx<<" "<<dz<<std::endl;
        //if (dx2 > 1e-5 || dz2 > 1e-5)
        //    std::cout<<" Error_al: "<<dx2<<" "<<dz2<<std::endl;
#endif
    }
    double x1 = xStart + direction * hh * dxdy1;
    double z1 = zStart + direction * hh * dzdy1;
    double _z1 = z1;
    z1 = double_mod(z1, twoPi);
    rk4Out<<iline<<", "<<double(it)<<", "<<0<<", "<<dxdy1<<", "<<dzdy1<<std::endl;
    //printf("  RES: %12.10e %12.10e\n", dxdy1, dzdy1);
    //printf("  bi-RES: %12.10e %12.10e\n", tmp1, tmp2);
    //std::cout<<"   RES= "<<dxdy1<<" "<<dzdy1<<std::endl;

    // RK4 Step 2
    double dxdy2, dzdy2;
    if (useSplineInterp)
    {
        SplineInterpolation spline2dx(xarray, zarray, dxdyh), spline2dz(xarray, zarray, dzdyh);
        dxdy2 = spline2dx.evaluate(x1, z1);
        dzdy2 = spline2dz.evaluate(x1, z1);
    }
    else
    {
        dxdy2 = bilinear_interp2(xarray, zarray, dxdyh, x1, z1);
        dzdy2 = bilinear_interp2(xarray, zarray, dzdyh, x1, z1);
    }
    double x2 = xStart + direction * hh * dxdy2;
    double z2 = zStart + direction * hh * dzdy2;
    double _z2 = z2;
    z2 = double_mod(z2, twoPi);
    rk4Out<<iline<<", "<<it+0.25<<", "<<0<<", "<<dxdy2<<", "<<dzdy2<<std::endl;

    //std::cout<<"   RES= "<<dxdy2<<" "<<dzdy2<<std::endl;

    // RK4 Step 3
    double dxdy3, dzdy3;
    if (useSplineInterp)
    {
        SplineInterpolation spline3dx(xarray, zarray, dxdyh), spline3dz(xarray, zarray, dzdyh);
        dxdy3 = spline3dx.evaluate(x2, z2);
        dzdy3 = spline3dz.evaluate(x2, z2);
    }
    else
    {
        dxdy3 = bilinear_interp2(xarray, zarray, dxdyh, x2, z2);
        dzdy3 = bilinear_interp2(xarray, zarray, dzdyh, x2, z2);
    }

    double x3 = xStart + direction * dxdy3;
    double z3 = zStart + direction * dzdy3;
    double _z3 = z3;
    z3 = double_mod(z3, twoPi);
    rk4Out<<iline<<", "<<it+0.5<<", "<<0<<", "<<dxdy3<<", "<<dzdy3<<std::endl;

    // RK4 Step 4
    double dxdy4, dzdy4;
    if (useSplineInterp)
    {
        SplineInterpolation spline4dx(xarray, zarray, dxdyn), spline4dz(xarray, zarray, dzdyn);
        dxdy4 = spline4dx.evaluate(x3, z3);
        dzdy4 = spline4dx.evaluate(x3, z3);
    }
    else
    {
        dxdy4 = bilinear_interp2(xarray, zarray, dxdyn, x3, z3);
        dzdy4 = bilinear_interp2(xarray, zarray, dzdyn, x3, z3);
    }

    rk4Out<<iline<<", "<<it+0.75<<", "<<0<<", "<<dxdy4<<", "<<dzdy4<<std::endl;

    // Compute final x and z
    double xEnd = xStart + direction * h6 * (dxdy1 + 2.0 * dxdy2 + 2.0 * dxdy3 + dxdy4);
    double zEnd = zStart + direction * h6 * (dzdy1 + 2.0 * dzdy2 + 2.0 * dzdy3 + dzdy4);
    //zEnd = double_mod(zEnd, twoPi);

    rk4Out<<iline<<", "<<it+0.99<<", "<<0<<", "<<xEnd<<", "<<zEnd<<std::endl;

    //printf("   STEP= %12.10e %12.10e\n", xEnd, zEnd);

    return {xEnd, zEnd};
}
