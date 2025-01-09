#include <vector>
#include "SplineInterpolation.h"

float interpolate2D(
    const std::vector<float>& xarray,
    const std::vector<float>& zarray,
    const std::vector<std::vector<float>>& dxdyp,
    float x, float z)
{
    // Step 1: Interpolate along z for each x
    size_t nx = xarray.size();
    std::vector<float> tempValues(nx);

    for (size_t i = 0; i < nx; ++i)
    {
        SplineInterpolation spline(zarray, dxdyp[i]);
        tempValues[i] = spline.evaluate(z);
    }

    // Step 2: Interpolate along x using the intermediate results
    SplineInterpolation finalSpline(xarray, tempValues);
    return finalSpline.evaluate(x);
}