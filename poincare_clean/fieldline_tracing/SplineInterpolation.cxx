#include <vector>
#include "SplineInterpolation.h"

double interpolate2D(
    const std::vector<double>& xarray,
    const std::vector<double>& zarray,
    const std::vector<std::vector<double>>& dxdyp,
    double x, double z)
{
    // Step 1: Interpolate along z for each x
    size_t nx = xarray.size();
    std::vector<double> tempValues(nx);

    for (size_t i = 0; i < nx; ++i)
    {
        SplineInterpolation_OLD spline(zarray, dxdyp[i]);
        tempValues[i] = spline.evaluate(z);
    }

    // Step 2: Interpolate along x using the intermediate results
    SplineInterpolation_OLD finalSpline(xarray, tempValues);
    return finalSpline.evaluate(x);
}


// Constructor for 1D interpolation
SplineInterpolation::SplineInterpolation(const std::vector<double>& x, const std::vector<double>& y) : x1D(x), y1D(y)
{
    if (x.size() != y.size() || x.empty())
        throw std::invalid_argument("1D input vectors must have the same non-zero size.");

    compute1DSplineCoefficients();
}

// Constructor for 2D interpolation
SplineInterpolation::SplineInterpolation(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& z)
    : x2D(x), y2D(y), z2D(z)
{
    if (x.empty() || y.empty() || z.empty() || z.size() != x.size() || z[0].size() != y.size())
        throw std::invalid_argument("2D input dimensions must be consistent and non-zero.");

    compute2DSplineCoefficients();
}

// Compute 1D spline coefficients
void SplineInterpolation::compute1DSplineCoefficients()
{
    size_t n = x1D.size();
    a1D = y1D;
    b1D.resize(n - 1);
    c1D.resize(n);
    d1D.resize(n - 1);

    std::vector<double> h(n - 1), alpha(n - 1);
    for (size_t i = 0; i < n - 1; ++i)
        h[i] = x1D[i + 1] - x1D[i];

    for (size_t i = 1; i < n - 1; ++i)
        alpha[i] = (3.0 / h[i] * (a1D[i + 1] - a1D[i])) - (3.0 / h[i - 1] * (a1D[i] - a1D[i - 1]));

    std::vector<double> l(n), mu(n), z(n);
    l[0] = 1.0;
    mu[0] = z[0] = 0.0;

    for (size_t i = 1; i < n - 1; ++i)
    {
        l[i] = 2.0 * (x1D[i + 1] - x1D[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1.0;
    z[n - 1] = c1D[n - 1] = 0.0;

    for (size_t j = n - 2; j != static_cast<size_t>(-1); --j)
    {
        c1D[j] = z[j] - mu[j] * c1D[j + 1];
        b1D[j] = (a1D[j + 1] - a1D[j]) / h[j] - h[j] * (c1D[j + 1] + 2.0 * c1D[j]) / 3.0;
        d1D[j] = (c1D[j + 1] - c1D[j]) / (3.0 * h[j]);
    }
}

// Evaluate 1D spline
double SplineInterpolation::evaluate(double x) const
{
    size_t i = findInterval(x1D, x);
    double dx = x - x1D[i];
    return a1D[i] + b1D[i] * dx + c1D[i] * dx * dx + d1D[i] * dx * dx * dx;
}

// Find interval for 1D and 2D
size_t SplineInterpolation::findInterval(const std::vector<double>& x, double val) const
{
    if (val < x.front() || val > x.back())
        throw std::out_of_range("Value is outside the interpolation range.");

    return std::upper_bound(x.begin(), x.end(), val) - x.begin() - 1;
}


void SplineInterpolation::compute2DSplineCoefficients()
{
    size_t nx = x2D.size();
    size_t ny = y2D.size();

    // Resize the 2D coefficient arrays
    a2D.resize(nx, std::vector<double>(ny));
    b2D.resize(nx, std::vector<double>(ny));
    c2D.resize(nx, std::vector<double>(ny));
    d2D.resize(nx, std::vector<double>(ny));

    // Compute row-wise splines
    for (size_t i = 0; i < nx; ++i)
    {
        SplineInterpolation rowSpline(y2D, z2D[i]);

        // Assign each component of the spline coefficients
        for (size_t j = 0; j < ny; ++j)
        {
            a2D[i][j] = rowSpline.a1D[j];
            b2D[i][j] = rowSpline.b1D[j];
            c2D[i][j] = rowSpline.c1D[j];
            d2D[i][j] = rowSpline.d1D[j];
        }
    }
}


// Evaluate 2D spline
double SplineInterpolation::evaluate(double x, double y) const
{
    size_t ix = findInterval(x2D, x);
    size_t iy = findInterval(y2D, y);

    double hx = x - x2D[ix];
    double hy = y - y2D[iy];

    return a2D[ix][iy] + b2D[ix][iy] * hy + c2D[ix][iy] * hy * hy + d2D[ix][iy] * hy * hy * hy;
}

