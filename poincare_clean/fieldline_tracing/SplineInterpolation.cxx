#include "SplineInterpolation.h"
//#include "alglib/src/interpolation.h"
//#include <tinysplinecxx.h>
#include <vector>

//using namespace tinyspline;

viskores::FloatDefault interp2Spline(const std::vector<viskores::FloatDefault>& x,
                                     const std::vector<viskores::FloatDefault>& y,
                                     const std::vector<std::vector<viskores::FloatDefault>>& Z,
                                     viskores::FloatDefault xi,
                                     viskores::FloatDefault yi)
{
#if 0
  int nx = x.size();
  int ny = y.size();

  // Flatten 2D data into a single vector
  std::vector<tsReal> ctrlp;
  for (int j = 0; j < ny; ++j)
  {
    for (int i = 0; i < nx; ++i)
    {
      ctrlp.push_back(x[i]);    // X-coordinate
      ctrlp.push_back(y[j]);    // Y-coordinate
      ctrlp.push_back(Z[j][i]); // Z-value
    }
  }
  BSpline spline(ctrlp.size(), 3, 3);

  // Create a B-spline surface (bicubic: degree 3)
  //tsBSpline spline = tsBSpline::interpolate_cubic_natural(ctrlp, nx, 3);
  //BSpline spline(nx, 3);  // `nx` control points, degree 3 (cubic)

  // Query point
  //DeBoorNet net;
  auto result = spline.evalAll({ xi, yi });
  return result[0];
  //std::vector<tsReal> result = net.result();

  //return result[2]; // Return interpolated Z value
#endif
  return 0.0f;
}

viskores::FloatDefault alglib_spline(const std::vector<viskores::FloatDefault>& x,
                                     const std::vector<viskores::FloatDefault>& y,
                                     const std::vector<std::vector<viskores::FloatDefault>>& f,
                                     viskores::FloatDefault xi,
                                     viskores::FloatDefault yi)
{
  return 0.0;
#if 0
    int nx = x.size();
    int ny = y.size();

    // Flatten 2D data into a single vector
    vector<tsReal> ctrlp;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            ctrlp.push_back(x[i]);   // X-coordinate
            ctrlp.push_back(y[j]);   // Y-coordinate
            ctrlp.push_back(Z[j][i]); // Z-value
        }
    }

    // Create a B-spline surface (bicubic: degree 3)
    tsBSpline spline(nx, ny, 3, TS_CLAMPED);
    spline.setControlPoints(ctrlp);

    // Query point
    vector<tsReal> query = {xi, yi};
    vector<tsReal> result = spline.eval(query).result();

    return result[2]; // Return interpolated Z value
#endif
}

viskores::FloatDefault interpolate2D_WRONG(const std::vector<viskores::FloatDefault>& xarray,
                                           const std::vector<viskores::FloatDefault>& zarray,
                                           const std::vector<std::vector<viskores::FloatDefault>>& dxdyp,
                                           viskores::FloatDefault x,
                                           viskores::FloatDefault z)
{
  // Step 1: Interpolate along z for each x
  size_t nx = xarray.size();
  std::vector<viskores::FloatDefault> tempValues(nx);

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
SplineInterpolation::SplineInterpolation(const std::vector<viskores::FloatDefault>& x, const std::vector<viskores::FloatDefault>& y)
  : x1D(x)
  , y1D(y)
{
  if (x.size() != y.size() || x.empty())
    throw std::invalid_argument("1D input vectors must have the same non-zero size.");

  compute1DSplineCoefficients();
}

// Constructor for 2D interpolation
SplineInterpolation::SplineInterpolation(const std::vector<viskores::FloatDefault>& x,
                                         const std::vector<viskores::FloatDefault>& y,
                                         const std::vector<std::vector<viskores::FloatDefault>>& z)
  : x2D(x)
  , y2D(y)
  , z2D(z)
{
  if (x.empty() || y.empty() || z.empty() || z.size() != x.size() || z[0].size() != y.size())
    throw std::invalid_argument("2D input dimensions must be consistent and non-zero.");

  compute2DSplineCoefficients();
}

void SplineInterpolation::compute1DSplineCoefficients()
{
  size_t n = x1D.size();
  if (n < 2)
    throw std::invalid_argument("At least two points are required for spline interpolation.");

  // Allocate memory for spline coefficients
  a1D = y1D;
  b1D.resize(n - 1);
  c1D.resize(n);
  d1D.resize(n - 1);

  std::vector<viskores::FloatDefault> h(n - 1), alpha(n - 2); // Fix: correct alpha size

  for (size_t i = 0; i < n - 1; ++i)
    h[i] = x1D[i + 1] - x1D[i];

  for (size_t i = 1; i < n - 1; ++i)
  {
    auto tmp = 3.0 * (a1D[i + 1] - a1D[i]);
    auto tmp2 = tmp / h[i];
    tmp = 3.0 * (a1D[i] - a1D[i - 1]);
    tmp2 = tmp / h[i - 1];
    alpha[i - 1] = (3.0 * (a1D[i + 1] - a1D[i]) / h[i]) - (3.0 * (a1D[i] - a1D[i - 1]) / h[i - 1]);
  }

  std::vector<viskores::FloatDefault> l(n), mu(n), z(n);
  l[0] = 1.0;
  mu[0] = z[0] = 0.0;

  for (size_t i = 1; i < n - 1; ++i)
  {
    l[i] = 2.0 * (x1D[i + 1] - x1D[i - 1]) - h[i - 1] * mu[i - 1];
    mu[i] = h[i] / l[i];
    z[i] = (alpha[i - 1] - h[i - 1] * z[i - 1]) / l[i];
  }

  l[n - 1] = 1.0;
  z[n - 1] = c1D[n - 1] = 0.0;

  for (int j = n - 2; j >= 0; --j) // Fix: signed integer
  {
    c1D[j] = z[j] - mu[j] * c1D[j + 1];
    b1D[j] = (a1D[j + 1] - a1D[j]) / h[j] - h[j] * (c1D[j + 1] + 2.0 * c1D[j]) / 3.0;
    d1D[j] = (c1D[j + 1] - c1D[j]) / (3.0 * h[j]);
  }
}

// Evaluate 1D spline
viskores::FloatDefault SplineInterpolation::evaluate(viskores::FloatDefault x) const
{
  size_t i = findInterval(x1D, x);
  viskores::FloatDefault dx = x - x1D[i];
  return a1D[i] + b1D[i] * dx + c1D[i] * dx * dx + d1D[i] * dx * dx * dx;
}

// Find interval for 1D and 2D
size_t SplineInterpolation::findInterval(const std::vector<viskores::FloatDefault>& x, viskores::FloatDefault val) const
{
  if (val < x.front() || val > x.back())
    throw std::out_of_range("Value is outside the interpolation range.");

  //return std::upper_bound(x.begin(), x.end(), val) - x.begin() - 1;
  size_t idx = std::upper_bound(x.begin(), x.end(), val) - x.begin() - 1;
  return (idx >= x.size() - 1) ? x.size() - 2 : idx;
}


void SplineInterpolation::compute2DSplineCoefficients()
{
  size_t nx = x2D.size();
  size_t ny = y2D.size();

  // Resize the 2D coefficient arrays
  a2D.resize(nx, std::vector<viskores::FloatDefault>(ny));
  b2D.resize(nx, std::vector<viskores::FloatDefault>(ny - 1)); // Fix: ny-1
  c2D.resize(nx, std::vector<viskores::FloatDefault>(ny));
  d2D.resize(nx, std::vector<viskores::FloatDefault>(ny - 1)); // Fix: ny-1

  // Compute row-wise splines
  for (size_t i = 0; i < nx; ++i)
  {
    SplineInterpolation rowSpline(y2D, z2D[i]);

    // Assign spline coefficients
    for (size_t j = 0; j < ny; ++j)
    {
      a2D[i][j] = rowSpline.a1D[j];
      c2D[i][j] = rowSpline.c1D[j]; // c1D has size ny
    }
    for (size_t j = 0; j < ny - 1; ++j) // Fix: only iterate up to ny-1
    {
      b2D[i][j] = rowSpline.b1D[j]; // b1D has size ny-1
      d2D[i][j] = rowSpline.d1D[j]; // d1D has size ny-1
    }
  }
}

viskores::FloatDefault SplineInterpolation::evaluate(viskores::FloatDefault x, viskores::FloatDefault y) const
{
  size_t ix = findInterval(x2D, x);
  size_t iy = findInterval(y2D, y);

  viskores::FloatDefault hx = x - x2D[ix];
  viskores::FloatDefault hy = y - y2D[iy];

  // Perform bicubic interpolation using correct weight functions
  viskores::FloatDefault f00 = a2D[ix][iy], f01 = a2D[ix][iy + 1];
  viskores::FloatDefault f10 = a2D[ix + 1][iy], f11 = a2D[ix + 1][iy + 1];

  viskores::FloatDefault b00 = b2D[ix][iy], b01 = b2D[ix][iy + 1];
  viskores::FloatDefault b10 = b2D[ix + 1][iy], b11 = b2D[ix + 1][iy + 1];

  viskores::FloatDefault c00 = c2D[ix][iy], c01 = c2D[ix][iy + 1];
  viskores::FloatDefault c10 = c2D[ix + 1][iy], c11 = c2D[ix + 1][iy + 1];

  viskores::FloatDefault d00 = d2D[ix][iy], d01 = d2D[ix][iy + 1];
  viskores::FloatDefault d10 = d2D[ix + 1][iy], d11 = d2D[ix + 1][iy + 1];

  viskores::FloatDefault fy0 = f00 + b00 * hy + c00 * hy * hy + d00 * hy * hy * hy;
  viskores::FloatDefault fy1 = f10 + b10 * hy + c10 * hy * hy + d10 * hy * hy * hy;

  return fy0 + (fy1 - fy0) * (hx / (x2D[ix + 1] - x2D[ix])); // Linear along x
}
