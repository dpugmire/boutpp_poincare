#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <viskores/cont/ArrayHandle.h>

// Bicubic Spline Interpolation using Eigen
viskores::FloatDefault interp2Spline(const std::vector<viskores::FloatDefault>& x,
                                     const std::vector<viskores::FloatDefault>& y,
                                     const std::vector<std::vector<viskores::FloatDefault>>& Z,
                                     viskores::FloatDefault xi,
                                     viskores::FloatDefault yi);

viskores::FloatDefault alglib_spline(const std::vector<viskores::FloatDefault>& x,
                                     const std::vector<viskores::FloatDefault>& y,
                                     const std::vector<std::vector<viskores::FloatDefault>>& f,
                                     viskores::FloatDefault xi,
                                     viskores::FloatDefault yi);

viskores::FloatDefault interpolate2D_WRONG(const std::vector<viskores::FloatDefault>& xarray,
                                           const std::vector<viskores::FloatDefault>& zarray,
                                           const std::vector<std::vector<viskores::FloatDefault>>& dxdyp,
                                           viskores::FloatDefault x,
                                           viskores::FloatDefault z);

viskores::FloatDefault bilinear_interp2(const std::vector<viskores::FloatDefault>& xarray,
                                        const std::vector<viskores::FloatDefault>& zarray,
                                        const std::vector<std::vector<viskores::FloatDefault>>& data,
                                        viskores::FloatDefault x,
                                        viskores::FloatDefault z);

class SplineInterpolation_OLD
{
public:
  // Constructor
  SplineInterpolation_OLD(const std::vector<viskores::FloatDefault>& x, const std::vector<viskores::FloatDefault>& y)
  {
    if (x.size() < 2 || y.size() < 2)
    {
      throw std::invalid_argument("Input vectors must contain at least two points.");
    }

    this->x = x;
    this->y = y;
    computeCoefficients();
  }

  // Evaluate the spline at a given point
  viskores::FloatDefault evaluate(viskores::FloatDefault t) const
  {
    if (t < x.front())
    {
      std::cout << "Interpolation point is outside the range of input data: " << t << " < " << x.front() << std::endl;
      return x.front();
    }
    if (t > x.back())
    {
      std::cout << "Interpolation point is outside the range of input data: " << t << " > " << x.back() << std::endl;
      return x.back();
    }

    // Find the interval [x[i], x[i+1]] containing t
    auto it = std::lower_bound(x.begin(), x.end(), t);
    size_t i = std::max(static_cast<size_t>(0), static_cast<size_t>(it - x.begin() - 1));

    viskores::FloatDefault h = t - x[i];
    return a[i] + b[i] * h + c[i] * h * h + d[i] * h * h * h;
  }

private:
  std::vector<viskores::FloatDefault> x, y;       // Input data points
  std::vector<viskores::FloatDefault> a, b, c, d; // Spline coefficients

  // Compute spline coefficients
  void computeCoefficients()
  {
    size_t n = x.size() - 1;
    std::vector<viskores::FloatDefault> h(n), alpha(n), l(n + 1), mu(n), z(n + 1);

    a = y;
    for (size_t i = 0; i < n; ++i)
    {
      h[i] = x[i + 1] - x[i];
    }

    for (size_t i = 1; i < n; ++i)
    {
      alpha[i] = (3.0f / h[i] * (a[i + 1] - a[i])) - (3.0f / h[i - 1] * (a[i] - a[i - 1]));
    }

    l[0] = 1.0f;
    mu[0] = 0.0f;
    z[0] = 0.0f;

    for (size_t i = 1; i < n; ++i)
    {
      l[i] = 2.0f * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
      mu[i] = h[i] / l[i];
      z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n] = 1.0f;
    z[n] = 0.0f;
    c.resize(n + 1, 0.0f);
    b.resize(n, 0.0f);
    d.resize(n, 0.0f);

    for (size_t j = n - 1; j < n; --j)
    {
      c[j] = z[j] - mu[j] * c[j + 1];
      b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0f * c[j]) / 3.0f;
      d[j] = (c[j + 1] - c[j]) / (3.0f * h[j]);
    }
  }
};

class SplineInterpolation
{
public:
  // Constructor for 1D interpolation
  SplineInterpolation(const std::vector<viskores::FloatDefault>& x, const std::vector<viskores::FloatDefault>& y);

  // Constructor for 2D interpolation
  SplineInterpolation(const std::vector<viskores::FloatDefault>& x,
                      const std::vector<viskores::FloatDefault>& y,
                      const std::vector<std::vector<viskores::FloatDefault>>& z);

  // Evaluate 1D interpolation
  viskores::FloatDefault evaluate(viskores::FloatDefault x) const;

  // Evaluate 2D interpolation
  viskores::FloatDefault evaluate(viskores::FloatDefault x, viskores::FloatDefault y) const;

private:
  // 1D data
  std::vector<viskores::FloatDefault> x1D, y1D;
  std::vector<viskores::FloatDefault> a1D, b1D, c1D, d1D;

  // 2D data
  std::vector<viskores::FloatDefault> x2D, y2D;
  std::vector<std::vector<viskores::FloatDefault>> z2D;

  // Helper for 1D spline coefficients
  void compute1DSplineCoefficients();

  // Helper for 1D evaluation
  size_t findInterval(const std::vector<viskores::FloatDefault>& x, viskores::FloatDefault val) const;

  // 2D spline coefficients
  std::vector<std::vector<viskores::FloatDefault>> a2D, b2D, c2D, d2D;

  // Helper for 2D spline coefficients
  void compute2DSplineCoefficients();
};


class Spline
{
public:
  // Constructor for 1D interpolation
  Spline(const std::vector<viskores::Vec3f>& pts);

  // Evaluate 1D interpolation
  viskores::Vec3f Evaluate(viskores::FloatDefault t) const;

private:
  std::vector<viskores::FloatDefault> T;
  std::vector<viskores::Vec3f> Points;
  std::vector<viskores::Vec3f> A, B, C, D;
  //std::vector<viskores::FloatDefault> X, Y, Z;
  //std::vector<viskores::FloatDefault> A[3], B[3], C[3], D[3];

  // Helper for 1D spline coefficients
  void ComputeSplineCoefficients();

  // Helper for 1D evaluation
  size_t FindInterval(viskores::FloatDefault val) const;
};
