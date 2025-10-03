#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <viskores/cont/ArrayHandle.h> // for viskores::Vec3f, FloatDefault


class SplineORIG
{
public:
  // Constructor for 1D interpolation
  SplineORIG(const std::vector<viskores::Vec3f>& pts);

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

class Spline
{
public:
  enum class BoundaryCondition
  {
    NotAKnot, // matches Octave/MATLAB spline
    Natural   // m(0)=m(n-1)=0
  };

  // Construct a cubic spline through 3D points (parameterized t=0..n-1)
  // bc selects boundary conditions; default is NotAKnot (Octave-compatible).
  explicit Spline(const std::vector<viskores::Vec3f>& pts, BoundaryCondition bc = BoundaryCondition::NotAKnot);

  // Evaluate at parameter t. Out-of-range t is clamped to [0, n-1].
  viskores::Vec3f Evaluate(viskores::FloatDefault t) const;

private:
  using T = viskores::FloatDefault;

  // parameter values (0..n-1)
  std::vector<T> Tvals;

  // Per-interval coefficients: s_i(dt) = A[i] + B[i]*dt + C[i]*dt^2 + D[i]*dt^3
  std::vector<viskores::Vec3f> A, B, C, D;

  BoundaryCondition BC = BoundaryCondition::NotAKnot;

  // Build everything
  void BuildSpline(const std::vector<viskores::Vec3f>& pts);

  // Find interval index (with clamping) for t
  std::size_t FindIntervalClamped(T t) const;

  // Small dense LU with partial pivoting
  struct LU
  {
    std::vector<T> A;     // row-major n x n
    std::vector<int> piv; // pivot rows
    std::size_t n = 0;
  };

  // Build/factor matrix for the chosen BC
  LU FactorMatrix(const std::vector<T>& h) const;
  // Solve A x = b in-place (b becomes x)
  void LUSolve(const LU& lu, std::vector<T>& b) const;
};




/// 1-D cubic spline y(x) with NOT-A-KNOT boundary (matches Octave/Matlab `spline`)
class CubicSpline1D
{
public:
  // Build from (x[i], y[i]) data. x must be strictly increasing.
  CubicSpline1D() = default;
  CubicSpline1D(const std::vector<double>& x, const std::vector<double>& y) { build(x, y); }

  void build(const std::vector<double>& x, const std::vector<double>& y);

  // Evaluate at xq. Extrapolates using the end cubic (matches Octave behavior).
  double eval(double xq) const;

  bool empty() const { return x_.empty(); }

private:
  // knots and cubic coefficients on each interval [x[i], x[i+1]):
  // s_i(t) = a[i] + b[i] t + c[i] t^2 + d[i] t^3,  t = (x - x[i])
  std::vector<double> x_, a_, b_, c_, d_;

  // Dense LU to solve the not-a-knot system for second derivatives m
  static void lu_factor(std::vector<double>& A, int n, std::vector<int>& piv);
  static void lu_solve(const std::vector<double>& LU, int n, const std::vector<int>& piv, std::vector<double>& b);

  int find_interval(double xq) const; // returns i in [0, n-2], clamps for extrapolation
};

/// Parametric 3-D cubic spline p(t) using three 1-D splines.
/// Parameterization can be uniform or chord-length.
class ParametricSpline3D
{
public:
  enum class Param
  {
    Uniform,
    ChordLength
  };

  ParametricSpline3D() = default;

  // pts are p0..p{n-1}. Parameterization chosen by 'p'.
  explicit ParametricSpline3D(const std::vector<viskores::Vec3f>& pts, Param p = Param::Uniform) { build(pts, p); }

  void build(const std::vector<viskores::Vec3f>& pts, Param p);

  // Evaluate at parameter t in [t0, tn-1]. Extrapolates at ends.
  viskores::Vec3f Evaluate(double t) const;

  // Get parameter vector (useful if you want to query at same t as Octave)
  const std::vector<double>& params() const { return t_; }

private:
  std::vector<double> t_;
  CubicSpline1D sx_, sy_, sz_;
};
