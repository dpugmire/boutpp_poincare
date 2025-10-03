#include "Spline.h"
#include <limits>

using T = viskores::FloatDefault;


SplineORIG::SplineORIG(const std::vector<viskores::Vec3f>& pts)
  : Points(pts)
{
  auto n = pts.size();

  this->T.resize(n);
  for (std::size_t i = 0; i < n; i++)
    this->T[i] = static_cast<viskores::FloatDefault>(i);


  this->ComputeSplineCoefficients();
}

void SplineORIG::ComputeSplineCoefficients()
{
  size_t n = this->Points.size();
  if (n < 2)
    throw std::invalid_argument("At least two points are required for spline interpolation.");

  std::vector<viskores::Vec3f> alpha(n - 2);
  std::vector<viskores::FloatDefault> h(n - 1);
  //std::vector<viskores::FloatDefault> h(n - 1), alpha(n - 2); // Fix: correct alpha size

  //for (size_t i = 0; i < n - 1; ++i)
  //    h[i] = x1D[i + 1] - x1D[i];
  for (size_t i = 0; i < n - 1; ++i)
    h[i] = this->T[i + 1] - this->T[i];


  this->A = this->Points;
  this->B.resize(n - 1);
  this->C.resize(n);
  this->D.resize(n - 1);

  for (size_t i = 1; i < n - 1; ++i)
  {
    //alpha[i - 1] = (3.0 * (a1D[i + 1] - a1D[i]) / h[i]) - (3.0 * (a1D[i] - a1D[i - 1]) / h[i - 1]);
    auto tmp = 3.0 * (this->A[i + 1] - this->A[i]);
    auto tmp2 = tmp / h[i];
    tmp = 3.0 * (this->A[i] - this->A[i - 1]);
    tmp2 = tmp / h[i];
    alpha[i - 1] = (3.0 * (this->A[i + 1] - this->A[i]) / h[i]) - (3.0 * (this->A[i] - this->A[i - 1]) / h[i - 1]);
  }

  std::vector<viskores::Vec3f> z(n);
  std::vector<viskores::FloatDefault> l(n), mu(n);
  l[0] = 1.0;
  mu[0] = 0.0;
  z[0] = { 0.0, 0.0, 0.0 };

  for (size_t i = 1; i < n - 1; ++i)
  {
    l[i] = 2.0 * (this->T[i + 1] - this->T[i - 1]) - h[i - 1] * mu[i - 1];
    //l[i] = 2.0 * (this->Points[i + 1] - this->Points[i - 1]) - h[i - 1] * mu[i - 1];
    mu[i] = h[i] / l[i];
    z[i] = (alpha[i - 1] - h[i - 1] * z[i - 1]) / l[i];
  }


  l[n - 1] = 1.0;
  z[n - 1] = this->C[n - 1] = { 0.0, 0.0, 0.0 };

  for (int j = n - 2; j >= 0; --j) // Fix: signed integer
  {
    this->C[j] = z[j] - mu[j] * this->C[j + 1];
    this->B[j] = (this->A[j + 1] - this->A[j]) / h[j] - h[j] * (this->C[j + 1] + 2.0 * this->C[j]) / 3.0;
    this->D[j] = (this->C[j + 1] - this->C[j]) / (3.0 * h[j]);
  }
}

viskores::Vec3f SplineORIG::Evaluate(viskores::FloatDefault t) const
{
  size_t i = this->FindInterval(t);
  viskores::FloatDefault dx = t - this->T[i];
  return this->A[i] + this->B[i] * dx + this->C[i] * dx * dx + this->D[i] * dx * dx * dx;
}


size_t SplineORIG::FindInterval(viskores::FloatDefault val) const
{
  if (val < this->T.front() || val > this->T.back())
    throw std::out_of_range("Value is outside the interpolation range.");

  //return std::upper_bound(x.begin(), x.end(), val) - x.begin() - 1;
  size_t idx = std::upper_bound(this->T.begin(), this->T.end(), val) - this->T.begin() - 1;
  return (idx >= this->T.size() - 1) ? this->T.size() - 2 : idx;
}


Spline::Spline(const std::vector<viskores::Vec3f>& pts, BoundaryCondition bc)
  : BC(bc)
{
  if (pts.size() < 2)
    throw std::invalid_argument("Spline needs at least two points.");
  this->BuildSpline(pts);
}

void Spline::BuildSpline(const std::vector<viskores::Vec3f>& pts)
{
  const std::size_t n = pts.size();

  // Parameterization: uniform t = 0..n-1
  this->Tvals.resize(n);
  for (std::size_t i = 0; i < n; ++i)
    this->Tvals[i] = static_cast<T>(i);

  // Interval lengths
  if (n == 2)
  {
    // Degenerate: single interval => linear interpolation
    this->A.resize(n);
    this->B.resize(1);
    this->C.assign(n, viskores::Vec3f{ 0, 0, 0 });
    this->D.assign(1, viskores::Vec3f{ 0, 0, 0 });

    this->A[0] = pts[0];
    this->A[1] = pts[1]; // unused by eval but keep size consistent
    const T h = T(1);
    this->B[0] = (pts[1] - pts[0]) / h;
    return;
  }

  std::vector<T> h(n - 1);
  for (std::size_t i = 0; i + 1 < n; ++i)
    h[i] = this->Tvals[i + 1] - this->Tvals[i]; // = 1 with uniform t

  // Factor the n×n system once for the chosen BC
  LU lu = this->FactorMatrix(h);

  // RHS builder for one component (second derivatives m[])
  auto solve_m = [&](int comp)
  {
    std::vector<T> rhs(n, T(0));

    // Interior equations (common to both BCs):
    // h[i-1] m[i-1] + 2(h[i-1]+h[i]) m[i] + h[i] m[i+1] =
    //   6( (y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1] ), i=1..n-2
    for (std::size_t i = 1; i + 1 < n; ++i)
    {
      const T y_im1 = pts[i - 1][comp];
      const T y_i = pts[i][comp];
      const T y_ip1 = pts[i + 1][comp];
      rhs[i] = 6.0 * ((y_ip1 - y_i) / h[i] - (y_i - y_im1) / h[i - 1]);
    }

    // Boundary rows depend on BC:
    //  - Not-a-knot: both rows are homogeneous (rhs already zero).
    //  - Natural: m0=0, mn-1=0 → rhs already zero.
    // So no action needed here; just solve.
    this->LUSolve(lu, rhs);
    return rhs; // this is m (second derivatives)
  };

  // Solve for m_x, m_y, m_z
  std::vector<T> mx = solve_m(0);
  std::vector<T> my = solve_m(1);
  std::vector<T> mz = solve_m(2);

  // Convert to piecewise cubic coefficients
  // s_i(dt) = A[i] + B[i] dt + C[i] dt^2 + D[i] dt^3, dt ∈ [0, h[i]]
  this->A.resize(n);
  this->B.resize(n - 1);
  this->C.resize(n);
  this->D.resize(n - 1);

  for (std::size_t i = 0; i < n; ++i)
  {
    this->A[i] = pts[i];
    this->C[i][0] = mx[i] * T(0.5);
    this->C[i][1] = my[i] * T(0.5);
    this->C[i][2] = mz[i] * T(0.5);
  }

  for (std::size_t i = 0; i + 1 < n; ++i)
  {
    const T hi = h[i];
    for (int c = 0; c < 3; ++c)
    {
      const T yi = pts[i][c];
      const T yip1 = pts[i + 1][c];
      const T mi = (c == 0 ? mx[i] : (c == 1 ? my[i] : mz[i]));
      const T mip1 = (c == 0 ? mx[i + 1] : (c == 1 ? my[i + 1] : mz[i + 1]));

      const T Bi = (yip1 - yi) / hi - ((2.0 * mi + mip1) * hi / 6.0);
      const T Di = (mip1 - mi) / (6.0 * hi);

      this->B[i][c] = Bi;
      this->D[i][c] = Di;
    }
  }
}

std::size_t Spline::FindIntervalClamped(T t) const
{
  if (t <= Tvals.front())
    return 0;
  if (t >= Tvals.back())
    return Tvals.size() - 2;

  auto it = std::upper_bound(Tvals.begin(), Tvals.end(), t);
  std::size_t idx = static_cast<std::size_t>(std::distance(Tvals.begin(), it)) - 1;
  if (idx >= Tvals.size() - 1)
    idx = Tvals.size() - 2;
  return idx;
}

viskores::Vec3f Spline::Evaluate(T t) const
{
  // Clamp to parameter domain
  if (t < Tvals.front())
    t = Tvals.front();
  if (t > Tvals.back())
    t = Tvals.back();

  const std::size_t i = this->FindIntervalClamped(t);
  const T dt = t - this->Tvals[i];
  return this->A[i] + this->B[i] * dt + this->C[i] * (dt * dt) + this->D[i] * (dt * dt * dt);
}

/* ---------------- Matrix construction + LU ---------------- */

Spline::LU Spline::FactorMatrix(const std::vector<T>& h) const
{
  const std::size_t n = h.size() + 1;
  LU lu;
  lu.n = n;
  lu.A.assign(n * n, T(0));
  lu.piv.resize(n);
  for (std::size_t i = 0; i < n; ++i)
    lu.piv[i] = int(i);

  auto Aat = [&](std::size_t r, std::size_t c) -> T& { return lu.A[r * n + c]; };

  if (this->BC == BoundaryCondition::NotAKnot)
  {
    // Row 0: h1*m0 - (h0+h1)*m1 + h0*m2 = 0
    if (n >= 3)
    {
      Aat(0, 0) = h[1];
      Aat(0, 1) = -(h[0] + h[1]);
      Aat(0, 2) = h[0];
    }
    else
    {
      // n==2 fallback to identity (shouldn’t get here; handled earlier)
      Aat(0, 0) = 1.0;
    }

    // Interior rows: i=1..n-2
    for (std::size_t i = 1; i + 1 < n; ++i)
    {
      Aat(i, i - 1) = h[i - 1];
      Aat(i, i) = 2.0 * (h[i - 1] + h[i]);
      Aat(i, i + 1) = h[i];
    }

    // Last row: h_{n-2}*m_{n-3} - (h_{n-3}+h_{n-2})*m_{n-2} + h_{n-3}*m_{n-1} = 0
    if (n >= 3)
    {
      Aat(n - 1, n - 3) = h[n - 2];
      Aat(n - 1, n - 2) = -(h[n - 3] + h[n - 2]);
      Aat(n - 1, n - 1) = h[n - 3];
    }
    else
    {
      Aat(n - 1, n - 1) = 1.0;
    }
  }
  else
  {
    // Natural spline
    // Row 0: m0 = 0
    Aat(0, 0) = 1.0;

    // Interior rows: i=1..n-2
    for (std::size_t i = 1; i + 1 < n; ++i)
    {
      Aat(i, i - 1) = h[i - 1];
      Aat(i, i) = 2.0 * (h[i - 1] + h[i]);
      Aat(i, i + 1) = h[i];
    }

    // Last row: mn-1 = 0
    Aat(n - 1, n - 1) = 1.0;
  }

  // LU factorization with partial pivoting
  for (std::size_t k = 0; k < n; ++k)
  {
    // Pivot
    std::size_t pivrow = k;
    T maxabs = std::abs(Aat(k, k));
    for (std::size_t r = k + 1; r < n; ++r)
    {
      T v = std::abs(Aat(r, k));
      if (v > maxabs)
      {
        maxabs = v;
        pivrow = r;
      }
    }
    if (maxabs <= std::numeric_limits<T>::epsilon())
      throw std::runtime_error("LU factorization failed (singular matrix).");

    if (pivrow != k)
    {
      for (std::size_t c = 0; c < n; ++c)
        std::swap(Aat(k, c), Aat(pivrow, c));
      std::swap(lu.piv[k], lu.piv[pivrow]);
    }

    // Eliminate
    for (std::size_t r = k + 1; r < n; ++r)
    {
      Aat(r, k) /= Aat(k, k);
      T mult = Aat(r, k);
      for (std::size_t c = k + 1; c < n; ++c)
        Aat(r, c) -= mult * Aat(k, c);
    }
  }

  return lu;
}

void Spline::LUSolve(const LU& lu, std::vector<T>& b) const
{
  const std::size_t n = lu.n;
  auto Aat = [&](std::size_t r, std::size_t c) -> const T& { return lu.A[r * n + c]; };

  // Apply row pivots to b
  for (std::size_t i = 0; i < n; ++i)
  {
    int pr = lu.piv[i];
    if (pr != int(i))
      std::swap(b[i], b[pr]);
  }

  // Forward solve Ly = Pb
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < i; ++j)
      b[i] -= Aat(i, j) * b[j];
    // diag(L)=1
  }

  // Backward solve Ux = y
  for (int i = int(n) - 1; i >= 0; --i)
  {
    for (std::size_t j = i + 1; j < n; ++j)
      b[std::size_t(i)] -= Aat(std::size_t(i), j) * b[j];
    b[std::size_t(i)] /= Aat(std::size_t(i), std::size_t(i));
  }
}


#include "Spline.h"
#include <cassert>
#include <limits>

/* ---------------- CubicSpline1D ---------------- */

void CubicSpline1D::build(const std::vector<double>& x, const std::vector<double>& y)
{
  const int n = static_cast<int>(x.size());
  if (n < 2 || y.size() != x.size())
    throw std::invalid_argument("CubicSpline1D: need >=2 points and matching sizes.");
  for (int i = 1; i < n; ++i)
    if (!(x[i] > x[i - 1]))
      throw std::invalid_argument("CubicSpline1D: x must be strictly increasing.");

  x_ = x;
  a_ = y; // a[i] = y[i]
  std::vector<double> h(n - 1);
  for (int i = 0; i < n - 1; ++i)
    h[i] = x_[i + 1] - x_[i];

  // Build dense system A * m = rhs for second derivatives m[0..n-1]
  const int N = n;
  std::vector<double> A(N * N, 0.0), rhs(N, 0.0);
  auto at = [&](int r, int c) -> double& { return A[r * N + c]; };

  // Not-a-knot boundary rows
  // Row 0: h1*m0 - (h0+h1)*m1 + h0*m2 = 0
  if (n >= 3)
  {
    at(0, 0) = h[1];
    at(0, 1) = -(h[0] + h[1]);
    at(0, 2) = h[0];
  }
  else
  {
    // n == 2: degenerate to identity; m[0]=m[1]=0
    at(0, 0) = 1.0;
    at(1, 1) = 1.0;
  }

  // Interior rows i=1..n-2:
  for (int i = 1; i <= n - 2; ++i)
  {
    at(i, i - 1) = h[i - 1];
    at(i, i) = 2.0 * (h[i - 1] + (i < n - 1 ? h[i] : 0.0));
    if (i < n - 1)
      at(i, i + 1) = h[i];
    if (i >= 1 && i <= n - 2)
    {
      rhs[i] = 6.0 * ((a_[i + 1] - a_[i]) / h[i] - (a_[i] - a_[i - 1]) / h[i - 1]);
    }
  }

  // Last row: h_{n-2}*m_{n-3} - (h_{n-3}+h_{n-2})*m_{n-2} + h_{n-3}*m_{n-1} = 0
  if (n >= 3)
  {
    at(n - 1, n - 3) = h[n - 2];
    at(n - 1, n - 2) = -(h[n - 3] + h[n - 2]);
    at(n - 1, n - 1) = h[n - 3];
  }

  // Factor + solve
  std::vector<int> piv(N);
  lu_factor(A, N, piv);
  lu_solve(A, N, piv, rhs);
  const std::vector<double>& m = rhs; // reuse storage

  // Convert to cubic coefficients on each interval
  b_.assign(n - 1, 0.0);
  c_.assign(n, 0.0);
  d_.assign(n - 1, 0.0);

  for (int i = 0; i < n; ++i)
    c_[i] = 0.5 * m[i];

  for (int i = 0; i < n - 1; ++i)
  {
    const double hi = h[i];
    const double yi = a_[i];
    const double yi1 = a_[i + 1];
    const double mi = m[i];
    const double mi1 = m[i + 1];

    b_[i] = (yi1 - yi) / hi - (2.0 * mi + mi1) * hi / 6.0;
    d_[i] = (mi1 - mi) / (6.0 * hi);
  }
}

int CubicSpline1D::find_interval(double xq) const
{
  // clamp to [0..n-2] so we extrapolate with end cubic
  if (xq <= x_.front())
    return 0;
  if (xq >= x_.back())
    return static_cast<int>(x_.size()) - 2;

  auto it = std::upper_bound(x_.begin(), x_.end(), xq);
  int i = static_cast<int>(std::distance(x_.begin(), it)) - 1;
  if (i < 0)
    i = 0;
  if (i > static_cast<int>(x_.size()) - 2)
    i = static_cast<int>(x_.size()) - 2;
  return i;
}

double CubicSpline1D::eval(double xq) const
{
  if (x_.empty())
    throw std::runtime_error("CubicSpline1D: not built.");
  const int i = find_interval(xq);
  const double t = xq - x_[i];
  return a_[i] + b_[i] * t + c_[i] * t * t + d_[i] * t * t * t;
}

/* ---- simple dense LU with partial pivoting ---- */

void CubicSpline1D::lu_factor(std::vector<double>& A, int n, std::vector<int>& piv)
{
  auto at = [&](int r, int c) -> double& { return A[r * n + c]; };
  for (int i = 0; i < n; ++i)
    piv[i] = i;

  for (int k = 0; k < n; ++k)
  {
    // pivot
    int pivrow = k;
    double maxabs = std::abs(at(k, k));
    for (int r = k + 1; r < n; ++r)
    {
      double v = std::abs(at(r, k));
      if (v > maxabs)
      {
        maxabs = v;
        pivrow = r;
      }
    }
    if (maxabs <= std::numeric_limits<double>::epsilon())
      throw std::runtime_error("LU: singular matrix.");

    if (pivrow != k)
    {
      for (int c = 0; c < n; ++c)
        std::swap(at(k, c), at(pivrow, c));
      std::swap(piv[k], piv[pivrow]);
    }

    // eliminate
    for (int r = k + 1; r < n; ++r)
    {
      at(r, k) /= at(k, k);
      const double m = at(r, k);
      for (int c = k + 1; c < n; ++c)
        at(r, c) -= m * at(k, c);
    }
  }
}

void CubicSpline1D::lu_solve(const std::vector<double>& LU, int n, const std::vector<int>& piv, std::vector<double>& b)
{
  auto at = [&](int r, int c) -> double { return LU[r * n + c]; };
  // apply pivots to b
  std::vector<double> bb = b;
  for (int i = 0; i < n; ++i)
    b[i] = bb[piv[i]];

  // forward solve (L y = b)
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < i; ++j)
      b[i] -= at(i, j) * b[j];
    // diag(L) = 1
  }
  // back solve (U x = y)
  for (int i = n - 1; i >= 0; --i)
  {
    for (int j = i + 1; j < n; ++j)
      b[i] -= at(i, j) * b[j];
    b[i] /= at(i, i);
  }
}

/* ---------------- ParametricSpline3D ---------------- */

static std::vector<double> make_params(const std::vector<viskores::Vec3f>& p, ParametricSpline3D::Param kind)
{
  const int n = static_cast<int>(p.size());
  std::vector<double> t(n, 0.0);
  if (n == 0)
    return t;

  if (kind == ParametricSpline3D::Param::Uniform)
  {
    for (int i = 0; i < n; ++i)
      t[i] = static_cast<double>(i);
    return t;
  }

  // chord length
  t[0] = 0.0;
  for (int i = 1; i < n; ++i)
  {
    const double dx = p[i][0] - p[i - 1][0];
    const double dy = p[i][1] - p[i - 1][1];
    const double dz = p[i][2] - p[i - 1][2];
    t[i] = t[i - 1] + std::sqrt(dx * dx + dy * dy + dz * dz);
  }
  // normalize to 0..(n-1) for a stable scale
  if (t.back() > 0)
  {
    const double scale = (n - 1) / t.back();
    for (double& v : t)
      v *= scale;
  }
  return t;
}

void ParametricSpline3D::build(const std::vector<viskores::Vec3f>& pts, Param p)
{
  if (pts.size() < 2)
    throw std::invalid_argument("ParametricSpline3D: need >=2 points.");
  t_ = make_params(pts, p);

  std::vector<double> x(pts.size()), y(pts.size()), z(pts.size());
  for (std::size_t i = 0; i < pts.size(); ++i)
  {
    x[i] = t_[i];
    y[i] = pts[i][0];
    z[i] = pts[i][1];
  }
  // we’ll map: t -> (X(t), Y(t), Z(t)) by three 1-D splines
  sx_.build(t_, y);
  sy_.build(t_, z);
  // if you have real Z for 3rd component, add it; otherwise keep 0s:
  std::vector<double> zz(pts.size(), 0.0);
  for (std::size_t i = 0; i < pts.size(); ++i)
    zz[i] = pts[i][2];
  sz_.build(t_, zz);
}

viskores::Vec3f ParametricSpline3D::Evaluate(double t) const
{
  viskores::Vec3f p;
  p[0] = sx_.eval(t);
  p[1] = sy_.eval(t);
  p[2] = sz_.eval(t);
  return p;
}
