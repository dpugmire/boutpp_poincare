#include "Spline.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

// If your vector type is different, adjust these accessors.
struct Vec3
{
  double x, y, z;
  Vec3()
    : x(0)
    , y(0)
    , z(0)
  {
  }
  Vec3(double X, double Y, double Z)
    : x(X)
    , y(Y)
    , z(Z)
  {
  }
};
static inline double norm(const viskores::Vec3f& v)
{
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// Helper to read Octave knot and eval files
static bool read_knots(const char* fname, std::vector<double>& t, std::vector<viskores::Vec3f>& P)
{
  std::ifstream fk(fname);
  if (!fk)
    return false;
  double ti, xi, yi, zi;
  while (fk >> ti >> xi >> yi >> zi)
  {
    t.push_back(ti);
    viskores::Vec3f pt(xi, yi, zi);
    P.emplace_back(pt);
  }
  return !t.empty();
}

int main()
{
  // 1) Load knot points used by Octave (for constructing the spline points)
  std::vector<double> t;
  std::vector<viskores::Vec3f> P;
  if (!read_knots("octave_param_knots3d.txt", t, P))
  {
    std::cerr << "Cannot open octave_param_knots3d.txt. Run Octave test first.\n";
    return 1;
  }
  if (t.size() < 2)
  {
    std::cerr << "Not enough control points.\n";
    return 1;
  }

  // 2) Build your 3D parametric spline with implicit parameterization t = 0..N-1
  //    Make sure your class name and API match this call.
  ParametricSpline3D sp3(P); // Constructor only takes points

  // 3) Compare against Octave dense evaluations
  std::ifstream fo("octave_param_spline3d.txt");
  if (!fo)
  {
    std::cerr << "Cannot open octave_param_spline3d.txt\n";
    return 1;
  }

  double max_err = 0.0;
  double max_t = std::numeric_limits<double>::quiet_NaN();
  std::cout << std::fixed << std::setprecision(6);

  double tq, xo, yo, zo;
  while (fo >> tq >> xo >> yo >> zo)
  {
    // Evaluate your spline at the same tq
    // Adjust method name if your API uses Evaluate() instead of eval().
    auto p = sp3.Evaluate(tq); // or sp3.eval(tq)

    // If your return type isn't {x,y,z}, adapt these:
    viskores::Vec3f cc(p[0], p[1], p[2]);

    viskores::Vec3f oc(xo, yo, zo);
    double err = norm(viskores::Vec3f(cc[0] - oc[0], cc[1] - oc[1], cc[2] - oc[2]));
    if (err > max_err)
    {
      max_err = err;
      max_t = tq;
    }

    std::cout << "t=" << tq << "  C++=(" << cc << ")"
              << "  Octave=(" << xo << ", " << yo << ", " << zo << ")"
              << "  |err|=" << err << "\n";
  }

  std::cout << "Max 3D error vs Octave: " << std::setprecision(12) << max_err << " at t=" << max_t << "\n";
  return 0;
}
