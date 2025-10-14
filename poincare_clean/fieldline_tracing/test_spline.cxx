#include "Spline.h"
#include <fstream>
#include <iostream>
#include <vector>

int main()
{
  // Same knots and values as Octave used to write octave_spline_data.txt
  std::vector<double> x = { 0, 1, 2, 3, 4 };
  std::vector<double> y = { 0, 1, 0, -1, 0 };

  CubicSpline1D sp(x, y);

  std::ifstream fin("octave_spline_data.txt");
  double xo, yo;
  double max_err = 0;
  while (fin >> xo >> yo)
  {
    double yc = sp.eval(xo);
    double err = std::abs(yc - yo);
    max_err = std::max(max_err, err);
    std::cout << xo << " " << yc << " (octave=" << yo << ", err=" << err << ")\n";
  }
  std::cout << "Max error vs Octave: " << max_err << "\n";
}
