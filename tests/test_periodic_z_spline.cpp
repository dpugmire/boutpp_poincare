#include "AparFieldModel.h"
#include "Interpolator.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kTwoPi = 2.0 * kPi;

int wrapZIndex(int iz, int nzG)
{
  int wrapped = iz % nzG;
  if (wrapped < 0)
    wrapped += nzG;
  return wrapped;
}

double wrapZ(double z)
{
  double out = std::fmod(z, kTwoPi);
  if (out < 0.0)
    out += kTwoPi;
  if (out >= kTwoPi)
    out -= kTwoPi;
  return out;
}

double periodicLagrange(const std::vector<double> &values, double dz, double z)
{
  const int nzG = static_cast<int>(values.size());
  const double zq = wrapZ(z);
  const int base =
      std::max(0, std::min(nzG - 1, static_cast<int>(std::floor(zq / dz))));

  const int izs[4] = {base - 1, base, base + 1, base + 2};
  const double zc0 = static_cast<double>(izs[0]) * dz;
  const double zc1 = static_cast<double>(izs[1]) * dz;
  const double zc2 = static_cast<double>(izs[2]) * dz;
  const double zc3 = static_cast<double>(izs[3]) * dz;

  const double v0 = values[wrapZIndex(izs[0], nzG)];
  const double v1 = values[wrapZIndex(izs[1], nzG)];
  const double v2 = values[wrapZIndex(izs[2], nzG)];
  const double v3 = values[wrapZIndex(izs[3], nzG)];

  return Interpolator::lagrange4(zq, zc0, zc1, zc2, zc3, v0, v1, v2, v3);
}

bool requireClose(const std::string &label, double actual, double expected)
{
  constexpr double tolerance = 1.0e-10;
  const double error = std::fabs(actual - expected);
  if (error <= tolerance)
    return true;

  std::cerr << label << " mismatch: actual=" << actual << " expected=" << expected
            << " error=" << error << "\n";
  return false;
}

AparData makeData(std::vector<double> &dxdy3dValues,
                  std::vector<double> &dzdy3dValues,
                  std::vector<double> &dxdy2dValues,
                  std::vector<double> &dzdy2dValues)
{
  AparData data;
  data.nx = 5;
  data.ny = 4;
  data.nz = 8;
  data.zperiod = 1;
  data.nzG = 8;
  data.nypf1 = 1;
  data.nypf2 = 2;
  data.zmin = 0.0;
  data.zmax = kTwoPi;
  data.dz_torus = kTwoPi / static_cast<double>(data.nzG);

  data.xarray.resize(data.nx);
  for (int ix = 0; ix < data.nx; ++ix)
    data.xarray[ix] = static_cast<double>(ix + 1);

  dxdy3dValues.resize(data.nzG);
  dzdy3dValues.resize(data.nzG);
  dxdy2dValues.resize(data.nzG);
  dzdy2dValues.resize(data.nzG);
  for (int iz = 0; iz < data.nzG; ++iz)
  {
    dxdy3dValues[iz] = 100.0 + 17.0 * iz + 3.0 * iz * iz;
    dzdy3dValues[iz] = -40.0 + 11.0 * iz - 2.0 * iz * iz;
    dxdy2dValues[iz] = 12.0 - 9.0 * iz + 4.0 * iz * iz;
    dzdy2dValues[iz] = -7.0 + 5.0 * iz + iz * iz * iz;
  }

  data.dxdy.assign(data.nx * data.ny * data.nzG, 0.0);
  data.dzdy.assign(data.nx * data.ny * data.nzG, 0.0);
  for (int ix = 0; ix < data.nx; ++ix)
  {
    for (int iy = 0; iy < data.ny; ++iy)
    {
      for (int iz = 0; iz < data.nzG; ++iz)
      {
        data.dxdy[data.idx3(ix, iy, iz)] = dxdy3dValues[iz];
        data.dzdy[data.idx3(ix, iy, iz)] = dzdy3dValues[iz];
      }
    }
  }

  data.dxdy_p1.assign(data.nx * data.nzG, 0.0);
  data.dzdy_p1.assign(data.nx * data.nzG, 0.0);
  data.dxdy_m1.assign(data.nx * data.nzG, 0.0);
  data.dzdy_m1.assign(data.nx * data.nzG, 0.0);
  for (int ix = 0; ix < data.nx; ++ix)
  {
    for (int iz = 0; iz < data.nzG; ++iz)
    {
      data.dxdy_p1[data.idx_xz(ix, iz)] = dxdy2dValues[iz];
      data.dzdy_p1[data.idx_xz(ix, iz)] = dzdy2dValues[iz];
      data.dxdy_m1[data.idx_xz(ix, iz)] = dxdy2dValues[iz];
      data.dzdy_m1[data.idx_xz(ix, iz)] = dzdy2dValues[iz];
    }
  }

  return data;
}

} // namespace

int main()
{
  std::vector<double> dxdy3dValues;
  std::vector<double> dzdy3dValues;
  std::vector<double> dxdy2dValues;
  std::vector<double> dzdy2dValues;
  const AparData data =
      makeData(dxdy3dValues, dzdy3dValues, dxdy2dValues, dzdy2dValues);
  const AparFieldModel model(data);

  bool ok = true;
  const double dz = data.dz_torus;
  const double zCases[4] = {0.2 * dz, kTwoPi - 0.2 * dz, kTwoPi + 0.2 * dz,
                            -0.2 * dz};

  for (double z : zCases)
  {
    XZDeriv deriv;
    model.evaluateStage({3.0, z}, 2, 1, 1, 0, deriv);
    ok = requireClose("3D dxdy", deriv.dxdy,
                      periodicLagrange(dxdy3dValues, dz, z)) &&
         ok;
    ok = requireClose("3D dzdy", deriv.dzdy,
                      periodicLagrange(dzdy3dValues, dz, z)) &&
         ok;

    model.evaluateStage({3.0, z}, data.nypf2, 0, 1, 2, deriv);
    ok = requireClose("2D twist dxdy", deriv.dxdy,
                      periodicLagrange(dxdy2dValues, dz, z)) &&
         ok;
    ok = requireClose("2D twist dzdy", deriv.dzdy,
                      periodicLagrange(dzdy2dValues, dz, z)) &&
         ok;
  }

  return ok ? 0 : 1;
}
