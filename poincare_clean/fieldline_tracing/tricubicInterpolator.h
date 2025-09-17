
#include <algorithm>
#include <cmath>
#include <vector>
#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/ArrayHandleCartesianProduct.h>
#include <viskores/cont/ArrayHandleUniformPointCoordinates.h>

// 1D cubic‐convolution (Catmull–Rom) kernel
// given four samples p0,p1,p2,p3 and relative t in [0,1]
inline viskores::FloatDefault cubicInterpolate(viskores::FloatDefault p0,
                                               viskores::FloatDefault p1,
                                               viskores::FloatDefault p2,
                                               viskores::FloatDefault p3,
                                               viskores::FloatDefault t)
{
  // Catmull–Rom basis: a = –0.5
  const viskores::FloatDefault a = -0.5;
  viskores::FloatDefault t2 = t * t;
  viskores::FloatDefault t3 = t2 * t;

  viskores::FloatDefault m0 = (p2 - p0) * 0.5;
  viskores::FloatDefault m1 = (p3 - p1) * 0.5;
  viskores::FloatDefault d0 = p1;
  viskores::FloatDefault d1 = p2;

  // Hermite form: h00, h10, h01, h11
  viskores::FloatDefault h00 = 2 * t3 - 3 * t2 + 1;
  viskores::FloatDefault h10 = t3 - 2 * t2 + t;
  viskores::FloatDefault h01 = -2 * t3 + 3 * t2;
  viskores::FloatDefault h11 = t3 - t2;

  return h00 * d0 + h10 * m0 + h01 * d1 + h11 * m1;
}

// tricubic interpolation at (x,y,z) in index‐space
// input data is densely packed: x fastest, then y, then z
// data.size() == nx*ny*nz
inline viskores::FloatDefault tricubicSample(const viskores::cont::ArrayHandle<viskores::FloatDefault>& data,
                                             const viskores::Id3& dims,
                                             const viskores::Vec3f& pt)
{
  auto portal = data.ReadPortal();
  viskores::Id nx = dims[0], ny = dims[1], nz = dims[2];
  viskores::FloatDefault x = pt[0], y = pt[1], z = pt[2];

  // clamp helper for boundaries
  auto clamp = [&](viskores::Id i, viskores::Id n) { return (i < 0 ? 0 : i >= n ? (n - 1) : i); };

  // base integer coords
  viskores::Id ix = static_cast<viskores::Id>(std::floor(x));
  viskores::Id iy = static_cast<viskores::Id>(std::floor(y));
  viskores::Id iz = static_cast<viskores::Id>(std::floor(z));

  // fractional offsets
  viskores::FloatDefault tx = x - ix;
  viskores::FloatDefault ty = y - iy;
  viskores::FloatDefault tz = z - iz;

  // 1D buffers for the 4×4×4 neighbourhood and intermediate results
  viskores::FloatDefault P[4 * 4 * 4]; // holds P[kk][jj][ii] flattened
  viskores::FloatDefault C[4 * 4];     // holds C[kk][jj] flattened
  viskores::FloatDefault D[4];         // holds D[kk]

  // 1) Gather 4×4×4 neighborhood into P
  //    P[kk][jj][ii] -> P[(kk*4 + jj)*4 + ii]
  //    data[z0][y0][x0] -> data[(z0*ny + y0)*nx + x0]
  for (viskores::Id kk = 0; kk < 4; ++kk)
  {
    viskores::Id z0 = clamp(iz - 1 + kk, nz);
    for (viskores::Id jj = 0; jj < 4; ++jj)
    {
      viskores::Id y0 = clamp(iy - 1 + jj, ny);
      for (viskores::Id ii = 0; ii < 4; ++ii)
      {
        viskores::Id x0 = clamp(ix - 1 + ii, nx);
        // flatten 3D (kk,jj,ii) to 1D:
        viskores::Id pIndex = (kk * 4 + jj) * 4 + ii;
        // flatten volume coords: (x0,y0,z0) -> 1D index
        viskores::Id dIndex = (z0 * ny + y0) * nx + x0;
        P[pIndex] = portal.Get(dIndex);
      }
    }
  }

  // 2) Interpolate in X for each (kk, jj) → C[kk][jj]
  //    C[kk][jj] -> C[kk*4 + jj]
  for (viskores::Id kk = 0; kk < 4; ++kk)
  {
    for (viskores::Id jj = 0; jj < 4; ++jj)
    {
      // base offset for P row
      viskores::Id baseP = (kk * 4 + jj) * 4;
      viskores::Id cIndex = kk * 4 + jj;
      C[cIndex] = cubicInterpolate(P[baseP + 0], P[baseP + 1], P[baseP + 2], P[baseP + 3], tx);
    }
  }

  // 3) Interpolate in Y for each kk → D[kk]
  //    C[kk][0..3] -> C[kk*4 + 0..3]
  for (viskores::Id kk = 0; kk < 4; ++kk)
  {
    D[kk] = cubicInterpolate(C[kk * 4 + 0], C[kk * 4 + 1], C[kk * 4 + 2], C[kk * 4 + 3], ty);
  }

  // 4) Interpolate in Z across D[0..3]
  return cubicInterpolate(D[0], D[1], D[2], D[3], tz);
}

inline viskores::Id findIndex(const viskores::cont::ArrayHandle<viskores::FloatDefault>& _coords, viskores::FloatDefault val)
{
  auto coords = _coords.ReadPortal();
  viskores::Id N = coords.GetNumberOfValues();
  // 1) Binary search for the largest index i with coords[i] <= val
  viskores::Id left = 0;
  viskores::Id right = N - 1;
  while (left <= right)
  {
    viskores::Id mid = left + (right - left) / 2;
    if (coords.Get(mid) <= val)
    {
      // mid is still ≤ val, so it might be our i
      left = mid + 1;
    }
    else
    {
      // coords[mid] > val, so the index we want is below mid
      right = mid - 1;
    }
  }
  // when loop ends, `right` is the last index where coords[right] <= val
  viskores::Id i = right;

  // 2) Clamp i into [1, N-3]
  if (i < 1)
    i = 1;
  else if (i > N - 3)
    i = N - 3;

  return i;
}

inline double cubicInterpolateNonUniform(viskores::FloatDefault x0,
                                         viskores::FloatDefault x1,
                                         viskores::FloatDefault x2,
                                         viskores::FloatDefault x3,
                                         viskores::FloatDefault p0,
                                         viskores::FloatDefault p1,
                                         viskores::FloatDefault p2,
                                         viskores::FloatDefault p3,
                                         viskores::FloatDefault x)
{
  // 1) Compute interval lengths
  viskores::FloatDefault h0 = x1 - x0;
  viskores::FloatDefault h1 = x2 - x1;
  viskores::FloatDefault h2 = x3 - x2;
  if (h0 <= 0 || h1 <= 0 || h2 <= 0)
    throw std::runtime_error("cubicInterpolateNonUniform: coordinates must be strictly increasing");

  // 2) Compute right‐hand sides for second‐derivative system
  viskores::FloatDefault rhs1 = 6.0 * ((p2 - p1) / h1 - (p1 - p0) / h0);
  viskores::FloatDefault rhs2 = 6.0 * ((p3 - p2) / h2 - (p2 - p1) / h1);

  // 3) Build and solve the 2×2 system:
  //     [2(h0+h1)   h1      ][d2_1] = [rhs1]
  //     [  h1     2(h1+h2)  ][d2_2]   [rhs2]
  viskores::FloatDefault a11 = 2.0 * (h0 + h1);
  viskores::FloatDefault a12 = h1;
  viskores::FloatDefault a21 = h1;
  viskores::FloatDefault a22 = 2.0 * (h1 + h2);
  viskores::FloatDefault det = a11 * a22 - a12 * a21;
  if (det == 0.0)
    throw std::runtime_error("cubicInterpolateNonUniform: degenerate knot spacing");
  viskores::FloatDefault d2_1 = (rhs1 * a22 - a12 * rhs2) / det;
  viskores::FloatDefault d2_2 = (a11 * rhs2 - rhs1 * a21) / det;

  // 4) Map x into local parameter t ∈ [0,1] on [x1,x2]
  viskores::FloatDefault t = (x - x1) / h1;

  // 5) Hermite form of the natural cubic on [x1, x2]
  viskores::FloatDefault A = 1.0 - t;
  viskores::FloatDefault B = t;
  viskores::FloatDefault h1_sq = h1 * h1;
  viskores::FloatDefault term1 = (A * A * A - A) * (h1_sq / 6.0) * d2_1;
  viskores::FloatDefault term2 = (B * B * B - B) * (h1_sq / 6.0) * d2_2;

  // 6) Combine the linear and curvature parts
  return A * p1 + B * p2 + term1 + term2;
}

// Clamp i into the valid index range [0, n-1]
inline viskores::Id clamp(viskores::Id i, viskores::Id n)
{
  if (i < 0)
    return 0;
  else if (i >= n)
    return n - 1;
  else
    return i;
}

inline viskores::FloatDefault tricubicSampleRectilinear(const viskores::cont::ArrayHandle<viskores::FloatDefault>& _data,
                                                        const viskores::cont::ArrayHandle<viskores::FloatDefault>& _xCoords,
                                                        const viskores::cont::ArrayHandle<viskores::FloatDefault>& _yCoords,
                                                        const viskores::cont::ArrayHandle<viskores::FloatDefault>& _zCoords,
                                                        const viskores::Vec3f& pt)
{
  auto data = _data.ReadPortal();
  auto xCoords = _xCoords.ReadPortal();
  auto yCoords = _yCoords.ReadPortal();
  auto zCoords = _zCoords.ReadPortal();
  auto x = pt[0], y = pt[1], z = pt[2];

  viskores::Id nx = _xCoords.GetNumberOfValues();
  viskores::Id ny = _yCoords.GetNumberOfValues();
  viskores::Id nz = _zCoords.GetNumberOfValues();

  // 1) find the base indices in each dimension
  viskores::Id iu = findIndex(_xCoords, x);
  viskores::Id iv = findIndex(_yCoords, y);
  viskores::Id iw = (nz >= 4 ? findIndex(_zCoords, z) : 0);

  // If we don’t have 4 Z‐levels, do 2D bicubic in X–Y only
  if (nz < 4)
  {
    // --- bicubic: gather a 4×4 patch in X–Y at single k = clamp(iw,0,ny−1) ---
    double P2d[4 * 4];
    for (viskores::Id jj = 0; jj < 4; ++jj)
    {
      viskores::Id j = clamp(iv - 1 + jj, ny);
      for (viskores::Id ii = 0; ii < 4; ++ii)
      {
        viskores::Id i = clamp(iu - 1 + ii, nx);
        // flatten (i,j, 0) → data index
        P2d[jj * 4 + ii] = data.Get((j * nx) + i);
      }
    }
    // 3) bicubic along X → C2[4]
    double C2[4];
    for (int jj = 0; jj < 4; ++jj)
    {
      auto x0 = xCoords.Get(iu - 1), x1 = xCoords.Get(iu), x2 = xCoords.Get(iu + 1), x3 = xCoords.Get(iu + 2);
      C2[jj] = cubicInterpolateNonUniform(x0, x1, x2, x3, P2d[jj * 4 + 0], P2d[jj * 4 + 1], P2d[jj * 4 + 2], P2d[jj * 4 + 3], x);
    }
    // 4) bicubic along Y on C2 → final
    auto y0 = yCoords.Get(iv - 1), y1 = yCoords.Get(iv), y2 = yCoords.Get(iv + 1), y3 = yCoords.Get(iv + 2);
    return cubicInterpolateNonUniform(y0, y1, y2, y3, C2[0], C2[1], C2[2], C2[3], y);
  }

  // --- otherwise do full tricubic in X–Y–Z ---
  double P[4 * 4 * 4];
  for (viskores::Id kk = 0; kk < 4; ++kk)
  {
    viskores::Id k = clamp(iw - 1 + kk, nz);
    for (viskores::Id jj = 0; jj < 4; ++jj)
    {
      viskores::Id j = clamp(iv - 1 + jj, ny);
      for (viskores::Id ii = 0; ii < 4; ++ii)
      {
        viskores::Id i = clamp(iu - 1 + ii, nx);
        auto pIndex = (kk * 4 + jj) * 4 + ii;
        auto dIndex = (k * ny + j) * nx + i;
        P[pIndex] = data.Get(dIndex);
      }
    }
  }

  // interpolate in X for each (kk,jj) → Cbuf[16]
  viskores::FloatDefault Cbuf[4 * 4];
  for (int kk = 0; kk < 4; ++kk)
    for (int jj = 0; jj < 4; ++jj)
    {
      auto x0 = xCoords.Get(iu - 1), x1 = xCoords.Get(iu), x2 = xCoords.Get(iu + 1), x3 = xCoords.Get(iu + 2);
      auto p0 = P[(kk * 4 + jj) * 4 + 0], p1 = P[(kk * 4 + jj) * 4 + 1], p2 = P[(kk * 4 + jj) * 4 + 2], p3 = P[(kk * 4 + jj) * 4 + 3];
      Cbuf[kk * 4 + jj] = cubicInterpolateNonUniform(x0, x1, x2, x3, p0, p1, p2, p3, x);
    }

  // interpolate in Y → D[4]
  viskores::FloatDefault D[4];
  for (int kk = 0; kk < 4; ++kk)
  {
    auto y0 = yCoords.Get(iv - 1), y1 = yCoords.Get(iv), y2 = yCoords.Get(iv + 1), y3 = yCoords.Get(iv + 2);
    D[kk] = cubicInterpolateNonUniform(y0, y1, y2, y3, Cbuf[kk * 4 + 0], Cbuf[kk * 4 + 1], Cbuf[kk * 4 + 2], Cbuf[kk * 4 + 3], y);
  }

  // interpolate in Z
  auto z0 = zCoords.Get(iw - 1), z1 = zCoords.Get(iw), z2 = zCoords.Get(iw + 1), z3 = zCoords.Get(iw + 2);
  return cubicInterpolateNonUniform(z0, z1, z2, z3, D[0], D[1], D[2], D[3], z);
}



inline viskores::FloatDefault tricubicSampleRectilinearOLD(const viskores::cont::ArrayHandle<viskores::FloatDefault>& _data,
                                                           const viskores::cont::ArrayHandle<viskores::FloatDefault>& _xCoords,
                                                           const viskores::cont::ArrayHandle<viskores::FloatDefault>& _yCoords,
                                                           const viskores::cont::ArrayHandle<viskores::FloatDefault>& _zCoords,
                                                           const viskores::Vec3f& pt)
{
  auto data = _data.ReadPortal();
  auto xCoords = _xCoords.ReadPortal();
  auto yCoords = _yCoords.ReadPortal();
  auto zCoords = _zCoords.ReadPortal();
  auto x = pt[0], y = pt[1], z = pt[2];

  viskores::Id nx = _xCoords.GetNumberOfValues();
  viskores::Id ny = _yCoords.GetNumberOfValues();
  viskores::Id nz = _zCoords.GetNumberOfValues();

  // 1) find the base indices
  viskores::Id iu = findIndex(_xCoords, x);
  viskores::Id iv = findIndex(_yCoords, y);
  viskores::Id iw = findIndex(_zCoords, z);

  // 2) gather 4×4×4 samples into P[64]
  //    P[kk][jj][ii] → P[(kk*4 + jj)*4 + ii]
  double P[4 * 4 * 4];
  for (viskores::Id kk = 0; kk < 4; ++kk)
  {
    viskores::Id k = iw - 1 + kk;
    for (viskores::Id jj = 0; jj < 4; ++jj)
    {
      viskores::Id j = iv - 1 + jj;
      for (viskores::Id ii = 0; ii < 4; ++ii)
      {
        viskores::Id i = iu - 1 + ii;
        viskores::Id pIndex = (kk * 4 + jj) * 4 + ii;
        viskores::Id dIndex = (k * ny + j) * nx + i;
        P[pIndex] = data.Get(dIndex);
      }
    }
  }

  // 3) interpolate in X for each (kk,jj) → C[16]
  //    C[kk][jj] → C[kk*4 + jj]
  viskores::FloatDefault Cbuf[4 * 4];
  for (int kk = 0; kk < 4; ++kk)
  {
    for (int jj = 0; jj < 4; ++jj)
    {
      // knot locations and sample values
      viskores::FloatDefault x0 = xCoords.Get(iu - 1), x1 = xCoords.Get(iu), x2 = xCoords.Get(iu + 1), x3 = xCoords.Get(iu + 2);
      viskores::FloatDefault p0 = P[(kk * 4 + jj) * 4 + 0];
      viskores::FloatDefault p1 = P[(kk * 4 + jj) * 4 + 1];
      viskores::FloatDefault p2 = P[(kk * 4 + jj) * 4 + 2];
      viskores::FloatDefault p3 = P[(kk * 4 + jj) * 4 + 3];

      Cbuf[kk * 4 + jj] = cubicInterpolateNonUniform(x0, x1, x2, x3, p0, p1, p2, p3, x);
    }
  }

  // 4) interpolate in Y for each kk → D[4]
  //    C[kk][0..3] → Cbuf[kk*4 + 0..3]
  viskores::FloatDefault D[4];
  for (int kk = 0; kk < 4; ++kk)
  {
    viskores::FloatDefault y0 = yCoords.Get(iv - 1), y1 = yCoords.Get(iv), y2 = yCoords.Get(iv + 1), y3 = yCoords.Get(iv + 2);
    D[kk] = cubicInterpolateNonUniform(y0, y1, y2, y3, Cbuf[kk * 4 + 0], Cbuf[kk * 4 + 1], Cbuf[kk * 4 + 2], Cbuf[kk * 4 + 3], y);
  }

  // 5) interpolate in Z
  viskores::FloatDefault z0 = zCoords.Get(iw - 1), z1 = zCoords.Get(iw), z2 = zCoords.Get(iw + 1), z3 = zCoords.Get(iw + 2);
  return cubicInterpolateNonUniform(z0, z1, z2, z3, D[0], D[1], D[2], D[3], z);
}

inline viskores::FloatDefault evalCubicVolume(const viskores::cont::ArrayHandle<viskores::FloatDefault>& data,
                                              const viskores::Id3& dims,
                                              const viskores::Vec3f& origin,
                                              const viskores::Vec3f& spacing,
                                              const viskores::Vec3f& pt)
{
  // map world→index space
  viskores::Vec3f pt_index;
  pt_index[0] = (pt[0] - origin[0]) / spacing[0];
  pt_index[1] = (pt[1] - origin[1]) / spacing[1];
  pt_index[2] = (pt[2] - origin[2]) / spacing[2];

  // call your tricubicSample implementation
  return tricubicSample(data, dims, pt_index);
}

inline viskores::FloatDefault CubicEval(const viskores::cont::DataSet& ds, const std::string& fieldName, const viskores::Vec3f& pt)
{
  using CoordsType = viskores::cont::ArrayHandleCartesianProduct<viskores::cont::ArrayHandle<viskores::FloatDefault>,
                                                                 viskores::cont::ArrayHandle<viskores::FloatDefault>,
                                                                 viskores::cont::ArrayHandle<viskores::FloatDefault>>;
  CoordsType coords;
  coords = ds.GetCoordinateSystem().GetData().AsArrayHandle<CoordsType>();
  auto xCoords = coords.GetFirstArray();
  auto yCoords = coords.GetSecondArray();
  auto zCoords = coords.GetThirdArray();
  viskores::cont::ArrayHandle<viskores::FloatDefault> array;
  ds.GetField(fieldName).GetData().AsArrayHandle(array);

  auto val = tricubicSampleRectilinear(array, xCoords, yCoords, zCoords, pt);
  return val;

  //  coords = ds.GetCoordinateSystem().GetData().AsArrayHandle<viskores::cont::ArrayHandleCartesianProduct>();
  /*
  viskores::cont::ArrayHandleUniformPointCoordinates coords;
  coords = ds.GetCoordinateSystem().GetData().AsArrayHandle<viskores::cont::ArrayHandleUniformPointCoordinates>();
  auto portal = coords.ReadPortal();
  auto dims = coords.GetDimensions();
  auto origin = portal.GetOrigin();
  auto spacing = portal.GetSpacing();

  viskores::cont::ArrayHandle<viskores::FloatDefault> array;
  ds.GetField(fieldName).GetData().AsArrayHandle(array);
  auto val = evalCubicVolume(array, dims, origin, spacing, pt);
  return val;
  */
}