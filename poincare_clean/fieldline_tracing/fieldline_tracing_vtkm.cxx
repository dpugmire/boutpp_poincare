#include "NetCDFLoader.h"
#include "RK4_vtkm.h"
#include <adios2.h>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <numeric>

#include "parse_args.h"
#include <viskores/Particle.h>
#include <viskores/cont/Algorithm.h>
#include <viskores/cont/ArrayHandleBasic.h>
#include <viskores/cont/ArrayHandleTransform.h>
#include <viskores/cont/CellLocatorRectilinearGrid.h>
#include <viskores/cont/CellLocatorUniformGrid.h>
//#include <viskores/cont/CubicHermiteSpline.h>
#include <viskores/cont/DataSetBuilderExplicit.h>
#include <viskores/cont/DataSetBuilderRectilinear.h>
#include <viskores/cont/DataSetBuilderUniform.h>
#include <viskores/cont/Initialize.h>
#include <viskores/exec/CellInterpolate.h>
#include <viskores/filter/flow/worklet/CellInterpolationHelper.h>
#include <viskores/filter/resampling/Probe.h>
#include <viskores/io/VTKDataSetWriter.h>
#include <viskores/worklet/WorkletMapField.h>

template <typename T>
void printArray(const std::string& name, const T& arr)
{
  std::cout << name << ": [";
  auto n = arr.GetNumberOfValues();
  auto portal = arr.ReadPortal();
  for (size_t i = 0; i < n; ++i)
  {
    std::cout << portal.Get(i);
    if (i < n - 1)
      std::cout << ", ";
  }
  std::cout << "]" << std::endl;
}

template <typename T>
void printArray(const std::string& name, const std::vector<T>& arr)
{
  std::cout << name << ": [";
  auto n = arr.size();
  for (size_t i = 0; i < n; ++i)
  {
    std::cout << arr[i];
    if (i < n - 1)
      std::cout << ", ";
  }
  std::cout << "]" << std::endl;
}

class Options
{
public:
  Options(const std::string& fileName, bool transpose = false)
    : loader(fileName)
  {
    std::cout << "Options from: " << fileName << std::endl;

    this->nx = this->loader.getDim("nx");
    this->ny = this->loader.getDim("ny");
    this->nz = this->loader.getDim("nz");

    this->nx_cfr = this->loader.getDim("nx_cfr");
    this->ny_cfr = this->loader.getDim("ny_cfr");

    this->ixseps1 = this->loader.readScalar("ixsep1");
    this->ixseps2 = this->loader.readScalar("ixsep2");
    this->ixsep = this->ixseps1;
    this->jyseps1_1 = this->loader.readScalar("jyseps1_1");
    this->jyseps1_2 = this->loader.readScalar("jyseps1_2");
    this->jyseps2_1 = this->loader.readScalar("jyseps2_1");
    this->jyseps2_2 = this->loader.readScalar("jyseps2_2");
    this->nypf1 = this->jyseps1_1;
    this->nypf2 = this->jyseps2_2;
    std::cout << "read ZP" << std::endl;
    this->zperiod = this->loader.readScalar("zperiod");
    std::cout << "********** set divertor correctly!!!!" << std::endl;
    this->divertor = 1;

    this->rxy = this->loader.read2DVariable("rxy", transpose);
    this->zxy = this->loader.read2DVariable("zxy", transpose);
    this->rxy_cfr = this->loader.read2DVariable("rxy_cfr", transpose);
    this->zxy_cfr = this->loader.read2DVariable("zxy_cfr", transpose);
    this->psixy = this->loader.read2DVariable("psixy", transpose);
    this->zShift = this->loader.read2DVariable("zShift", transpose);

    this->zShift_cfr = this->loader.read2DVariable("zShift_cfr", transpose);
    this->shiftAngle = this->loader.read1DVariable("shiftAngle");
    this->dxdy = this->loader.read3DVariable("dxdy", transpose);
    this->dzdy = this->loader.read3DVariable("dzdy", transpose);

    this->dxdy_m1 = this->loader.read2DVariable("dxdy_m1", transpose);
    this->dxdy_p1 = this->loader.read2DVariable("dxdy_p1", transpose);
    this->dzdy_m1 = this->loader.read2DVariable("dzdy_m1", transpose);
    this->dzdy_p1 = this->loader.read2DVariable("dzdy_p1", transpose);
    this->nzG = this->nz; // * this->zperiod;

    std::cout << "*********************************************************************************"
                 "**********"
              << std::endl;
    std::cout << "*********************************************************************************"
                 "**********"
              << std::endl;
    std::cout << " Fix me. zarray length issue..." << std::endl;
    std::cout << "      Right now using dz = (zmax - zmin) / (nzG-1)" << std::endl;
    std::cout << "*********************************************************************************"
                 "**********"
              << std::endl;
    std::cout << "*********************************************************************************"
                 "**********"
              << std::endl;
    /*
        this->dz = (this->zmax - this->zmin) / this->nzG;
        this->init_array(this->ziarray, this->nzG+1);
        this->zarray.resize(this->nzG+1);
        for (int i = 0; i <= this->nzG; i++)
            this->zarray[i] = this->ziarray[i] * this->dz;
        */
    this->init_array(this->ziarray, this->nzG);
    this->dz = (this->zmax - this->zmin) / (this->nzG - 1);
    this->zarray.resize(this->nzG);
    for (int i = 0; i < this->nzG; i++)
      this->zarray[i] = this->ziarray[i] * this->dz;

    this->init_array(this->xiarray, this->nx);
    this->init_array(this->yiarray, this->ny);
    this->init_array(this->xiarray_cfr, this->ixsep);
    this->init_array(this->yiarray_cfr, this->nypf1 + 1, this->nypf2 + 2);

    this->init_array(this->yiarray_cfr, this->nypf1 + 1 - 1, this->nypf2 + 2 - 1);

    /*
    writeArray1DToFile(this->xiarray, "xiarray");
    writeArray1DToFile(this->xiarray_cfr, "xiarray_cfr");
    writeArray1DToFile(this->yiarray_cfr, "yiarray_cfr");
    writeArray2DToFile(this->zShift_cfr, "zshift_cfr");
    writeArray2DToFile(this->rxy_cfr, "rxy_cfr");
    writeArray2DToFile(this->zxy_cfr, "zxy_cfr");
    writeArray3DToFile(this->dxdy, "dxdy_0");
    writeArray3DToFile(this->dzdy, "dzdy_0");
    */

    // Find the index of the maximum value in the last row of rxy
    auto& lastRow = this->rxy.back();
    auto maxIt = std::max_element(lastRow.begin(), lastRow.end());
    this->jyomp = std::distance(lastRow.begin(), maxIt);

    for (const auto& row : psixy)
      this->xarray.push_back(row[this->jyomp]);

    auto minmaxIt = std::minmax_element(this->xarray.begin(), this->xarray.end());
    this->xMin = *minmaxIt.first;
    this->xMax = *minmaxIt.second;

    //VTKm stuff.
    viskores::cont::DataSetBuilderRectilinear builderRect;
    viskores::cont::DataSetBuilderUniform builder;

    std::vector<viskores::FloatDefault> yarray;
    yarray.reserve(this->ny);
    for (int i = 0; i < this->ny; i++)
      yarray.push_back(static_cast<viskores::FloatDefault>(i));

    this->Grid2D = builderRect.Create(this->xarray, yarray);

    viskores::Id3 dims(this->nx, this->ny, 1);
    viskores::Vec3f origin(0.0, 0.0, 0.0), spacing(1.0, 1.0, 1.0);
    //this->Grid2D = builder.Create(dims, origin, spacing);
    dims = { this->nx_cfr, this->ny_cfr, 1 };
    this->Grid2D_cfr = builder.Create(dims, origin, spacing);
    dims = { this->nx, this->nz, 1 };
    this->Grid2D_xz = builder.Create(dims, origin, spacing);

    //Add fields.
    this->AddField(this->psixy, "psixy", this->Grid2D);
    this->AddField(this->zShift, "zShift", this->Grid2D);
    this->AddField(this->rxy, "rxy", this->Grid2D);
    this->AddField(this->zxy, "zxy", this->Grid2D);

    this->AddField(this->zShift_cfr, "zShift_cfr", this->Grid2D_cfr);
    this->AddField(this->rxy_cfr, "rxy_cfr", this->Grid2D_cfr);
    this->AddField(this->zxy_cfr, "zxy_cfr", this->Grid2D_cfr);


    this->AddField(this->dxdy_m1, "dxdy_m1", this->Grid2D_xz);
    this->AddField(this->dxdy_p1, "dxdy_p1", this->Grid2D_xz);
    this->AddField(this->dzdy_m1, "dzdy_m1", this->Grid2D_xz);
    this->AddField(this->dzdy_p1, "dzdy_p1", this->Grid2D_xz);

    this->Grid3D = builderRect.Create(this->xarray, yarray, this->zarray);
    //viskores::Id3 dims3d(this->xarray.size(), yarray.size(), this->zarray.size());
    //viskores::Vec3f origin3d(this->xarray[0], yarray[0], this->zarray[0]);
    //viskores::Vec3f spacing3d(this->xarray[1] - origin3d[0], yarray[1] - origin3d[1], this->zarray[1] - origin3d[2]);
    //this->Grid3D = builder.Create(dims3d, origin3d, spacing3d);

    this->AddField(this->dxdy, "dxdy", this->Grid3D);
    this->AddField(this->dzdy, "dzdy", this->Grid3D);

#if 0
        this->Grid3D.PrintSummary(std::cout);
        std::cout<<"dims: "<<this->dxdy[0][0].size()<<" "<<this->dxdy[0].size()<<" "<<this->dxdy.size()<<std::endl;
        std::cout<<" npts / ncells: "<<this->Grid3D.GetNumberOfPoints()<<" "<<this->Grid3D.GetNumberOfCells()<<std::endl;

        viskores::io::VTKDataSetWriter writer("/Users/dpn/grid2D.vtk");
        writer.WriteDataSet(this->Grid2D);
        writer = viskores::io::VTKDataSetWriter("/Users/dpn/grid2D_cfr.vtk");
        writer.WriteDataSet(this->Grid2D_cfr);
        writer = viskores::io::VTKDataSetWriter("/Users/dpn/grid2D_xz.vtk");
        writer.WriteDataSet(this->Grid2D_xz);
        writer = viskores::io::VTKDataSetWriter("/Users/dpn/grid3D.vtk");
        writer.WriteDataSet(this->Grid3D);
#endif

    this->XArray = viskores::cont::make_ArrayHandle(this->xarray, viskores::CopyFlag::On);
    this->XiArray = viskores::cont::make_ArrayHandle(this->xiarray, viskores::CopyFlag::On);
    this->YArray = viskores::cont::make_ArrayHandle(yarray, viskores::CopyFlag::On);
    this->ZArray = viskores::cont::make_ArrayHandle(this->zarray, viskores::CopyFlag::On);
    this->ZiArray = viskores::cont::make_ArrayHandle(this->ziarray, viskores::CopyFlag::On);
    this->ShiftAngle = viskores::cont::make_ArrayHandle(this->shiftAngle, viskores::CopyFlag::On);
    this->YiArray_cfr = viskores::cont::make_ArrayHandle(this->yiarray_cfr, viskores::CopyFlag::On);
    this->ComputeCenter();
    this->ComputeXPoint();
    this->ComputeTheta();
    this->Theta_cfr = viskores::cont::make_ArrayHandle(this->theta_cfr, viskores::CopyFlag::On);
    std::cout << "Center: " << this->center << std::endl;
    std::cout << "XPoint: " << this->xpoint << std::endl;
  }

  void ComputeCenter()
  {
    //Compute the center of the 2D plane.
    auto jStart = this->nypf1;
    auto jEnd = (this->ny - this->nypf1) - 1;
    const viskores::Id i = 0;
    viskores::Range rRange, zRange;
    const auto rxyPortal =
      this->Grid2D.GetPointField("rxy").GetData().AsArrayHandle<viskores::cont::ArrayHandle<viskores::FloatDefault>>().ReadPortal();
    const auto zxyPortal =
      this->Grid2D.GetPointField("zxy").GetData().AsArrayHandle<viskores::cont::ArrayHandle<viskores::FloatDefault>>().ReadPortal();
    for (viskores::Id j = jStart; j <= jEnd; j++)
    {
      const auto rv = rxyPortal.Get(i + j * this->nx);
      const auto zv = zxyPortal.Get(i + j * this->nx);
      rRange.Include(rv);
      zRange.Include(zv);
    }
    this->center[0] = 0.5 * (rRange.Min + rRange.Max);
    this->center[1] = 0.5 * (zRange.Min + zRange.Max);
  }

  // Call this with isIxsepOneBased=true if ixsep came from MATLAB; else false.
  void ComputeXPoint(bool isIxsepOneBased = true)
  {
    const viskores::Id nx = this->nx, ny = this->ny, nypf1 = this->nypf1;
    if (nx < 2 || ny < 2)
      throw std::runtime_error("Grid too small for X-point");

    // Choose the two columns that straddle the separatrix, honoring the indexing base.
    const viskores::Id i0 = isIxsepOneBased ? (this->ixsep - 1) : this->ixsep;
    const viskores::Id i1 = i0 + 1;
    if (i0 < 0 || i1 >= nx)
      throw std::runtime_error("ixsep out of range for separatrix pair");

    // The four j positions used by the MATLAB (converted to 0-based):
    const viskores::Id j0 = nypf1 - 1;      // nypf1
    const viskores::Id j1 = nypf1;          // nypf1+1
    const viskores::Id j2 = ny - nypf1 - 1; // ny - nypf1
    const viskores::Id j3 = ny - nypf1;     // ny - nypf1 + 1
    if (j0 < 0 || j3 >= ny)
      throw std::runtime_error("nypf1/ny inconsistent with grid size");

    const auto rxyPortal =
      this->Grid2D.GetPointField("rxy").GetData().AsArrayHandle<viskores::cont::ArrayHandle<viskores::FloatDefault>>().ReadPortal();
    const auto zxyPortal =
      this->Grid2D.GetPointField("zxy").GetData().AsArrayHandle<viskores::cont::ArrayHandle<viskores::FloatDefault>>().ReadPortal();

    auto rAt = [&](viskores::Id i, viskores::Id j) { return rxyPortal.Get(i + j * nx); };
    auto zAt = [&](viskores::Id i, viskores::Id j) { return zxyPortal.Get(i + j * nx); };

    auto mid = [&](viskores::Id j)
    {
      const auto rx = static_cast<viskores::FloatDefault>(0.5) * (rAt(i0, j) + rAt(i1, j));
      const auto rz = static_cast<viskores::FloatDefault>(0.5) * (zAt(i0, j) + zAt(i1, j));
      return std::pair<viskores::FloatDefault, viskores::FloatDefault>{ rx, rz };
    };

    const auto [x0, y0] = mid(j0);
    const auto [x1, y1] = mid(j1);
    const auto [x2, y2] = mid(j2);
    const auto [x3, y3] = mid(j3);

    this->xpoint[0] = static_cast<viskores::FloatDefault>(0.25) * (x0 + x1 + x2 + x3);
    this->xpoint[1] = static_cast<viskores::FloatDefault>(0.25) * (y0 + y1 + y2 + y3);
  }

  void ComputeTheta()
  {
    const viskores::Id nx = this->nx;
    const viskores::Id ny = this->ny;
    const viskores::Id nypf1 = this->nypf1; // MATLAB's nypf1 (note: 0-based index nypf1 maps to MATLAB nypf1+1)
    if (nx < 1 || ny < 1)
      throw std::runtime_error("Grid too small in ComputeTheta");

    // --- Base vectors (Z components are 0, we stay in 2D) ---
    // u = center - xpoint (reference direction for theta=0 after shift)
    const viskores::FloatDefault ux = static_cast<viskores::FloatDefault>(this->center[0] - this->xpoint[0]);
    const viskores::FloatDefault uy = static_cast<viskores::FloatDefault>(this->center[1] - this->xpoint[1]);

    // Access fields
    const auto rxyPortal =
      this->Grid2D.GetPointField("rxy").GetData().AsArrayHandle<viskores::cont::ArrayHandle<viskores::FloatDefault>>().ReadPortal();
    const auto zxyPortal =
      this->Grid2D.GetPointField("zxy").GetData().AsArrayHandle<viskores::cont::ArrayHandle<viskores::FloatDefault>>().ReadPortal();

    auto rAt = [&](viskores::Id i, viskores::Id j) -> viskores::FloatDefault
    { return static_cast<viskores::FloatDefault>(rxyPortal.Get(i + j * nx)); };
    auto zAt = [&](viskores::Id i, viskores::Id j) -> viskores::FloatDefault
    { return static_cast<viskores::FloatDefault>(zxyPortal.Get(i + j * nx)); };

    // --- theta(j) = atan2( |u x v|, u·v ) / pi, where v = center - (rxy(0,j), zxy(0,j)) ---
    this->theta.resize(static_cast<std::size_t>(ny), 0.0);
    const viskores::FloatDefault pi = 3.141592653589793238462643383279502884;

    for (viskores::Id j = 0; j < ny; ++j)
    {
      const viskores::FloatDefault vx = static_cast<viskores::FloatDefault>(this->center[0]) - rAt(0, j); // rxy(1, j+1) in MATLAB
      const viskores::FloatDefault vy = static_cast<viskores::FloatDefault>(this->center[1]) - zAt(0, j); // zxy(1, j+1)

      const viskores::FloatDefault dot = ux * vx + uy * vy;
      // 2D cross product magnitude equals |ux*vy - uy*vx| (z-component of 3D cross)
      const viskores::FloatDefault crs = std::abs(ux * vy - uy * vx);

      const viskores::FloatDefault ang = std::atan2(crs, dot); // in [0, pi]
      this->theta[static_cast<std::size_t>(j)] = ang / pi;     // scale to [0,1], like MATLAB
    }

    // --- Unwrap / reflect like the MATLAB code ---
    // [c, itheta] = max(theta);
    // theta(itheta:ny) = 2 - theta(itheta:ny);
    {
      auto it = std::max_element(this->theta.begin(), this->theta.end());
      const std::size_t itheta = static_cast<std::size_t>(std::distance(this->theta.begin(), it));
      for (std::size_t j = itheta; j < theta.size(); ++j)
        this->theta[j] = 2.0 - this->theta[j];
    }

    // recompute max; if (itheta != ny) then theta(itheta:ny) = 4 - theta(itheta:ny);
    {
      auto it = std::max_element(theta.begin(), theta.end());
      const std::size_t itheta = static_cast<std::size_t>(std::distance(theta.begin(), it));
      if (itheta != static_cast<std::size_t>(ny - 1))
      {
        for (std::size_t j = itheta; j < theta.size(); ++j)
          this->theta[j] = 4.0 - this->theta[j];
      }
    }

    // --- Shift so that theta at (1, nypf1+1) equals zero (MATLAB) ---
    // MATLAB uses 1-based index (nypf1+1). Our 0-based array uses index nypf1.
    if (nypf1 < 0 || nypf1 >= ny)
      throw std::runtime_error("nypf1 out of range in ComputeTheta");
    const viskores::FloatDefault ref = theta[static_cast<std::size_t>(nypf1)];
    for (viskores::FloatDefault& t : this->theta)
      t -= ref;

    printArray("theta", this->theta);

    //compute theta_cfr.
    // ----- Closed Flux Region (CFR) indices and theta_cfr -----
    if (this->ixsep < 0 || this->ixsep >= nx)
      throw std::runtime_error("ixsep out of range for CFR");
    if (this->nypf1 < 0 || this->nypf2 < this->nypf1 || this->nypf2 >= ny)
      throw std::runtime_error("nypf1/nypf2 out of range for CFR");

    // xiarray_cfr = 1:ixsep  (MATLAB, 1-based)  -->  0..ixsep (0-based, inclusive)
    std::vector<viskores::FloatDefault> xiarrayCfr;
    xiarrayCfr.reserve(static_cast<std::size_t>(this->ixsep + 1));
    for (viskores::Id i = 0; i <= this->ixsep; ++i)
      xiarrayCfr.push_back(i);

    // yiarray_cfr = nypf1+1:nypf2+1  (MATLAB)  -->  nypf1..nypf2 (0-based, inclusive)
    std::vector<viskores::FloatDefault> yiarrayCfr;
    yiarrayCfr.reserve(static_cast<std::size_t>(this->nypf2 - this->nypf1 + 1));
    for (viskores::Id j = this->nypf1; j <= this->nypf2; ++j)
      yiarrayCfr.push_back(j);

    // theta_cfr = theta(yiarray_cfr); theta_cfr(end) = 2.0;
    std::vector<viskores::FloatDefault> thetaCfr;
    thetaCfr.reserve(yiarrayCfr.size());
    for (auto j : yiarrayCfr)
      thetaCfr.push_back(this->theta[static_cast<std::size_t>(j)]);
    if (!thetaCfr.empty())
      thetaCfr.back() = static_cast<viskores::FloatDefault>(2.0);

    // (Optional) stash as members for later use:
    this->xiarray_cfr = std::move(xiarrayCfr); // std::vector<viskores::Id>
    this->yiarray_cfr = std::move(yiarrayCfr); // std::vector<viskores::Id>
    this->theta_cfr = std::move(thetaCfr);     // std::vector<viskores::FloatDefault>
    printArray("theta_cfr: ", this->theta_cfr);
  }

  void AddField(const std::vector<std::vector<double>>& vals, const std::string& fieldName, viskores::cont::DataSet& ds)
  {
    viskores::Id n_y = vals[0].size();
    viskores::Id n_x = vals.size();

    std::vector<double> field;
    field.reserve(n_x * n_y);
    for (viskores::Id j = 0; j < n_y; ++j)
      for (viskores::Id i = 0; i < n_x; ++i)
        field.push_back(vals[i][j]);
    ds.AddPointField(fieldName, field);
  }

  void AddField(const std::vector<std::vector<std::vector<double>>>& vals, const std::string& fieldName, viskores::cont::DataSet& ds)
  {
    viskores::Id n_z = vals[0][0].size();
    viskores::Id n_y = vals[0].size();
    viskores::Id n_x = vals.size();
    std::cout << "Add Field3D: " << fieldName << " " << n_x << " " << n_y << " " << n_z << std::endl;


    std::vector<double> field;
    field.reserve(n_x * n_y * n_z);
    for (viskores::Id k = 0; k < n_z; ++k)
      for (viskores::Id j = 0; j < n_y; ++j)
        for (viskores::Id i = 0; i < n_x; ++i)
          field.push_back(vals[i][j][k]);

    ds.AddPointField(fieldName, field);
  }

  void init_array(std::vector<double>& arr, size_t n) { this->init_array(arr, 0, n); }

  void init_array(std::vector<double>& arr, size_t n0, size_t n1)
  {
    size_t n = n1 - n0;
    arr.resize(n);
    for (size_t i = 0; i < n; i++)
      arr[i] = static_cast<double>(i + n0);
  }

  //members.
  int nx = 260;
  int ny = 128;
  int nz = 256;
  int nx_cfr = 195;
  int ny_cfr = 97;
  double xMin, xMax;
  int jyomp;
  double zmin = 0.0;
  double zmax = 2.0 * M_PI;
  double dz;
  int zperiod = -1;
  int jyseps1_1 = -1;
  int jyseps1_2 = -1;
  int jyseps2_1 = -1;
  int jyseps2_2 = -1;
  int ixseps1 = -1;
  int ixseps2 = -1;
  int ixsep = -1;
  int nypf1 = -1;
  int nypf2 = -1;

  int direction = 1;
  int divertor = 0;

  int nzG;
  std::vector<viskores::FloatDefault> xiarray, xarray, xiarray_cfr;
  std::vector<viskores::FloatDefault> yiarray, yiarray_cfr;
  std::vector<viskores::FloatDefault> ziarray, zarray;

  std::vector<double> dy;
  double dy0 = 0.0;
  std::vector<std::vector<std::vector<double>>> dxdy, dzdy;
  std::vector<std::vector<double>> rxy, zxy, psixy, zShift;
  std::vector<std::vector<double>> rxy_cfr, zxy_cfr, zShift_cfr;
  std::vector<std::vector<double>> dxdy_p1, dzdy_p1, dxdy_m1, dzdy_m1;
  std::vector<double> shiftAngle;

  NetCDFLoader loader;

  //vtkm data.
  viskores::cont::DataSet Grid2D, Grid2D_cfr, Grid2D_xz;
  viskores::cont::DataSet Grid3D;
  viskores::cont::ArrayHandle<viskores::FloatDefault> XArray, YArray, ZArray;
  viskores::cont::ArrayHandle<viskores::FloatDefault> XiArray, YiArray, ZiArray;
  viskores::cont::ArrayHandle<viskores::FloatDefault> ShiftAngle;
  viskores::cont::ArrayHandle<viskores::FloatDefault> YiArray_cfr, Theta_cfr;
  viskores::Vec2f center, xpoint;
  std::vector<viskores::FloatDefault> theta, theta_cfr;
};

// Interpolate psixy[i][j] along i (x-index). i in [0..nx-1], j in [0..ny-1].
template <typename T>
inline T IndexInterp(const std::vector<std::vector<T>>& psixy, double xind, int jy)
{
  static_assert(std::is_arithmetic<T>::value, "T must be arithmetic");

  const size_t nx = psixy.size();
  if (nx == 0)
    throw std::invalid_argument("psixy has zero rows");
  if (jy < 0 || psixy[0].empty() || static_cast<size_t>(jy) >= psixy[0].size())
    throw std::out_of_range("jy out of range");

  // Clamp to edges
  if (xind <= 0.0)
    return psixy[0][static_cast<size_t>(jy)];
  if (xind >= double(nx - 1))
    return psixy[nx - 1][static_cast<size_t>(jy)];

  const size_t i0 = static_cast<size_t>(std::floor(xind));
  const size_t i1 = i0 + 1;
  const double t = xind - double(i0); // in [0,1)

  using U = typename std::common_type<T, double>::type; // widen for math
  const U v0 = static_cast<U>(psixy[i0][static_cast<size_t>(jy)]);
  const U v1 = static_cast<U>(psixy[i1][static_cast<size_t>(jy)]);
  const U v = v0 + static_cast<U>(t) * (v1 - v0);
  return static_cast<T>(v);
}

void SaveOutput(const cli::Options<viskores::FloatDefault>& cliOpts,
                const viskores::cont::ArrayHandle<viskores::Id>& validIds,
                const viskores::cont::ArrayHandle<viskores::Id>& validResultIndices,
                const viskores::cont::ArrayHandle<viskores::Vec4f>& validResultRZThetaPsi)
{
  adios2::ADIOS adios(MPI_COMM_WORLD);
  auto io = adios.DeclareIO("punctures");
  int localNum = static_cast<int>(validResultIndices.GetNumberOfValues());
  int totalNum = 0;
  MPI_Allreduce(&localNum, &totalNum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  int rankOffset = 0;
  MPI_Exscan(&localNum, &rankOffset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (cliOpts.rank == 0)
    rankOffset = 0;

  std::vector<viskores::FloatDefault> dataR(localNum), dataZ(localNum), dataPsi(localNum), dataTheta(localNum);
  auto ptPortal = validResultRZThetaPsi.ReadPortal();
  for (viskores::Id i = 0; i < localNum; i++)
  {
    dataR[i] = ptPortal.Get(i)[0];
    dataZ[i] = ptPortal.Get(i)[1];
    dataPsi[i] = ptPortal.Get(i)[2];
    dataTheta[i] = ptPortal.Get(i)[3];
  }

  auto varID =
    io.DefineVariable<viskores::Id>("ID", { static_cast<size_t>(totalNum) }, { static_cast<size_t>(rankOffset) }, { static_cast<size_t>(localNum) });
  auto varIdx =
    io.DefineVariable<viskores::Id>("Idx", { static_cast<size_t>(totalNum) }, { static_cast<size_t>(rankOffset) }, { static_cast<size_t>(localNum) });
  auto rVar = io.DefineVariable<viskores::FloatDefault>(
    "R", { static_cast<size_t>(totalNum) }, { static_cast<size_t>(rankOffset) }, { static_cast<size_t>(localNum) });
  auto zVar = io.DefineVariable<viskores::FloatDefault>(
    "Z", { static_cast<size_t>(totalNum) }, { static_cast<size_t>(rankOffset) }, { static_cast<size_t>(localNum) });
  auto psiVar = io.DefineVariable<viskores::FloatDefault>(
    "Psi", { static_cast<size_t>(totalNum) }, { static_cast<size_t>(rankOffset) }, { static_cast<size_t>(localNum) });
  auto thetaVar = io.DefineVariable<viskores::FloatDefault>(
    "Theta", { static_cast<size_t>(totalNum) }, { static_cast<size_t>(rankOffset) }, { static_cast<size_t>(localNum) });

  viskores::cont::ArrayHandleBasic<viskores::Id> idBasic(validIds);
  viskores::cont::ArrayHandleBasic<viskores::Id> idxBasic(validResultIndices);

  auto engine = io.Open(cliOpts.GetOutputFileName(), adios2::Mode::Write);
  engine.BeginStep();
  engine.Put(varID, idBasic.GetReadPointer());
  engine.Put(varIdx, idxBasic.GetReadPointer());
  engine.Put(rVar, dataR.data());
  engine.Put(zVar, dataZ.data());
  engine.Put(psiVar, dataPsi.data());
  engine.Put(thetaVar, dataTheta.data());
  engine.EndStep();
  engine.Close();
}

int main(int argc, char* argv[])
{
  cli::Options<viskores::FloatDefault> cliOpts;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &cliOpts.numRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &cliOpts.rank);
  std::cout << "Rank: " << cliOpts.rank << " of " << cliOpts.numRanks << std::endl;

  cli::parseArgs<viskores::FloatDefault>(argc, argv, cliOpts);
  std::cout << "Options: " << std::endl;
  printArray("xind: ", cliOpts.xind);
  printArray("IDs: ", cliOpts.IDs);
  std::cout << std::endl << std::endl;
  std::cout << "apar:: " << cliOpts.apar << std::endl;
  std::cout << "maxPunc:: " << cliOpts.maxpunc << std::endl;

  if (cliOpts.numRanks == 1)
  {
    auto opts = viskores::cont::InitializeOptions::DefaultAnyDevice;
    auto config = viskores::cont::Initialize(argc, argv, opts);
    //viskores::cont::GetRuntimeDeviceTracker().ForceDevice(viskores::cont::DeviceAdapterTagKokkos{});
    viskores::cont::GetRuntimeDeviceTracker().ForceDevice(viskores::cont::DeviceAdapterTagTBB{});
    std::cout << "Using TBB" << std::endl;
  }
  else
    viskores::cont::GetRuntimeDeviceTracker().ForceDevice(viskores::cont::DeviceAdapterTagSerial{});

  std::string fname = cliOpts.apar;

  bool transpose = true;
  if (fname.find("python") != std::string::npos)
    transpose = true;

  Options opts(fname, transpose);

  int divertor = 1; //single null
  double xind = 0.0f;

  std::vector<viskores::Vec3f> points;
  int idx = 0;
  for (const auto& iline : cliOpts.xind)
  {
    xind = static_cast<viskores::FloatDefault>(iline);
    int yyy = opts.jyomp;
    int yStart = opts.jyomp;
    int zzz = 0;
    auto zStart = opts.zarray[zzz];

    auto xStart = IndexInterp(opts.psixy, xind, opts.jyomp);
    viskores::Vec3f p0(xStart, static_cast<viskores::FloatDefault>(yStart), zStart);
    points.push_back(p0);
  }
  std::cout << "numPts= " << points.size() << std::endl;
  std::cout << "numPuncs= " << cliOpts.maxpunc << std::endl;

  const auto& grid3D = opts.Grid3D;
  const auto& grid2D = opts.Grid2D;
  viskores::cont::CellLocatorRectilinearGrid locator2D, locator3D;
  locator3D.SetCoordinates(grid3D.GetCoordinateSystem());
  locator3D.SetCellSet(grid3D.GetCellSet());
  locator3D.Update();
  locator2D.SetCoordinates(grid2D.GetCoordinateSystem());
  locator2D.SetCellSet(grid2D.GetCellSet());
  locator2D.Update();
  viskores::Id maxSteps = cliOpts.maxpunc * 100;

  RK4Worklet worklet(cliOpts.maxpunc, maxSteps);
  worklet.ds3D = grid3D;
  worklet.ds2D = grid2D;
  worklet.grid2DBounds = grid2D.GetCoordinateSystem().GetBounds();
  worklet.nypf1 = opts.nypf1;
  worklet.nypf2 = opts.nypf2;
  worklet.ixsep1 = opts.ixseps1;
  worklet.ixsep2 = opts.ixseps2;
  worklet.divertor = opts.divertor;
  worklet.ny = opts.ny;
  worklet.direction = 1;

  BoutppField boutppField(
    opts.Grid3D, opts.Grid2D, opts.XiArray, opts.XArray, opts.YArray, opts.ZiArray, opts.ZArray, opts.ShiftAngle, opts.YiArray_cfr, opts.Theta_cfr);

  viskores::cont::Invoker invoker;
  auto inPts = viskores::cont::make_ArrayHandle<viskores::Vec3f>(points, viskores::CopyFlag::On);
  auto xinds = viskores::cont::make_ArrayHandle<viskores::FloatDefault>(cliOpts.xind, viskores::CopyFlag::On);
  auto IDs = viskores::cont::make_ArrayHandle<viskores::Id>(cliOpts.IDs, viskores::CopyFlag::On);
  viskores::cont::ArrayHandle<viskores::FloatDefault> dxdyField, dzdyField, rxyField, zShiftField;
  viskores::cont::ArrayHandle<viskores::Vec3f> puncturesCart, resultIndex;
  grid3D.GetField("dxdy").GetData().AsArrayHandle<viskores::FloatDefault>(dxdyField);
  grid3D.GetField("dzdy").GetData().AsArrayHandle<viskores::FloatDefault>(dzdyField);
  grid2D.GetField("rxy").GetData().AsArrayHandle<viskores::FloatDefault>(rxyField);
  grid2D.GetField("zShift").GetData().AsArrayHandle<viskores::FloatDefault>(zShiftField);


  auto start = std::chrono::high_resolution_clock::now();

  viskores::cont::ArrayHandle<viskores::Vec4f> punctures;
  punctures.Allocate(inPts.GetNumberOfValues() * cliOpts.maxpunc);

  viskores::cont::ArrayHandle<viskores::Id> punctureIDs, punctureIndex;
  punctureIDs.AllocateAndFill(inPts.GetNumberOfValues() * cliOpts.maxpunc, -1);
  punctureIndex.AllocateAndFill(inPts.GetNumberOfValues() * cliOpts.maxpunc, -1);

  invoker(worklet, IDs, inPts, boutppField, grid3D.GetCellSet(), grid2D.GetCellSet(), punctureIDs, punctureIndex, punctures);

  viskores::cont::ArrayHandle<viskores::Id> validIndex, validIDs;
  viskores::cont::ArrayHandle<viskores::Vec4f> validPunctures;

  struct GreaterThanEqZero
  {
    VISKORES_EXEC_CONT bool operator()(viskores::Id x) const { return x >= 0; }
  };
  viskores::cont::Algorithm::CopyIf(punctures, punctureIndex, validPunctures, GreaterThanEqZero());
  viskores::cont::Algorithm::CopyIf(punctureIDs, punctureIndex, validIDs, GreaterThanEqZero());
  viskores::cont::Algorithm::CopyIf(punctureIndex, punctureIndex, validIndex, GreaterThanEqZero());

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;
  std::cout << " NUM particles= " << cliOpts.xind.size() << " numPunc= " << cliOpts.maxpunc << std::endl;

  for (viskores::Id i = 0; i < validPunctures.GetNumberOfValues(); i++)
  {
    auto pt = validPunctures.ReadPortal().Get(i);
    auto id = validIDs.ReadPortal().Get(i);
    auto idx = validIndex.ReadPortal().Get(i);
    cliOpts.puncSplineOut << id << ", " << idx << ", " << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << std::endl;
  }
  SaveOutput(cliOpts, validIDs, validIndex, validPunctures);

  MPI_Finalize();
  return 0;
}
