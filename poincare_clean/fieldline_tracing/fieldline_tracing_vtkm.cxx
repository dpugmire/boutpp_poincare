#include "NetCDFLoader.h"
#include "RK4_vtkm.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>

#include "parse_args.h"
#include <viskores/Particle.h>
#include <viskores/cont/Algorithm.h>
#include <viskores/cont/ArrayHandleTransform.h>
#include <viskores/cont/CellLocatorRectilinearGrid.h>
#include <viskores/cont/CellLocatorUniformGrid.h>
#include <viskores/cont/CubicHermiteSpline.h>
#include <viskores/cont/DataSetBuilderExplicit.h>
#include <viskores/cont/DataSetBuilderRectilinear.h>
#include <viskores/cont/DataSetBuilderUniform.h>
#include <viskores/cont/Initialize.h>
#include <viskores/exec/CellInterpolate.h>
#include <viskores/filter/flow/worklet/CellInterpolationHelper.h>
#include <viskores/filter/resampling/Probe.h>
#include <viskores/io/VTKDataSetWriter.h>
#include <viskores/worklet/WorkletMapField.h>

std::ofstream trajsplineFid("./trajspline.v.txt", std::ofstream::out);
std::ofstream puncFid("./punc.v.txt");
std::ofstream punc_ip_Fid("./punc_ip.v.txt");
std::ofstream puncFid2("./punc2.v.txt");
std::ofstream puncSplineFid("./puncspline.v.txt");
std::ofstream rawPunc("./rawpunc.v.txt", std::ofstream::out);
std::ofstream trajOut("./traj.v.txt", std::ofstream::out);
std::ofstream tanOut("./tan.v.txt", std::ofstream::out);
std::ofstream stepOut("./steps.v.txt", std::ofstream::out);
std::ofstream rk4Out("./rk4.v.txt", std::ofstream::out);

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

    std::cout << "Create rectilinear: " << this->xarray.size() << " " << yarray.size() << " " << this->zarray.size() << std::endl;
    this->Grid3D = builderRect.Create(this->xarray, yarray, this->zarray);
    //viskores::Id3 dims3d(this->xarray.size(), yarray.size(), this->zarray.size());
    //viskores::Vec3f origin3d(this->xarray[0], yarray[0], this->zarray[0]);
    //viskores::Vec3f spacing3d(this->xarray[1] - origin3d[0], yarray[1] - origin3d[1], this->zarray[1] - origin3d[2]);
    //this->Grid3D = builder.Create(dims3d, origin3d, spacing3d);

    std::cout << "Create rectilinear: " << this->xarray.size() << " " << yarray.size() << " " << this->zarray.size() << std::endl;
    this->AddField(this->dxdy, "dxdy", this->Grid3D);
    this->AddField(this->dzdy, "dzdy", this->Grid3D);
#ifndef VISKORES_HIP
    std::cout << std::setprecision(15) << "Grid3d bounds: " << this->Grid3D.GetCoordinateSystem().GetBounds() << std::endl;
    std::cout << "  dz= " << this->dz << " z0: " << this->zarray[0] << " " << this->zarray[1] << " ... " << this->zarray[this->zarray.size() - 2]
              << " " << this->zarray[this->zarray.size() - 1] << std::endl;
#endif


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

    this->ComputeCenter();
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
    std::cout << "******************* Center of 2D plane: " << this->center << std::endl;
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

    std::cout << "  **** field size: " << field.size() << std::endl;
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
  std::vector<double> xiarray, xarray, xiarray_cfr;
  std::vector<double> yiarray, yiarray_cfr;
  std::vector<double> ziarray, zarray;

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
  viskores::Vec2f center;
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

int main(int argc, char* argv[])
{
  cli::Options<viskores::FloatDefault> cliOpts;
  cli::parseArgs<viskores::FloatDefault>(argc, argv, cliOpts);
  std::cout << "Options: " << std::endl;
  printArray("xind: ", cliOpts.xind);
  std::cout << std::endl << std::endl;
  std::cout << "apar:: " << cliOpts.apar << std::endl;
  std::cout << "maxPunc:: " << cliOpts.maxpunc << std::endl;

  if (true)
  {
    auto opts = viskores::cont::InitializeOptions::DefaultAnyDevice;
    auto config = viskores::cont::Initialize(argc, argv, opts);
    //viskores::cont::GetRuntimeDeviceTracker().ForceDevice(viskores::cont::DeviceAdapterTagKokkos{});
    viskores::cont::GetRuntimeDeviceTracker().ForceDevice(viskores::cont::DeviceAdapterTagTBB{});
    std::cout << "Using TBB" << std::endl;
  }
  else
    viskores::cont::GetRuntimeDeviceTracker().ForceDevice(viskores::cont::DeviceAdapterTagSerial{});

  trajsplineFid << "ID, STEP, X, Y, Z, REGION\n";
  trajOut << "ID, STEP, X, Y, Z" << std::endl;
  tanOut << "ID, STEP, X, Y, Z, VX, VY, VZ" << std::endl;
  rawPunc << "ID, STEP, X, Y, Z\n";
  puncFid << "ID, STEP, X, Y, Z\n";
  puncFid2 << "ID, STEP, X, Y, Z\n";
  puncSplineFid << "Xind0, STEP, X, Y, Z, THETA, PSI\n";
  punc_ip_Fid << "ID, STEP, X, Y, Z\n";
  stepOut << "ID, STEP, X, Y, Z\n";
  rk4Out << "ID, STEP, X, Y, Z\n";
  auto TRAJ_FID = fopen("./pt_traj.c.txt", "w");
  fprintf(TRAJ_FID, "IT, xind, yEnd, zind, REG, zEnd\n");

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
  worklet.grid3DBounds = grid3D.GetCoordinateSystem().GetBounds();
  worklet.grid2DBounds = grid2D.GetCoordinateSystem().GetBounds();
  worklet.nypf1 = opts.nypf1;
  worklet.nypf2 = opts.nypf2;
  worklet.ixsep1 = opts.ixseps1;
  worklet.ixsep2 = opts.ixseps2;
  worklet.divertor = opts.divertor;
  worklet.ny = opts.ny;
  worklet.direction = 1;
  worklet.Center = opts.center;
  //printArray("ZiArray", opts.ZiArray);
  //printArray("ZArray", opts.ZArray);
  BoutppField boutppField(opts.Grid3D, opts.Grid2D, opts.XiArray, opts.XArray, opts.YArray, opts.ZiArray, opts.ZArray, opts.ShiftAngle);

  viskores::cont::Invoker invoker;
  auto inPts = viskores::cont::make_ArrayHandle<viskores::Vec3f>(points, viskores::CopyFlag::On);
  auto xinds = viskores::cont::make_ArrayHandle<viskores::FloatDefault>(cliOpts.xind, viskores::CopyFlag::On);
  viskores::cont::printSummary_ArrayHandle(xinds, std::cout, true);
  viskores::cont::ArrayHandle<viskores::FloatDefault> dxdyField, dzdyField, rxyField, zShiftField;
  viskores::cont::ArrayHandle<viskores::Vec3f> puncturesCart, resultIndex;
  grid3D.GetField("dxdy").GetData().AsArrayHandle<viskores::FloatDefault>(dxdyField);
  grid3D.GetField("dzdy").GetData().AsArrayHandle<viskores::FloatDefault>(dzdyField);
  grid2D.GetField("rxy").GetData().AsArrayHandle<viskores::FloatDefault>(rxyField);
  grid2D.GetField("zShift").GetData().AsArrayHandle<viskores::FloatDefault>(zShiftField);
  std::cout << "NumPoints= " << grid3D.GetNumberOfPoints() << std::endl;
  grid3D.PrintSummary(std::cout);


  auto start = std::chrono::high_resolution_clock::now();

  viskores::cont::ArrayHandle<viskores::Vec3f> result;
  viskores::cont::ArrayHandle<viskores::Vec2f> resultPsiTheta;
  viskores::cont::ArrayHandle<viskores::FloatDefault> resultStep;
  result.Allocate(inPts.GetNumberOfValues() * cliOpts.maxpunc);
  resultPsiTheta.Allocate(inPts.GetNumberOfValues() * cliOpts.maxpunc);

  resultStep.Allocate(inPts.GetNumberOfValues() * cliOpts.maxpunc);

  viskores::cont::ArrayHandle<viskores::FloatDefault> resultIndices;
  resultIndices.Allocate(inPts.GetNumberOfValues() * cliOpts.maxpunc);

  viskores::cont::ArrayHandle<bool> validPuncs;
  validPuncs.AllocateAndFill(inPts.GetNumberOfValues() * cliOpts.maxpunc, false);

  invoker(
    worklet, inPts, xinds, boutppField, grid3D.GetCellSet(), grid2D.GetCellSet(), validPuncs, result, resultPsiTheta, resultStep, resultIndices);

  viskores::cont::ArrayHandle<viskores::Vec3f> validResult;
  viskores::cont::ArrayHandle<viskores::Vec2f> validResultPsiTheta;
  viskores::cont::ArrayHandle<viskores::FloatDefault> validResultIndices;
  viskores::cont::ArrayHandle<viskores::FloatDefault> validStep;

  viskores::cont::Algorithm::CopyIf(result, validPuncs, validResult);
  viskores::cont::Algorithm::CopyIf(resultPsiTheta, validPuncs, validResultPsiTheta);

  viskores::cont::Algorithm::CopyIf(resultStep, validPuncs, validStep);

  viskores::cont::Algorithm::CopyIf(resultIndices, validPuncs, validResultIndices);
  std::cout << "******** " << result.GetNumberOfValues() << " ----> " << validResult.GetNumberOfValues() << std::endl;


  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;
  std::cout << " NUM particles= " << cliOpts.xind.size() << " numPunc= " << cliOpts.maxpunc << std::endl;

  //output all punctures.
  std::cout << "numValid puncs: " << validPuncs.GetNumberOfValues() << " # valid results " << validResult.GetNumberOfValues() << std::endl;
  //viskores::cont::printSummary_ArrayHandle(validPuncs, std::cout, true);

  for (viskores::Id i = 0; i < validResult.GetNumberOfValues(); i++)
  {
    auto pt = validResult.ReadPortal().Get(i);
    auto psiTheta = validResultPsiTheta.ReadPortal().Get(i);
    //psiTheta[1] = ComputeTheta(opts.center, viskores::Vec2f(pt[1], pt[2]));
    auto step = validStep.ReadPortal().Get(i);
    auto x0 = validResultIndices.ReadPortal().Get(i);
    puncSplineFid << x0 << ", " << step << ", " << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << psiTheta[1] << ", " << psiTheta[0] << std::endl;
  }

  return 0;
}
