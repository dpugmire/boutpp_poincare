#include "Evaluation.h"
#include "NetCDFLoader.h"
#include "RK4.h"
#include "RK4_vtkm.h"
#include "SplineInterpolation.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>

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
std::ofstream trajspline("./trajspline.v.txt", std::ofstream::out);
std::ofstream puncFid("./punc.v.txt");
std::ofstream punc_ip_Fid("./punc_ip.v.txt");
std::ofstream puncFid2("./punc2.v.txt");
std::ofstream puncSplineFid("./puncspline.v.txt");
std::ofstream rawPunc("./rawpunc.v.txt", std::ofstream::out);
std::ofstream trajOut("./traj.v.txt", std::ofstream::out);
std::ofstream tanOut("./tan.v.txt", std::ofstream::out);
std::ofstream stepOut("./steps.v.txt", std::ofstream::out);
std::ofstream rk4Out("./rk4.v.txt", std::ofstream::out);

double double_mod(double val, double mod_base)
{
  double result = std::fmod(val, mod_base);
  if (result < 0)
    result += mod_base;
  return result;
}

void min_max_values(const std::vector<double>& arr, double& vmin, double& vmax)
{
  auto [minIt, maxIt] = std::minmax_element(arr.begin(), arr.end());
  vmin = *minIt;
  vmax = *maxIt;
}

class Point
{
public:
  Point() = default;
#if 0
    template <typename T, typename U, typename V>
    Point(T _x, U _y, V _z, int _xi, int _yi, int _iter)
    : x(static_cast<double>(_x)), y(static_cast<double>(_y)), z(static_cast<double>(_z)), xi(_xi), yi(_yi), iter(_iter) {}
    Point(const Point& pt) : x(pt.x), y(pt.y), z(pt.z), xi(pt.xi), yi(pt.yi), iter(pt.iter) {}

    Point& operator=(const Point& pt)
    {
        if (this != &pt)
        {
            x = pt.x;
            y = pt.y;
            z = pt.z;
            xi = pt.xi;
            yi = pt.yi;
            iter = pt.iter;
        }
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& pt)
    {
        os << "(" << pt.x << ", " << pt.y << ", " << pt.z << ")";
        return os;
    }
#endif

  double traj1 = 0.0, traj2 = 0.0, traj3 = 0.0, traj4 = 0.0, traj5 = 0.0, traj6 = 0.0, traj7 = 0.0;
  /*
    int xi = 0;
    int yi = 0;
    int iter = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    */
};

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
    this->nzG = this->nz * this->zperiod;

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

    writeArray1DToFile(this->xiarray_cfr, "xiarray_cfr");
    writeArray1DToFile(this->yiarray_cfr, "yiarray_cfr");
    writeArray2DToFile(this->zShift_cfr, "zshift_cfr");
    writeArray2DToFile(this->rxy_cfr, "rxy_cfr");
    writeArray2DToFile(this->zxy_cfr, "zxy_cfr");
    writeArray3DToFile(this->dxdy, "dxdy_0");
    writeArray3DToFile(this->dzdy, "dzdy_0");

//zs_cfr.
#if 0
        size_t numRows = this->ixsep, numCols = this->nypf2 - this->nypf1 + 1;
        this->zs_cfr.resize(numRows, std::vector<double>(numCols, 0.0));
        for (size_t i = 0; i < numRows; ++i)
            for (size_t j = 0; j < numCols; ++j)
                this->zs_cfr[i][j] = this->zShift[i][this->nypf1 + j];

        // Add an additional column to zs_cfr
        std::cout<<"Compute nu!!! "<<std::endl;
        for (size_t i = 0; i < numRows; ++i)
        {
            // Compute the value for the additional column
            double nuVal1 = 0.0; /* Replace with nu(xiarray_cfr[i], nypf1 + 1) */
            double nuVal2 = 0.0; /* Replace with nu(xiarray_cfr[i], nypf2) */
            double newVal = 0.5 * (nuVal1 + nuVal2) * dy0 + zs_cfr[i][numCols - 1];

            // Append the value to the row
            this->zs_cfr[i].push_back(newVal);
        }
#endif

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

  int GetRegion(const viskores::Vec3f& p) const
  {
    int region = -1;
    if (p[0] < static_cast<double>(this->ixsep + 0.5))
    {
      region = 0; //Closed flux surface
      std::cout << "Check this +1 stuff.." << std::endl;
      if (p[1] < this->nypf1 + 1 || p[1] > nypf2 - 1)
      {
        region = 2; //PFR
        std::cout << " We hit the +1 stuff..." << std::endl;
      }
    }
    else
      region = 1; //SOL

    // Check for divertor starting points
    if (this->direction == 1 && p[1] == this->ny - 1)
      region = 14;
    else if (this->direction == -1 && p[1] == 0)
      region = 13;

    return region;
  }

  int GetRegion(double xind, int yStart) const
  {
    int region = -1;
    if (xind < static_cast<double>(this->ixsep + 0.5))
    {
      region = 0; //Closed flux surface
      std::cout << "Check this +1 stuff.." << std::endl;
      if (yStart < this->nypf1 + 1 || yStart > nypf2 - 1)
      {
        region = 2; //PFR
        std::cout << " We hit the +1 stuff..." << std::endl;
      }
    }
    else
      region = 1; //SOL

    // Check for divertor starting points
    if (this->direction == 1 && yStart == this->ny - 1)
      region = 14;
    else if (this->direction == -1 && yStart == 0)
      region = 13;

    return region;
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
  int zperiod = 1;
  int jyseps1_1 = 16;
  int jyseps1_2 = -1;
  int jyseps2_1 = -1;
  int jyseps2_2 = 112;
  int ixseps1 = 195;
  int ixseps2 = 260;
  int ixsep = 195;
  int nypf1 = 16;
  int nypf2 = 112;

  int direction = 1;

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
};



double INTERP(const std::vector<double>& X, const std::vector<double>& Y, double val)
{
  // Ensure input vectors have the same size
  /*
    if (X.size() != Y.size())
        throw std::invalid_argument("X and Y must have the same size.");
    */

  // Find the appropriate interval for the value
  int idx = -1;
  for (size_t i = 0; i < X.size() - 1; ++i)
  {
    if (X[i] <= val && val <= X[i + 1])
    {
      idx = i;
      break;
    }
  }

  // Handle out-of-bounds value
  if (idx == -1)
    throw std::out_of_range("val is out of the interpolation range.");

  // Perform linear interpolation
  double x1 = X[idx];
  double x2 = X[idx + 1];
  double y1 = Y[idx];
  double y2 = Y[idx + 1];
  double y_i = y1 + (y2 - y1) * (val - x1) / (x2 - x1);

  return y_i;
}

void dumpTrajSamples(int iline, const std::vector<viskores::Vec3f>& points)
{
  int n = points.size();
  std::vector<double> tivals, xvals, yvals, zvals;

  for (int i = 0; i < n; i++)
  {
    tivals.push_back(double(i));
    xvals.push_back(points[i][0]);
    yvals.push_back(points[i][1]);
    zvals.push_back(points[i][2]);
  }

  //Replace with Hermite...
  SplineInterpolation splineX(tivals, xvals), splineY(tivals, yvals), splineZ(tivals, zvals);

  int region = -1;
  std::vector<double> tvals;
  double tmin = 0.0, tmax = double(n - 1);
  int id = 0;
  double x0 = -1.0;
  for (double t = tmin; t < tmax; t += 0.0005, id++)
  {
    double x = splineX.evaluate(t);
    double y = splineY.evaluate(t);
    double z = splineZ.evaluate(t);
    trajsplineFid << iline << ", " << t << ", " << x << ", " << y << ", " << z << ", " << region << std::endl;
    if (id > 0 && x * x0 < 0.0 && y > 0.0)
      puncFid << iline << ", " << t << ", " << x << ", " << y << ", " << z << std::endl;
    x0 = x;
  }
}

bool signChange(const double& x0, const double& x1)
{
  return (x0 * x1 < 0 && x0 != 0 && x1 != 0);
}

std::vector<viskores::Vec3f> ConvertToXYZSpace(const Options& opts, int iline, int id, int region, const std::vector<Point>& pts)
{
  std::vector<viskores::Vec3f> ptsXYZ;

  auto _fid = fopen("./problem.c.txt", "w");
  auto trajvals_fid = fopen("./trajvals.c.txt", "w");
  fprintf(trajvals_fid, "iter, xind, yend, zind, zend, x3d, y3d, z3d\n");

  for (const auto& pt : pts)
  {
    auto xind = pt.traj2;
    auto zind = pt.traj4;

    if (id == 58)
      std::cout << " ***** Issue." << std::endl;
    //double x = pt.x, y = pt.y, z = pt.z;
    //int xi = static_cast<int>(x), yi = static_cast<int>(y), zi = static_cast<int>(z);
    int yi = int(pt.traj3);

    // Get the column vector rxy[:, yi]
    std::vector<double> rxy_column(opts.rxy.size());
    for (size_t i = 0; i < opts.rxy.size(); ++i)
      rxy_column[i] = opts.rxy[i][yi];

    // Interpolate values
    double rxyvalue = INTERP(opts.xiarray, rxy_column, pt.traj2);

    std::vector<double> zShift_column(opts.zShift.size());
    for (size_t i = 0; i < opts.zShift.size(); ++i)
      zShift_column[i] = opts.zShift[i][yi];
    const auto& __xi = opts.xiarray;
    double zsvalue = INTERP(opts.xiarray, zShift_column, pt.traj2);

    double zvalue = INTERP(opts.ziarray, opts.zarray, pt.traj4);

    // Compute x3d_tmp and y3d_tmp
    double x3d_tmp = rxyvalue * std::cos(zsvalue);
    double y3d_tmp = rxyvalue * std::sin(zsvalue);

    // Compute x3d and y3d
    double x3d = x3d_tmp * std::cos(zvalue) - y3d_tmp * std::sin(zvalue);
    double y3d = x3d_tmp * std::sin(zvalue) + y3d_tmp * std::cos(zvalue);

    // Get the column vector zxy[:, yi]
    std::vector<double> zxy_column(opts.zxy.size());
    for (size_t i = 0; i < opts.zxy.size(); ++i)
      zxy_column[i] = opts.zxy[i][static_cast<size_t>(yi)];
    double z3d = INTERP(opts.xiarray, zxy_column, pt.traj2);

    trajOut << iline << ", " << id << ", " << x3d << ", " << y3d << ", " << z3d << ", " << region << ", " << yi << ", " << zsvalue << ", " << zvalue
            << std::endl;
    //fprintf(_fid, "*** istep= %d traj= %d %12.10e %d %12.10e\n", id, (int)pt.traj1, pt.traj2, (int)pt.traj3, pt.traj4);
    fprintf(_fid, "%d  traj4, zvalue= %12.10e %12.10e\n", id, pt.traj4, zvalue);
    //printf("*** istep= %d r= %12.10e zs= %12.10e z= %12.10e t23= %12.10e %d\n", id, rxyvalue, zsvalue, zvalue, pt.traj2, int(pt.traj3));
    ptsXYZ.push_back({ x3d, y3d, z3d });

    fprintf(trajvals_fid, "%d, %10.8f, %d, %10.8f, %10.8f, %10.8f, %10.8f, %10.8f\n", id, pt.traj2, (int)pt.traj3, pt.traj4, pt.traj7, x3d, y3d, z3d);

    id++;
  }

  fclose(trajvals_fid);
  dumpTrajSamples(iline, ptsXYZ);

  return ptsXYZ;
}

std::vector<viskores::Vec3f> FindPunctures(const std::vector<viskores::Id>& ids, const std::vector<viskores::Vec3f>& ptsXYZ, bool dumpPts)
{
  std::size_t n = ptsXYZ.size();
  viskores::FloatDefault t = 0.0;
  std::vector<viskores::FloatDefault> xVals, yVals, zVals, tVals;
  for (const auto& pt : ptsXYZ)
  {
    xVals.push_back(pt[0]);
    yVals.push_back(pt[1]);
    zVals.push_back(pt[2]);
    tVals.push_back(t);
    t += 1.0;
  }

  //SplineInterpolation splineX(tVals, xVals), splineY(tVals, yVals), splineZ(tVals, zVals);
  std::vector<viskores::Vec3f> punctures;
  Spline spl(ptsXYZ);

  for (std::size_t i = 1; i < n; i++)
  {
    const auto& p0 = ptsXYZ[i - 1];
    const auto& p1 = ptsXYZ[i];
    if (p0[0] * p1[0] > 0.0) // check for change sign in X.
      continue;

    if (p0[1] < 0.0 || p1[1] < 0.0)
      continue;

    //crosses the x=0 plane between i-1 and i.
    viskores::FloatDefault t0 = static_cast<viskores::FloatDefault>(i - 1);
    viskores::FloatDefault t1 = static_cast<viskores::FloatDefault>(i);

    auto x0 = spl.Evaluate(t0)[0];
    auto x1 = spl.Evaluate(t1)[0];
    int stop_idx = 0;
    for (int cnt = 0; cnt < 100; cnt++)
    {
      viskores::FloatDefault tMid = (t0 + t1) / 2.0;
      auto xMid = spl.Evaluate(tMid)[0];
      if (xMid * x0 < 0.0)
      {
        t1 = tMid;
        x1 = xMid;
      }
      else
      {
        t0 = tMid;
        x0 = xMid;
      }
      viskores::FloatDefault diff = viskores::Abs(t0 - t1);
      if (viskores::Abs(t0 - t1) < viskores::Epsilon<viskores::FloatDefault>())
        break;
      stop_idx++;
    }
    viskores::FloatDefault tVal = (t0 + t1) / 2.0;
    //viskores::Vec3f pt(splineX.evaluate(tVal), splineY.evaluate(tVal), splineZ.evaluate(tVal));
    auto pt = spl.Evaluate(tVal);
    if (viskores::Abs(pt[0]) > 10.0 * viskores::Epsilon<viskores::FloatDefault>())
    {
      //VISKORES_ASSERT(false);
      std::cout << "********* Intersection point didn't converge: " << pt[0] << " idx= " << stop_idx << " t: " << t0 << " " << t1 << std::endl;
    }
    if (dumpPts)
      puncSplineFid << ids[i] << ", " << tVal << ", " << pt[0] << ", " << pt[1] << ", " << pt[2] << std::endl;
    punctures.push_back(pt);
  }

  return punctures;
}
class EvaluateSplineWorklet : public viskores::worklet::WorkletMapField
{
public:
  using ControlSignature = void(ExecObject spline, FieldOut points, FieldOut params);
  using ExecutionSignature = void(InputIndex, _1, _2, _3);
  using InputDomain = _2;

  EvaluateSplineWorklet(viskores::FloatDefault tmin, viskores::FloatDefault tmax, viskores::Id n)
    : TMin(tmin)
    , TMax(tmax)
  {
    this->dT = (this->TMax - this->TMin) / static_cast<viskores::FloatDefault>(n - 1);
  }

  template <typename CubicSplineType, typename ResultType>
  VISKORES_EXEC void operator()(const viskores::Id& idx, const CubicSplineType& spline, ResultType& pos, viskores::FloatDefault& param) const
  {
    param = static_cast<viskores::FloatDefault>(idx) * this->dT;
    auto res = spline.Evaluate(param, pos);

    if (res != viskores::ErrorCode::Success)
      this->RaiseError("Spline evaluation failed.");
  }

private:
  viskores::FloatDefault TMin, TMax, dT;
};

class ComputePuncturesWorklet : public viskores::worklet::WorkletMapField
{
public:
  ComputePuncturesWorklet() = default;

  using ControlSignature = void(FieldIn t0, FieldIn t1, ExecObject cubicSpline, FieldOut param, FieldOut pos);
  using ExecutionSignature = void(_1, _2, _3, _4, _5);
  using InputDomain = _1;

  template <typename CubicSplineType, typename ResultType>
  VISKORES_EXEC void operator()(const viskores::FloatDefault& _t0,
                                const viskores::FloatDefault& _t1,
                                const CubicSplineType& spline,
                                viskores::FloatDefault& param,
                                ResultType& pos) const
  {
    viskores::FloatDefault t0 = _t0, t1 = _t1;
    viskores::Vec3f p0, p1;
    auto res = spline.Evaluate(t0, p0);
    if (res != viskores::ErrorCode::Success)
      this->RaiseError("Spline evaluation failed.");

    res = spline.Evaluate(t1, p1);
    if (res != viskores::ErrorCode::Success)
      this->RaiseError("Spline evaluation failed.");
    if (p0[0] * p1[0] > 0)
    {
#ifndef VISKORES_HIP
      std::cout << " t0/1= " << t0 << " " << t1 << " p0 " << p0[0] << " p1 " << p1[0] << std::endl;
#endif
      this->RaiseError("Points not on either side of puncture");
    }

    viskores::Vec3f pMid;
    viskores::FloatDefault tMid;
    for (int cnt = 0; cnt < 100; cnt++)
    {
      tMid = (t0 + t1) / 2.0;
      res = spline.Evaluate(tMid, pMid);
      if (res != viskores::ErrorCode::Success)
        this->RaiseError("Spline evaluation failed.");

      if (pMid[0] * p0[0] < 0.0)
      {
        t1 = tMid;
        p1 = pMid;
      }
      else
      {
        t0 = tMid;
        p0 = pMid;
      }
      viskores::FloatDefault diff = viskores::Abs(t0 - t1);
      if (viskores::Abs(t0 - t1) < viskores::Epsilon<viskores::FloatDefault>())
        break;
    }
    param = (t0 + t1) / 2.0f;
    res = spline.Evaluate(param, pos);
    if (res != viskores::ErrorCode::Success)
      this->RaiseError("Spline evaluation failed.");
    //std::cout << "wPUNC " << param << " " << pos << std::endl;
  }
};

struct SubOneFunctor
{
  VISKORES_EXEC viskores::FloatDefault operator()(viskores::FloatDefault x) const { return x - 1.0f; }
};


//./fieldline_tracing_vtkm 180 200 10 100
int main(int argc, char* argv[])
{
  bool doVTKm = true;

  double x0 = (double)std::stof(argv[1]);
  double x1 = (double)std::stof(argv[2]);
  int nlines = std::stoi(argv[3]);
  viskores::Id maxPuncs = std::stoi(argv[4]);
  std::string stuffFileName = argv[5];
  std::cout << "Running lines: " << x0 << " " << x1 << " #= " << nlines << std::endl;
  double dx = (x1 - x0) / (double)(nlines - 1);

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
  trajOut << "ID, STEP, X, Y, Z, REGION, YI, ZSVALUE, ZVALUE" << std::endl;
  tanOut << "ID, STEP, X, Y, Z, VX, VY, VZ" << std::endl;
  rawPunc << "ID, STEP, X, Y, Z\n";
  puncFid << "ID, STEP, X, Y, Z\n";
  puncFid2 << "ID, STEP, X, Y, Z\n";
  puncSplineFid << "ID, STEP, X, Y, Z\n";
  punc_ip_Fid << "ID, STEP, X, Y, Z\n";
  stepOut << "ID, STEP, X, Y, Z\n";
  rk4Out << "ID, STEP, X, Y, Z\n";
  auto TRAJ_FID = fopen("./pt_traj.c.txt", "w");
  fprintf(TRAJ_FID, "IT, xind, yEnd, zind, REG, zEnd\n");

  std::string fname = "/Users/dpn/proj/bout++/poincare/boutpp_poincare/poincare_clean/stuff.nc";
  //fname = "/lustre/orion/csc143/proj-shared/pugmire/stuff.nc";
  fname = "/Users/dpn/proj/bout++/poincare/boutpp_poincare/poincare_clean/python.27.Aug/stuff_python.nc";
  fname = "/Users/dpn/proj/bout++/poincare/boutpp_poincare/poincare_clean/utils/stuff.nc";

  fname = stuffFileName;

  bool transpose = true;
  if (stuffFileName.find("python") != std::string::npos)
    transpose = true;

  Options opts(fname, transpose);

  int divertor = 1; //single null
  double xind = 0.0f;

  std::vector<viskores::FloatDefault> LINES; // = {0, 50, 100, 150, 200, 250};
  for (viskores::FloatDefault x = x0; x < x1; x += dx)
    LINES.push_back(x);
  //LINES = { 0, 1, 2 };
  int nturns = 20;

  //LINES = {149};
  //LINES = { 0, 50, 100, 150, 200, 250 };
  //LINES = { 190 };
  //LINES = { 150.0, 150.1, 150.2, 150.3 };
  //LINES = { 50.0, 100.0, 150.0 };

  std::vector<Point> Points;

  //std::ofstream trajOut("traj.c.txt", std::ofstream::out);
  //trajOut<<"XI, ITER, X, Y, Z"<<std::endl;

  /*
    double _xStart = -0.135524, _yStart = 55, _zStart = 0;
    auto step0 = RK4_FLT1(_xStart, _yStart, _zStart, opts.dxdy, opts.dzdy, opts.xarray, opts.zarray, 0, opts.dxdy_p1, opts.dzdy_p1, 1, opts.nypf1, opts.nypf2);

    _yStart--;
    auto step2 = RK4_FLT1(_xStart, _yStart, _zStart, opts.dxdy, opts.dzdy, opts.xarray, opts.zarray, 0, opts.dxdy_p1, opts.dzdy_p1, 1, opts.nypf1, opts.nypf2);
*/

  auto start = std::chrono::high_resolution_clock::now();

  if (doVTKm)
  {
    std::vector<viskores::Vec3f> points;
    int idx = 0;
    for (const auto& iline : LINES)
    {
      xind = static_cast<viskores::FloatDefault>(iline);
      std::cout << "idx: " << idx++ << " xind= " << xind << std::endl;
      int yyy = opts.jyomp;
      int yStart = opts.jyomp;
      int zzz = 0;
      auto zStart = opts.zarray[zzz];
      //std::cout << "**** xind= " << xind << " opts.jyomp= " << opts.jyomp << " zStart= " << zStart << std::endl;
      std::cout << "xind= " << (int)xind << " jyomp= " << opts.jyomp << " dims: " << opts.psixy.size() << " " << opts.psixy[0].size() << std::endl;
      double xStart = opts.psixy[static_cast<int>(xind)][opts.jyomp];
      //double xStart2 = scalarField2DEval(opts.Grid2D, "psixy", viskores::Vec3f(xind, (double)yStart, 0));

      //int yind = yStart;
      //double zind = INTERP(opts.zarray, opts.ziarray, zStart);
      //auto zind2 = scalarField1DEval(opts.ZiArray, opts.ZArray, { zStart });

      //viskores::Particle p({ xStart, static_cast<viskores::FloatDefault>(yStart), zStart }, 0);
      viskores::Vec3f p0(xStart, static_cast<viskores::FloatDefault>(yStart), zStart);
      points.push_back(p0);
    }
    std::cout << "numPts= " << points.size() << std::endl;
    std::cout << "numPuncs= " << maxPuncs << std::endl;

    const auto& grid3D = opts.Grid3D;
    const auto& grid2D = opts.Grid2D;
    viskores::cont::CellLocatorRectilinearGrid locator2D, locator3D;
    locator3D.SetCoordinates(grid3D.GetCoordinateSystem());
    locator3D.SetCellSet(grid3D.GetCellSet());
    locator3D.Update();
    locator2D.SetCoordinates(grid2D.GetCoordinateSystem());
    locator2D.SetCellSet(grid2D.GetCellSet());
    locator2D.Update();
    viskores::Id maxSteps = maxPuncs * 20;

    RK4Worklet worklet(maxPuncs, maxSteps);
    worklet.ds3D = grid3D;
    worklet.ds2D = grid2D;
    worklet.grid3DBounds = grid3D.GetCoordinateSystem().GetBounds();
    worklet.grid2DBounds = grid2D.GetCoordinateSystem().GetBounds();
    worklet.nypf1 = opts.nypf1;
    worklet.nypf2 = opts.nypf2;
    worklet.ixsep1 = opts.ixseps1;
    worklet.ixsep2 = opts.ixseps2;
    BoutppField boutppField(opts.Grid3D, opts.Grid2D, opts.XiArray, opts.XArray, opts.YArray, opts.ZiArray, opts.ZArray, opts.ShiftAngle);

    viskores::cont::Invoker invoker;
    auto inPts = viskores::cont::make_ArrayHandle<viskores::Vec3f>(points, viskores::CopyFlag::On);
    auto xinds = viskores::cont::make_ArrayHandle<viskores::FloatDefault>(LINES, viskores::CopyFlag::On);
    viskores::cont::ArrayHandle<viskores::FloatDefault> dxdyField, dzdyField, rxyField, zShiftField;
    viskores::cont::ArrayHandle<viskores::Vec3f> result, tangent;
    grid3D.GetField("dxdy").GetData().AsArrayHandle<viskores::FloatDefault>(dxdyField);
    grid3D.GetField("dzdy").GetData().AsArrayHandle<viskores::FloatDefault>(dzdyField);
    grid2D.GetField("rxy").GetData().AsArrayHandle<viskores::FloatDefault>(rxyField);
    grid2D.GetField("zShift").GetData().AsArrayHandle<viskores::FloatDefault>(zShiftField);


    auto start = std::chrono::high_resolution_clock::now();

    viskores::cont::ArrayHandle<viskores::Id> puncIndices, pointIds;
    viskores::cont::ArrayHandle<bool> validSteps, validPuncs;
    result.Allocate(inPts.GetNumberOfValues() * maxSteps);
    tangent.Allocate(inPts.GetNumberOfValues() * maxSteps);
    validSteps.AllocateAndFill(inPts.GetNumberOfValues() * maxSteps, false);
    validPuncs.AllocateAndFill(inPts.GetNumberOfValues() * maxPuncs, false);
    puncIndices.Allocate(inPts.GetNumberOfValues() * maxPuncs);
    pointIds.AllocateAndFill(inPts.GetNumberOfValues() * maxSteps, -1);
    //invoker(worklet, inPts, locator3D, grid3D.GetCellSet(), locator2D, grid2D.GetCellSet(), dxdyField, dzdyField, rxyField, zShiftField, puncIndices, result);

    invoker(worklet, inPts, boutppField, grid3D.GetCellSet(), grid2D.GetCellSet(), puncIndices, pointIds, result, tangent, validPuncs, validSteps);

    viskores::cont::ArrayHandle<viskores::Vec3f> validResult, validTangent;
    viskores::cont::ArrayHandle<viskores::Id> validPuncIndices, validPointIds;
    viskores::cont::Algorithm::CopyIf(result, validSteps, validResult);
    viskores::cont::Algorithm::CopyIf(tangent, validSteps, validTangent);
    viskores::cont::Algorithm::CopyIf(puncIndices, validPuncs, validPuncIndices);
    viskores::cont::Algorithm::CopyIf(pointIds, validSteps, validPointIds);
    //viskores::cont::printSummary_ArrayHandle(validPuncIndices, std::cout, true);
    //viskores::cont::printSummary_ArrayHandle(validResult, std::cout, true);
    // do with splines.
    std::vector<viskores::Vec3f> _points;
    std::vector<viskores::Id> _ids;
    for (viskores::Id i = 0; i < validResult.GetNumberOfValues(); i++)
    {
      auto pt = validResult.ReadPortal().Get(i);
      _points.push_back(pt);
      _ids.push_back(validPointIds.ReadPortal().Get(i));
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;
    std::cout << " NUM particles= " << LINES.size() << " numPunc= " << maxPuncs << std::endl;
    auto puncs = FindPunctures(_ids, _points, true);
    return 0;


    viskores::cont::CubicHermiteSpline trajSpline;

    //compute punctures.
    viskores::cont::printSummary_ArrayHandle(validPuncIndices, std::cout, true);

    //Puncture occors between puncIndices and puncIndices+1.
    //auto puncParams1 = viskores::cont::make_ArrayHandleCast<viskores::FloatDefault>(validPuncIndices);
    //auto puncParams0 = viskores::cont::make_ArrayHandleTransform(puncParams1, SubOneFunctor{});
    viskores::cont::ArrayHandle<viskores::FloatDefault> puncParams0, puncParams1;
    puncParams0.AllocateAndFill(validPuncIndices.GetNumberOfValues(), -1.0f);
    puncParams1.AllocateAndFill(validPuncIndices.GetNumberOfValues(), -1.0f);



    //viskores::cont::printSummary_ArrayHandle(puncParams0, std::cout, true);
    //viskores::cont::printSummary_ArrayHandle(puncParams1, std::cout, true);

    //std::cout << "rs= " << validResult.GetNumberOfValues() << std::endl;
    std::vector<viskores::FloatDefault> knots(validResult.GetNumberOfValues());
    std::iota(knots.begin(), knots.end(), 0.0f);
    trajSpline.SetData(validResult);
    //trajSpline.SetKnots(knots);
    //trajSpline.SetTangents(validTangent);
    viskores::cont::ArrayHandle<viskores::Vec3f> puncturePoints;
    viskores::cont::ArrayHandle<viskores::FloatDefault> punctureParams;

    auto _knots = trajSpline.GetKnots().ReadPortal();
    auto _nk = _knots.GetNumberOfValues();
    auto _nn = validPuncIndices.GetNumberOfValues();
    auto _idxPortal = validPuncIndices.ReadPortal();
    auto _portal0 = puncParams0.WritePortal();
    auto _portal1 = puncParams1.WritePortal();
    for (viskores::Id i = 0; i < _idxPortal.GetNumberOfValues(); i++)
    {
      viskores::Id idx = _idxPortal.Get(i);
      if (idx == 0)
        continue;
      //_portal0.Set(idx, _knots.Get(idx - 1));
      _portal0.Set(i, _knots.Get(idx - 1));
      _portal1.Set(i, _knots.Get(idx));
      auto k0 = _knots.Get(idx - 1);
      auto k1 = _knots.Get(idx);
      std::cout << "i: " << i << " idx= " << idx << " " << k0 << " " << k1 << std::endl;
      if (k0 > k1)
        std::cout << "**** ERROR: k0 > k1" << std::endl;
    }

    viskores::cont::printSummary_ArrayHandle(puncParams0, std::cout, true);
    viskores::cont::printSummary_ArrayHandle(puncParams1, std::cout, true);


    invoker(ComputePuncturesWorklet{}, puncParams0, puncParams1, trajSpline, punctureParams, puncturePoints);


    auto portal = puncturePoints.ReadPortal();
    auto portalp = punctureParams.ReadPortal();
    auto portali = validPuncIndices.ReadPortal();
    auto portalt = puncParams0.ReadPortal();
    for (viskores::Id i = 0; i < portal.GetNumberOfValues(); i++)
    {
      auto pt = portal.Get(i);
      puncSplineFid << 0 << ", " << portalp.Get(i) << ", " << pt[0] << ", " << pt[1] << ", " << pt[2] << std::endl;
      //std::cout << " PUNC: " << i << " " << portalt.Get(i) << " " << portali.Get(i) << " " << portal.Get(i) << std::endl;
    }

    //evaluate the points...
    viskores::Id nEval = 50000;
    auto k0 = trajSpline.GetKnots().ReadPortal().Get(0);
    auto _n = trajSpline.GetKnots().GetNumberOfValues();
    auto k1 = trajSpline.GetKnots().ReadPortal().Get(_n - 1);
    //auto k0 = knots[0];
    //auto k1 = knots[knots.size() - 1];
    //auto _n = knots.size();
    EvaluateSplineWorklet evalWorklet(k0, k1, nEval);
    viskores::cont::ArrayHandle<viskores::Vec3f> evalPoints;
    viskores::cont::ArrayHandle<viskores::FloatDefault> evalParam;
    evalPoints.Allocate(nEval);
    invoker(evalWorklet, trajSpline, evalPoints, evalParam);
    portal = evalPoints.ReadPortal();
    auto portal_p = evalParam.ReadPortal();
    for (viskores::Id i = 0; i < portal.GetNumberOfValues(); i++)
    {
      auto pt = portal.Get(i);
      trajsplineFid << "0, " << portal_p.Get(i) << ", " << pt[0] << ", " << pt[1] << ", " << pt[2] << ", 0 " << std::endl;
      std::cout << " dumping: " << pt << std::endl;
    }

    //    auto end = std::chrono::high_resolution_clock::now();
    //    std::chrono::duration<double> elapsed = end - start;
    //    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
  }

  for (const auto& iline : LINES)
  {
    xind = static_cast<double>(iline);
    int yyy = opts.jyomp;
    int yStart = opts.jyomp;
    int zzz = 0;
    auto zStart = opts.zarray[zzz];
    std::cout << "**** xind= " << xind << " opts.jyomp= " << opts.jyomp << " zStart= " << zStart << std::endl;
    double xStart = opts.psixy[static_cast<int>(xind)][opts.jyomp];
    double xStart2 = scalarField2DEval(opts.Grid2D, "psixy", viskores::Vec3f(xind, (double)yStart, 0));

    //int yind = yStart;
    double zind = INTERP(opts.zarray, opts.ziarray, zStart);
    auto zind2 = scalarField1DEval(opts.ZiArray, opts.ZArray, { zStart });

    //viskores::Particle p({ xStart, static_cast<viskores::FloatDefault>(yStart), zStart }, 0);
    viskores::Vec3f p0(xStart, static_cast<viskores::FloatDefault>(yStart), zStart);

#ifndef VISKORES_HIP
    std::cout << std::setprecision(12);
    std::cout << "   xStart=  " << xStart << std::endl;
    std::cout << "   xStart2= " << xStart2 << std::endl;
    std::cout << "   zind= " << zind << std::endl;
    std::cout << "   zind2= " << zind << std::endl;
#endif

    int region = opts.GetRegion(xind, p0[1]);
    int region_ = opts.GetRegion(p0);
    int iturn = 0, it = 0;

    auto zindFID = fopen("./zind.c.txt", "w");
    std::cout << "Region= " << region << std::endl;
    int step = 0;
    while (region < 10 && iturn < nturns)
    {
      // Start field-line tracing.
      // iy is not used -- just a loop...
      for (int iy_ = 0; iy_ < opts.ny - 1; iy_++)
      {
        //if (iturn > 5) break;

        //trajOut<<iline<<", "<<iy<<", "<<it<<", "<<iturn<<", "<<xStart<<", "<<p0[1]<<", "<<p0[2]<<std::endl;
        if (it == 0)
        {
          Point _p;
          _p.traj1 = it;
          _p.traj2 = xind;
          _p.traj3 = p0[1];
          _p.traj4 = zind;
          _p.traj5 = region;
          _p.traj7 = p0[2];
          fprintf(TRAJ_FID, "%d, %12.8f, %d, %12.8f, %d, %12.8f\n", (int)_p.traj1, _p.traj2, (int)_p.traj3, _p.traj4, (int)_p.traj5, _p.traj7);
          Points.push_back(_p);
          //Points.push_back({xind, p0[1], zind, iline, iy, it});
          step++;
        }

        if (p0[1] + 1 == opts.dxdy[0].size())
        {
          std::cout << "Overflow of some kind... Need to track this down." << std::endl;
          break;
        }

        if (step == 58)
        {
          std::cout << "Hey man, stop here...." << std::endl;
        }

        //double xEnd, zEnd, yEnd;
        viskores::Vec3f p1;
        if (region == 0 && p0[1] >= opts.nypf1 && p0[1] < opts.nypf2 + 1)
        {
          bool dumpFiles = false;
          //auto step = RK4_FLT1(xStart, p0[1], p0[2], opts.dxdy, opts.dzdy, opts.xarray, opts.zarray, region, opts.dxdy_p1, opts.dzdy_p1, 1, opts.nypf1, opts.nypf2, rk4Out, iline, it, dumpFiles);
          std::cout << "Begin Step: iturn= " << iturn << " iy= " << iy_ << std::endl;
          p1 = RK4_FLT1_vtkm(p0,
                             opts.Grid2D,
                             opts.Grid2D_cfr,
                             opts.Grid2D_xz,
                             opts.Grid3D,
                             opts.XArray,
                             opts.ZArray,
                             region,
                             1,
                             opts.nypf1,
                             opts.nypf2,
                             rk4Out,
                             iline,
                             it,
                             dumpFiles);
          p1[1] += 1.0;
        }
        std::cout << "  *** vRK4: end= " << p1 << std::endl;
        stepOut << iline << ", " << it << ", " << p1[0] << ", " << p1[1] << ", " << p1[2] << std::endl;


        // Check where the field line ends
        if (p1[0] > opts.xMax)
        {
          std::cout << "  Starting xind= " << xind << " line= " << iline << " reaches outer boundary" << std::endl;
          region = 12;
        }
        else if (p1[0] < opts.xMin)
        {
          std::cout << "  Starting xind= " << xind << " line= " << iline << " reaches inner boundary" << std::endl;
          region = 11;
        }
        else
        {
          //xind = INTERP(opts.xarray, opts.xiarray, xEnd);
          // interpoloated index from the x value (p1[0]).
          xind = scalarField1DEval(opts.XArray, opts.XiArray, p1[0]);
          std::cout << "XIND: " << p1 << " --> " << xind << std::endl;
          if (xind > static_cast<double>(opts.ixseps1) + 0.5)
          {
            region = 1;
            std::cout << "  Starting xind= " << xind << " line= " << iline << " enters the SOL." << std::endl;
          }
        }

        //Twist-shift at branch cut.
        if (p0[1] == opts.nypf2 - 1 && region == 0)
        {
          std::cout << "Branch cut: " << p0[1] << " " << opts.nypf2 << std::endl;
          //double shiftAngle = INTERP(opts.xiarray, opts.shiftAngle, xind);
          double shiftAngle = scalarField1DEval(opts.XiArray, opts.ShiftAngle, xind);
          p1[2] = p1[2] + shiftAngle;
          p1[1] = opts.nypf1;
        }

        Point _p;
        _p.traj1 = it;
        _p.traj2 = xind;
        _p.traj3 = p1[1];
        _p.traj5 = region;
        _p.traj7 = p1[2];

        double zEnd_no_mod = p1[2];
        //Relabel toroidal location.
        if (p1[2] < opts.zmin || p1[2] > opts.zmax)
          p1[2] = double_mod(p1[2], opts.zmax);
        //zind = INTERP(opts.zarray, opts.ziarray, zEnd);
        zind = scalarField1DEval(opts.ZArray, opts.ZiArray, { p1[2] });
        _p.traj4 = zind;
        if (step == 58)
        {
          std::cout << " ****** issue" << std::endl;
          ConvertToXYZSpace(opts, iline, step, region, { _p });
        }
        Points.push_back(_p);
        std::cout << "x/zind: " << p1 << " --> " << xind << " " << zind << std::endl;

        //std::cout<<"********** it= "<<it<<" pt1= "<<xEnd<<" "<<yEnd<<" "<<zEnd<<" zind "<<zind<<std::endl;
        fprintf(zindFID, "%d %12.10f --> %12.10f\n", it, p1[2], zind);
        //Points.push_back({xind, yEnd, zind, iline, iy, it});
        if (p1[2] > opts.zmax)
        {
          std::cout << "We have a problem now..." << std::endl;
          double diff = std::abs(p1[2] - opts.zmax);
          std::cout << "diff is " << diff << std::endl;
          std::cout << "*******" << std::endl;
          throw std::runtime_error("Meow");
        }

        it = it + 1;
        p0 = p1;
        //throw std::runtime_error("Meow");

        fprintf(TRAJ_FID, "%d, %12.8f, %d, %12.8f, %d, %12.8f\n", it, _p.traj2, (int)_p.traj3, _p.traj4, (int)_p.traj5, _p.traj7);
        step++;
      }
      iturn++;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

    //Convert to XYZ space.
    int id = 0;
    auto PointsXYZ = ConvertToXYZSpace(opts, iline, id, region, Points);
    std::vector<viskores::Id> _lines;
    _lines.push_back(iline);
    auto punctures = FindPunctures(_lines, PointsXYZ, true);
    return 0;
  }
}
