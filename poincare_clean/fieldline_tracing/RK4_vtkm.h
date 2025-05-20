#pragma once

#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/CellLocatorRectilinearGrid.h>
#include <viskores/cont/DataSet.h>
#include <viskores/exec/CellInterpolate.h>
#include <viskores/worklet/WorkletMapField.h>

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

void writeArray1DToFile(std::vector<double>& array, const std::string& fname);

void writeArray2DToFile(std::vector<std::vector<double>>& array, const std::string& fname);

void writeArray3DToFile(const std::vector<std::vector<std::vector<double>>>& array, const std::string& fname);

viskores::Vec3f RK4_FLT1_vtkm(const viskores::Vec3f& pStart,
                              const viskores::cont::DataSet& grid2D,
                              const viskores::cont::DataSet& grid2D_cfr,
                              const viskores::cont::DataSet& grid2D_xz,
                              const viskores::cont::DataSet& grid3D,
                              const viskores::cont::ArrayHandle<viskores::FloatDefault>& xarray,
                              const viskores::cont::ArrayHandle<viskores::FloatDefault>& zarray,
                              int region,
                              int direction,
                              int nypf1,
                              int nypf2,
                              std::ofstream& rk4Out,
                              int iline,
                              int it,
                              bool dumpFiles);

class BoutppFieldExecutionObject
{
  using LocatorType = viskores::exec::CellLocatorRectilinearGrid;
  using ArrayType = viskores::cont::ArrayHandle<viskores::FloatDefault>;
  using ArrayPortalType = typename ArrayType::ReadPortalType;

public:
  BoutppFieldExecutionObject(const LocatorType& locator3D,
                             const LocatorType& locator2D,
                             const ArrayPortalType& _dxdy,
                             const ArrayPortalType& _dzdy,
                             const ArrayPortalType& _rxy,
                             const ArrayPortalType& _zxy,
                             const ArrayPortalType& _zshift,
                             const ArrayPortalType& xiarray,
                             const ArrayPortalType& xarray,
                             const ArrayPortalType& yarray,
                             const ArrayPortalType& ziarray,
                             const ArrayPortalType& zarray,
                             const ArrayPortalType& shiftAngle)

    : Locator3D(locator3D)
    , Locator2D(locator2D)
    , dxdy(_dxdy)
    , dzdy(_dzdy)
    , rxy(_rxy)
    , zxy(_zxy)
    , zShift(_zshift)
    , XiArray(xiarray)
    , XArray(xarray)
    , YArray(yarray)
    , ZiArray(ziarray)
    , ZArray(zarray)
    , ShiftAngle(shiftAngle)
  {
  }

  LocatorType Locator2D, Locator3D;
  ArrayPortalType dxdy, dzdy;
  ArrayPortalType rxy, zxy, zShift;
  ArrayPortalType XiArray, XArray, YArray, ZiArray, ZArray, ShiftAngle;
};

class BoutppField : public viskores::cont::ExecutionObjectBase
{
  using ArrayType = viskores::cont::ArrayHandle<viskores::FloatDefault>;
  using Structured2DType = viskores::cont::CellSetStructured<2>;
  using Structured3DType = viskores::cont::CellSetStructured<3>;
  //using Structured3DType = viskores::exec::CellSetStructured<3>;
  //using Structured2DType = viskores::exec::CellSetStructured<2>;
  //using CellSet2DExecType = viskores::exec::CellSet<Structured2DType>;
  //using CellSet3DExecType = viskores::exec::CellSet<Structured3DType>;
  //using X2Type = viskores::exec::ConnectivityStructured<viskores::TopologyElementTagCell{}, viskores::TopologyElementTagPoint{}, 2>;
  //using X3Type = viskores::exec::ConnectivityStructured<viskores::TopologyElementTagCell{}, viskores::TopologyElementTagPoint{}, 3>;


public:
  BoutppField(const viskores::cont::DataSet& grid3D,
              const viskores::cont::DataSet& grid2D,
              const ArrayType& xiarray,
              const ArrayType& xarray,
              const ArrayType& yarray,
              const ArrayType& ziarray,
              const ArrayType& zarray,
              const ArrayType& shiftAngle)
    : Grid3D(grid3D)
    , Grid2D(grid2D)
    , XiArray(xiarray)
    , XArray(xarray)
    , YArray(yarray)
    , ZiArray(ziarray)
    , ZArray(zarray)
    , ShiftAngle(shiftAngle)
  {
  }

  VISKORES_CONT BoutppFieldExecutionObject PrepareForExecution(viskores::cont::DeviceAdapterId device, viskores::cont::Token& token) const
  {
    viskores::cont::CellLocatorRectilinearGrid locator3D;
    viskores::cont::CellLocatorRectilinearGrid locator2D;
    locator3D.SetCoordinates(this->Grid3D.GetCoordinateSystem());
    locator3D.SetCellSet(this->Grid3D.GetCellSet());
    locator2D.SetCoordinates(this->Grid2D.GetCoordinateSystem());
    locator2D.SetCellSet(this->Grid2D.GetCellSet());
    locator3D.Update();
    locator2D.Update();

    //get the fields.
    viskores::cont::ArrayHandle<viskores::FloatDefault> dxdyField, dzdyField, rxyField, zxyField, zShiftField;
    this->Grid3D.GetField("dxdy").GetData().AsArrayHandle<viskores::FloatDefault>(dxdyField);
    this->Grid3D.GetField("dzdy").GetData().AsArrayHandle<viskores::FloatDefault>(dzdyField);
    this->Grid2D.GetField("rxy").GetData().AsArrayHandle<viskores::FloatDefault>(rxyField);
    this->Grid2D.GetField("zxy").GetData().AsArrayHandle<viskores::FloatDefault>(zxyField);
    this->Grid2D.GetField("zShift").GetData().AsArrayHandle<viskores::FloatDefault>(zShiftField);

    return BoutppFieldExecutionObject(locator3D.PrepareForExecution(device, token),
                                      locator2D.PrepareForExecution(device, token),
                                      dxdyField.PrepareForInput(device, token),
                                      dzdyField.PrepareForInput(device, token),
                                      rxyField.PrepareForInput(device, token),
                                      zxyField.PrepareForInput(device, token),
                                      zShiftField.PrepareForInput(device, token),
                                      this->XiArray.PrepareForInput(device, token),
                                      this->XArray.PrepareForInput(device, token),
                                      this->YArray.PrepareForInput(device, token),
                                      this->ZiArray.PrepareForInput(device, token),
                                      this->ZArray.PrepareForInput(device, token),
                                      this->ShiftAngle.PrepareForInput(device, token));
  }

  viskores::cont::DataSet Grid3D;
  viskores::cont::DataSet Grid2D;
  ArrayType XiArray, XArray;
  ArrayType YArray;
  ArrayType ZiArray, ZArray;
  ArrayType ShiftAngle;
};

extern std::ofstream stepOut, trajOut, tanOut;
class RK4Worklet : public viskores::worklet::WorkletMapField
{
public:
  RK4Worklet(viskores::Id maxPuncs, viskores::Id maxSteps)
    : MaxPunctures(maxPuncs)
    , MaxSteps(maxSteps)
  {
  }

  viskores::FloatDefault xMin, xMax;
  viskores::Bounds grid3DBounds, grid2DBounds;
  int nypf1, nypf2;
  int ixsep1, ixsep2;

  using ControlSignature = void(FieldIn points,
                                ExecObject boutppField,
                                WholeCellSetIn<> cells3D,
                                WholeCellSetIn<> cells2D,
                                WholeArrayOut puncIndices,
                                WholeArrayOut pointId,
                                WholeArrayOut result,
                                WholeArrayOut tangent,
                                WholeArrayInOut validPuncs,
                                WholeArrayInOut validSteps);
  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10);
  using InputDomain = _1;

  template <typename BoutppFieldType,
            typename CellSet3DType,
            typename CellSet2DType,
            typename IndexType,
            typename ResultType,
            typename ValidPuncType,
            typename ValidStepType>
  VISKORES_EXEC void operator()(const viskores::Id& idx,
                                const viskores::Vec3f& pStart,
                                BoutppFieldType& boutppField,
                                const CellSet3DType& cells3D,
                                const CellSet2DType& cells2D,
                                IndexType& puncIndex,
                                IndexType& pointId,
                                ResultType& result,
                                ResultType& tangent,
                                ValidPuncType& validPuncs,
                                ValidStepType& validSteps) const
  {
    constexpr double twoPi = 2.0 * M_PI;
    viskores::Id numPunc = 0;
    viskores::Id numSteps = 0;
    auto p0 = pStart;

    viskores::Id stepOffset = idx * this->MaxSteps;
    viskores::Id puncOffset = idx * this->MaxPunctures;
    result.Set(stepOffset, p0);
    //stepOut << idx << ", " << 0 << ", " << p0[0] << ", " << p0[1] << ", " << p0[2] << std::endl;
    int region = 0;

    viskores::Vec3f pCart0, pCart1, vCart0;

    //pCart0 = this->ConvertToCartesian(p0, boutppField, cells2D, 0, true);
    pCart0 = this->ConvertToCartesian(p0, vCart0, boutppField, cells2D, 0, true, vCart0);
    //vCart0 = this->EvaluateTangent(p0, boutppField, cells2D);
    tangent.Set(stepOffset, vCart0);
    tanOut << "0, 0, " << pCart0[0] << ", " << pCart0[1] << ", " << pCart0[2] << ", " << vCart0[0] << ", " << vCart0[1] << ", " << vCart0[2]
           << std::endl;

    for (viskores::Id step = 1; region == 0 && step < this->MaxSteps; step++)
    {
      viskores::Vec3f p1, tan0, tan1;
      this->RK4Step(p0, boutppField, cells3D, p1, tan0, tan1);
      auto xind = this->LinearInterpolate(boutppField.XArray, boutppField.XiArray, p1[0]);
      if (p1[0] > this->grid2DBounds.X.Max)
      {
        region = 12;
        std::cout << "p1= " << p1 << " *** region 12" << std::endl;
        VISKORES_ASSERT(false);
      }
      else if (p1[0] < this->grid2DBounds.X.Min)
      {
        region = 11;
        std::cout << "p1= " << p1 << " *** region 11" << std::endl;
        VISKORES_ASSERT(false);
      }
      else
      {
        if (xind > static_cast<viskores::FloatDefault>(this->ixsep1) + 0.5)
          region = 1;
        //        std::cout << "p1= " << p1 << std::endl;
        //std::cout << "Region 1: xind= " << xind << std::endl;
        //VISKORES_ASSERT(false);
        /*
                xind = scalarField1DEval(opts.XArray, opts.XiArray, p1[0]);
                std::cout<<"   vINTERP:  xind "<<xind<<std::endl;
                if (xind > static_cast<double>(opts.ixsep1) + 0.5)
                    region = 1;
                */
      }

      //do zshift stuff.
      //Twist-shift at branch cut.
      if (p0[1] == this->nypf2 - 1 && region == 0)
      {
        auto sa = this->LinearInterpolate(boutppField.XiArray, boutppField.ShiftAngle, xind);
        //double shiftAngle = 0; //scalarField1DEval(opts.XiArray, opts.ShiftAngle, xind);
        p1[2] = p1[2] + sa;
        p1[1] = this->nypf1;
      }

      p1[2] = this->floatMod(p1[2], twoPi);
      stepOut << idx << ", " << step << ", " << p1[0] << ", " << p1[1] << ", " << p1[2] << std::endl;
      //convert to cartesian.
      pCart0 = pCart1;
      viskores::Vec3f vCart1;
      pCart1 = this->ConvertToCartesian(p1, tan1, boutppField, cells2D, step, true, vCart1);

      //trajOut << "0, " << step << " " << tmp[0] << ", " << tmp[1] << ", " << tmp[2] << ", 0, 0, 0, 0" << std::endl;

      result.Set(stepOffset + step, pCart1);
      //vCart1 = this->EvaluateTangent(p1, boutppField, cells2D);
      tangent.Set(stepOffset + step, vCart1);
      tanOut << "0, " << step << ", " << pCart1[0] << ", " << pCart1[1] << ", " << pCart1[2] << ", " << vCart1[0] << ", " << vCart1[1] << ", "
             << vCart1[2] << std::endl;
      validSteps.Set(stepOffset + step, true);
      pointId.Set(stepOffset + step, idx);

      //point crosses the X=0 plane.
      if (step > 1 && pCart0[0] * pCart1[0] < 0)
      {
        //std::cout << "SET_PUNC: " << puncOffset << " " << numPunc << " " << (step - 1) << std::endl;
        puncIndex.Set(puncOffset + numPunc, step - 1);
        validPuncs.Set(puncOffset + numPunc, true);
        numPunc++;
      }

      p0 = p1;
      if (numPunc == this->MaxPunctures)
        break;
    }
  }

private:
  //x/zind: [-0.135435384332,55,0.000486567073233] --> 149.790391597 0.019747086487

  template <typename BoutppFieldType, typename CellSetType>
  viskores::Vec3f EvaluateTangent(const viskores::Vec3f& pt0, BoutppFieldType& boutppField, const CellSetType& cells) const
  {
    constexpr double twoPi = 2.0 * M_PI;
    auto pt0Cart = this->ConvertToCartesian(pt0, boutppField, cells, -1, false);

    auto vec2d = this->Evaluate(pt0, boutppField, cells);
    viskores::Vec3f vec0(vec2d[0], 1, vec2d[1]);
    auto pt1 = pt0 + vec0;
    pt1[2] = this->floatMod(pt1[2], twoPi);
    auto pt1Cart = this->ConvertToCartesian(pt1, boutppField, cells, -1, false);

    auto tangent = pt1Cart - pt0Cart;
    tangent = tangent * 0.75;
    return tangent;
  }

  template <typename BoutppFieldType, typename CellSetType>
  viskores::Vec3f ConvertToCartesian(const viskores::Vec3f& pt,
                                     const viskores::Vec3f& vec,
                                     BoutppFieldType& boutppField,
                                     const CellSetType& cells,
                                     int step,
                                     bool outputIt,
                                     viskores::Vec3f& vecCart) const
  {
    viskores::Id yi = static_cast<viskores::Id>(pt[1]);
    auto xind = this->LinearInterpolate(boutppField.XArray, boutppField.XiArray, pt[0]);
    auto zind = this->LinearInterpolate(boutppField.ZArray, boutppField.ZiArray, pt[2]);

    viskores::Vec3f ptXY(pt[0], pt[1], 0.0);
    auto rxy = this->Evaluate(ptXY, boutppField, boutppField.rxy, cells);
    auto zxy = this->Evaluate(ptXY, boutppField, boutppField.zxy, cells);
    auto zvalue = this->LinearInterpolate(boutppField.ZiArray, boutppField.ZArray, zind);
    auto zsvalue = this->Evaluate(ptXY, boutppField, boutppField.zShift, cells);

    // Compute x3d and y3d
    viskores::FloatDefault x3d = rxy * std::cos(zsvalue) * std::cos(zvalue) - rxy * std::sin(zsvalue) * std::sin(zvalue);
    viskores::FloatDefault y3d = rxy * std::cos(zsvalue) * std::sin(zvalue) + rxy * std::sin(zsvalue) * std::cos(zvalue);
    viskores::FloatDefault z3d = zxy;
    if (outputIt)
      trajOut << "0, " << step << ", " << x3d << ", " << y3d << ", " << z3d << ", 0, " << yi << ", " << zsvalue << ", " << zvalue << std::endl;

    // Compute the tangent vector if requested
    if (true)
    {
      viskores::FloatDefault vr = vec[0];     // Radial component
      viskores::FloatDefault vtheta = vec[1]; // Angular component
      viskores::FloatDefault vz = vec[2];     // Z component

      vtheta *= rxy; // This is the critical fix

      // Compute Cartesian vector components
      viskores::FloatDefault cos_z = std::cos(zvalue);
      viskores::FloatDefault sin_z = std::sin(zvalue);
      viskores::FloatDefault cos_zs = std::cos(zsvalue);
      viskores::FloatDefault sin_zs = std::sin(zsvalue);
      viskores::FloatDefault cos_zs_z = std::cos(zsvalue + zvalue);
      viskores::FloatDefault sin_zs_z = std::sin(zsvalue + zvalue);

      vecCart[0] = vr * cos_zs * cos_z - vr * sin_zs * sin_z - rxy * vtheta * sin_zs_z;
      vecCart[1] = vr * cos_zs * sin_z + vr * sin_zs * cos_z + rxy * vtheta * cos_zs_z;
      vecCart[2] = vz; // Z component is directly mapped
      //vecCart = vecCart * 0.50;
      //std::cout << " *********** vecCart= " << vecCart << std::endl;
    }

    return viskores::Vec3f(x3d, y3d, z3d);
  }

  template <typename BoutppFieldType, typename CellSetType>
  viskores::Vec3f ConvertToCartesianWithTangent(const viskores::Vec3f& pt,
                                                const viskores::Vec3f& tangent,
                                                BoutppFieldType& boutppField,
                                                const CellSetType& cells,
                                                viskores::Vec3f& tangentCartian,
                                                int step,
                                                bool outputIt) const
  {
    constexpr double twoPi = 2.0 * M_PI;
    //return this->ConvertToCartesian(pt, boutppField, cells, step, outputIt);

    viskores::Id yi = static_cast<viskores::Id>(pt[1]);
    auto xind = this->LinearInterpolate(boutppField.XArray, boutppField.XiArray, pt[0]);
    auto zind = this->LinearInterpolate(boutppField.ZArray, boutppField.ZiArray, pt[2]);

    viskores::Vec3f ptXY(pt[0], pt[1], 0.0);
    auto rxy = this->Evaluate(ptXY, boutppField, boutppField.rxy, cells);
    auto zxy = this->Evaluate(ptXY, boutppField, boutppField.zxy, cells);
    auto zvalue = this->LinearInterpolate(boutppField.ZiArray, boutppField.ZArray, zind);
    auto zsvalue = this->Evaluate(ptXY, boutppField, boutppField.zShift, cells);

    // Compute Cartesian coordinates
    viskores::FloatDefault x3d = rxy * std::cos(zsvalue) * std::cos(zvalue) - rxy * std::sin(zsvalue) * std::sin(zvalue);
    viskores::FloatDefault y3d = rxy * std::cos(zsvalue) * std::sin(zvalue) + rxy * std::sin(zsvalue) * std::cos(zvalue);
    viskores::FloatDefault z3d = zxy;

    // Jacobian matrix elements for tangent transformation
    viskores::FloatDefault dx_dr = std::cos(zsvalue) * std::cos(zvalue) - std::sin(zsvalue) * std::sin(zvalue);
    viskores::FloatDefault dx_dzs = -rxy * std::sin(zsvalue) * std::cos(zvalue) - rxy * std::cos(zsvalue) * std::sin(zvalue);
    viskores::FloatDefault dx_dz = -rxy * std::cos(zsvalue) * std::sin(zvalue) - rxy * std::sin(zsvalue) * std::cos(zvalue);

    viskores::FloatDefault dy_dr = std::sin(zsvalue) * std::cos(zvalue) + std::cos(zsvalue) * std::sin(zvalue);
    viskores::FloatDefault dy_dzs = rxy * std::cos(zsvalue) * std::cos(zvalue) - rxy * std::sin(zsvalue) * std::sin(zvalue);
    viskores::FloatDefault dy_dz = -rxy * std::sin(zsvalue) * std::sin(zvalue) + rxy * std::cos(zsvalue) * std::cos(zvalue);

    viskores::FloatDefault dz_dr = 0.0;
    viskores::FloatDefault dz_dzs = 0.0;
    viskores::FloatDefault dz_dz = 1.0;

    // Transform the tangent vector using the Jacobian
    viskores::FloatDefault vx3d = dx_dr * tangent[0] + dx_dzs * tangent[1] + dx_dz * tangent[2];
    viskores::FloatDefault vy3d = dy_dr * tangent[0] + dy_dzs * tangent[1] + dy_dz * tangent[2];
    viskores::FloatDefault vz3d = dz_dr * tangent[0] + dz_dzs * tangent[1] + dz_dz * tangent[2];
    tangentCartian = viskores::Vec3f(dx_dr, dy_dr, dz_dr);

    if (outputIt)
    {
      trajOut << "1, " << step << ", " << x3d << ", " << y3d << ", " << z3d << ", " << vx3d << ", " << vy3d << ", " << vz3d << ", " << yi << ", "
              << zsvalue << ", " << zvalue << std::endl;
    }

    return viskores::Vec3f(x3d, y3d, z3d);
  }

  template <typename FieldType1, typename FieldType2>
  viskores::FloatDefault LinearInterpolate(const FieldType1& x, const FieldType2& y, const viskores::FloatDefault val) const
  {
    // Perform binary search to find the interval.
    viskores::Id size = x.GetNumberOfValues();
    viskores::Id low = 0;
    viskores::Id high = size - 1;

    while (low <= high)
    {
      viskores::Id mid = low + (high - low) / 2;

      if (val >= x.Get(mid) && val <= x.Get(mid + 1))
      {
        auto x1 = x.Get(mid);
        auto x2 = x.Get(mid + 1);
        auto y1 = y.Get(mid);
        auto y2 = y.Get(mid + 1);

        auto result = y1 + (y2 - y1) * (val - x1) / (x2 - x1);
        return result;
      }
      else if (val < x.Get(mid))
        high = mid - 1;
      else
        low = mid + 1;
    }

    // Handle edge cases if the value is outside the range
    auto x1 = x.Get(0);
    auto x2 = x.Get(size - 1);
    auto y1 = y.Get(0);
    auto y2 = y.Get(size - 1);

    auto result = y1 + (y2 - y1) * (val - x1) / (x2 - x1);
    return result;
  }

  viskores::FloatDefault floatMod(viskores::FloatDefault val, viskores::FloatDefault mod_base) const
  {
    viskores::FloatDefault result = std::fmod(val, mod_base);
    if (result < 0)
      result += mod_base;
    return result;
  }

  template <typename BoutppType, typename CellSetType>
  VISKORES_EXEC void RK4Step(const viskores::Vec3f& pStart,
                             BoutppType& boutppField,
                             const CellSetType& cells,
                             viskores::Vec3f& pEnd,
                             viskores::Vec3f& tangent0,
                             viskores::Vec3f& tangent1) const
  {
    constexpr double twoPi = 2.0 * M_PI;
    constexpr double h = 1.0;
    constexpr double hh = h / 2.0;
    constexpr double h6 = h / 6.0;
    constexpr double direction = 1.0;
    constexpr viskores::Vec3f yPlusH(0, h, 0);

    //std::cout << "viskores::RK begin: " << pStart << std::endl;
    //Step 1
    auto res1 = this->Evaluate(pStart, boutppField, cells);
    tangent0 = { res1[0], 1, res1[1] };

    viskores::Vec3f p1;
    p1[0] = pStart[0] + direction * hh * res1[0];
    p1[1] = pStart[1];
    p1[2] = pStart[2] + direction * hh * res1[1];
    p1[2] = this->floatMod(p1[2], twoPi);
    //std::cout << "step1: " << res1[0] << " " << res1[1] << " :: " << p1[0] << " " << p1[2] << std::endl;

    //Step 2
    auto res2 = this->Evaluate(p1, boutppField, cells);
    res2 += this->Evaluate(p1 + yPlusH, boutppField, cells);
    res2[0] /= 2.0;
    res2[1] /= 2.0;
    viskores::Vec3f p2;
    p2[0] = pStart[0] + direction * hh * res2[0];
    p2[1] = pStart[1];
    p2[2] = pStart[2] + direction * hh * res2[1];
    p2[2] = this->floatMod(p2[2], twoPi);
    //std::cout << "step2: " << res2[0] << " " << res2[1] << " :: " << p2[0] << " " << p2[2] << std::endl;

    //Step 3
    auto res3 = this->Evaluate(p2, boutppField, cells);
    res3 += this->Evaluate(p2 + yPlusH, boutppField, cells);
    res3[0] /= 2.0;
    res3[1] /= 2.0;

    viskores::Vec3f p3;
    p3[0] = pStart[0] + direction * h * res3[0];
    p3[1] = pStart[1];
    p3[2] = pStart[2] + direction * h * res3[1];
    p3[2] = this->floatMod(p3[2], twoPi);
    //std::cout << "step3: " << res3[0] << " " << res3[1] << " :: " << p3[0] << " " << p3[2] << std::endl;

    //Step 4
    auto res4 = this->Evaluate(p3 + yPlusH, boutppField, cells);

    //viskores::Vec3f pEnd;
    pEnd[0] = pStart[0] + direction * h6 * (res1[0] + 2.0 * res2[0] + 2.0 * res3[0] + res4[0]);
    pEnd[1] = pStart[1] + h;
    pEnd[2] = pStart[2] + direction * h6 * (res1[1] + 2.0 * res2[1] + 2.0 * res3[1] + res4[1]);
    tangent1 = pEnd - pStart;
  }

  template <typename BoutppFieldType, typename CellSetType>
  VISKORES_EXEC bool Locate(const viskores::Vec3f& pt,
                            BoutppFieldType& boutppField,
                            const CellSetType& cells,
                            viskores::UInt8& cellShape,
                            viskores::IdComponent& numVerts,
                            viskores::VecVariable<viskores::Id, 8>& ptIndices,
                            viskores::Vec3f& parametric) const
  {
    viskores::Id cellId;
    boutppField.Locator3D.FindCell(pt, cellId, parametric);
    if (cellId == -1)
    {
      std::cout << "Locator failed: " << pt << std::endl;
      VISKORES_ASSERT(false);
      return false;
    }

    cellShape = cells.GetCellShape(cellId).Id;
    numVerts = cells.GetNumberOfIndices(cellId);
    ptIndices = cells.GetIndices(cellId);
    return true;
  }

  template <typename BoutppFieldType, typename ArrayType, typename CellSetType>
  VISKORES_EXEC viskores::FloatDefault Evaluate(const viskores::Vec3f& pt,
                                                BoutppFieldType& boutppField,
                                                const ArrayType& field,
                                                const CellSetType& cells) const
  {
    viskores::UInt8 cellShape;
    viskores::Vec3f parametric;
    viskores::IdComponent numVerts;
    viskores::VecVariable<viskores::Id, 8> ptIndices;

    this->Locate(pt, boutppField, cells, cellShape, numVerts, ptIndices, parametric);
    viskores::VecVariable<viskores::FloatDefault, 8> ptVals;
    for (viskores::IdComponent i = 0; i < numVerts; i++)
      ptVals.Append(field.Get(ptIndices[i]));

    viskores::FloatDefault val;
    auto status = viskores::exec::CellInterpolate(ptVals, parametric, cellShape, val);

    if (status != viskores::ErrorCode::Success)
      this->RaiseError("Evaluate failed.");
    return val;
  }

  template <typename BoutppFieldType, typename CellSetType>
  VISKORES_EXEC viskores::Vec2f Evaluate(const viskores::Vec3f& pt, BoutppFieldType& boutppField, const CellSetType& cells) const
  {
    viskores::UInt8 cellShape;
    viskores::Vec3f parametric;
    viskores::IdComponent numVerts;
    viskores::VecVariable<viskores::Id, 8> ptIndices;

    this->Locate(pt, boutppField, cells, cellShape, numVerts, ptIndices, parametric);

#if 0
      viskores::Id cellId;
      boutppField.Locator3D.FindCell(pt, cellId, parametric);
      if (cellId == -1)
      {
        std::cout << "Locator failed: " << pt << std::endl;
        VISKORES_ASSERT(false);
        this->RaiseError("Locator.FindCell failed.");
        }

        std::cout<<"cellId: "<<cellId<<" nvals= "<<boutppField.dxdy.GetNumberOfValues()<<std::endl;
        auto cellShape = cells.GetCellShape(cellId);
        viskores::IdComponent numVerts = cells.GetNumberOfIndices(cellId);
        viskores::VecVariable<viskores::Id, 8> ptIndices = cells.GetIndices(cellId);
#endif

    viskores::VecVariable<viskores::FloatDefault, 8> dxdyVals, dzdyVals;
    for (viskores::IdComponent i = 0; i < numVerts; i++)
    {
      dxdyVals.Append(boutppField.dxdy.Get(ptIndices[i]));
      dzdyVals.Append(boutppField.dzdy.Get(ptIndices[i]));
    }

    viskores::Vec2f result;
    auto status1 = viskores::exec::CellInterpolate(dxdyVals, parametric, cellShape, result[0]);
    auto status2 = viskores::exec::CellInterpolate(dzdyVals, parametric, cellShape, result[1]);
    if (status1 != viskores::ErrorCode::Success || status2 != viskores::ErrorCode::Success)
      this->RaiseError("CellInterpolate failed.");
    return result;
  }

  viskores::Id MaxPunctures;
  viskores::Id MaxSteps;
};