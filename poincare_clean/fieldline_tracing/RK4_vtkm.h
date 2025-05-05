#pragma once

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CellLocatorRectilinearGrid.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/exec/CellInterpolate.h>
#include <vtkm/worklet/WorkletMapField.h>

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

void writeArray1DToFile(std::vector<double>& array, const std::string& fname);

void writeArray2DToFile(std::vector<std::vector<double>>& array, const std::string& fname);

void writeArray3DToFile(const std::vector<std::vector<std::vector<double>>>& array, const std::string& fname);

vtkm::Vec3f RK4_FLT1_vtkm(const vtkm::Vec3f& pStart,
                          const vtkm::cont::DataSet& grid2D,
                          const vtkm::cont::DataSet& grid2D_cfr,
                          const vtkm::cont::DataSet& grid2D_xz,
                          const vtkm::cont::DataSet& grid3D,
                          const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& xarray,
                          const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& zarray,
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
  using LocatorType = vtkm::exec::CellLocatorRectilinearGrid;
  using ArrayType = vtkm::cont::ArrayHandle<vtkm::FloatDefault>;
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
    , ZiArray(ziarray)
    , ZArray(zarray)
    , ShiftAngle(shiftAngle)
  {
  }

  LocatorType Locator2D, Locator3D;
  ArrayPortalType dxdy, dzdy;
  ArrayPortalType rxy, zxy, zShift;
  ArrayPortalType XiArray, XArray, ZiArray, ZArray, ShiftAngle;
};

class BoutppField : public vtkm::cont::ExecutionObjectBase
{
  using ArrayType = vtkm::cont::ArrayHandle<vtkm::FloatDefault>;
  using Structured2DType = vtkm::cont::CellSetStructured<2>;
  using Structured3DType = vtkm::cont::CellSetStructured<3>;
  //using Structured3DType = vtkm::exec::CellSetStructured<3>;
  //using Structured2DType = vtkm::exec::CellSetStructured<2>;
  //using CellSet2DExecType = vtkm::exec::CellSet<Structured2DType>;
  //using CellSet3DExecType = vtkm::exec::CellSet<Structured3DType>;
  //using X2Type = vtkm::exec::ConnectivityStructured<vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}, 2>;
  //using X3Type = vtkm::exec::ConnectivityStructured<vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}, 3>;


public:
  BoutppField(const vtkm::cont::DataSet& grid3D,
              const vtkm::cont::DataSet& grid2D,
              const ArrayType& xiarray,
              const ArrayType& xarray,
              const ArrayType& ziarray,
              const ArrayType& zarray,
              const ArrayType& shiftAngle)
    : Grid3D(grid3D)
    , Grid2D(grid2D)
    , XiArray(xiarray)
    , XArray(xarray)
    , ZiArray(ziarray)
    , ZArray(zarray)
    , ShiftAngle(shiftAngle)
  {
  }

  VTKM_CONT BoutppFieldExecutionObject PrepareForExecution(vtkm::cont::DeviceAdapterId device, vtkm::cont::Token& token) const
  {
    vtkm::cont::CellLocatorRectilinearGrid locator3D;
    vtkm::cont::CellLocatorRectilinearGrid locator2D;
    locator3D.SetCoordinates(this->Grid3D.GetCoordinateSystem());
    locator3D.SetCellSet(this->Grid3D.GetCellSet());
    locator2D.SetCoordinates(this->Grid2D.GetCoordinateSystem());
    locator2D.SetCellSet(this->Grid2D.GetCellSet());
    locator3D.Update();
    locator2D.Update();

    //get the fields.
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> dxdyField, dzdyField, rxyField, zxyField, zShiftField;
    this->Grid3D.GetField("dxdy").GetData().AsArrayHandle<vtkm::FloatDefault>(dxdyField);
    this->Grid3D.GetField("dzdy").GetData().AsArrayHandle<vtkm::FloatDefault>(dzdyField);
    this->Grid2D.GetField("rxy").GetData().AsArrayHandle<vtkm::FloatDefault>(rxyField);
    this->Grid2D.GetField("zxy").GetData().AsArrayHandle<vtkm::FloatDefault>(zxyField);
    this->Grid2D.GetField("zShift").GetData().AsArrayHandle<vtkm::FloatDefault>(zShiftField);

    return BoutppFieldExecutionObject(locator3D.PrepareForExecution(device, token),
                                      locator2D.PrepareForExecution(device, token),
                                      dxdyField.PrepareForInput(device, token),
                                      dzdyField.PrepareForInput(device, token),
                                      rxyField.PrepareForInput(device, token),
                                      zxyField.PrepareForInput(device, token),
                                      zShiftField.PrepareForInput(device, token),
                                      this->XiArray.PrepareForInput(device, token),
                                      this->XArray.PrepareForInput(device, token),
                                      this->ZiArray.PrepareForInput(device, token),
                                      this->ZArray.PrepareForInput(device, token),
                                      this->ShiftAngle.PrepareForInput(device, token));
  }

  vtkm::cont::DataSet Grid3D;
  vtkm::cont::DataSet Grid2D;
  ArrayType XiArray, XArray;
  ArrayType ZiArray, ZArray;
  ArrayType ShiftAngle;
};

extern std::ofstream stepOut, trajOut;
class RK4Worklet : public vtkm::worklet::WorkletMapField
{
public:
  RK4Worklet(vtkm::Id maxPuncs, vtkm::Id maxSteps)
    : MaxPunctures(maxPuncs)
    , MaxSteps(maxSteps)
  {
  }

  vtkm::FloatDefault xMin, xMax;
  vtkm::Bounds grid3DBounds, grid2DBounds;
  int nypf1, nypf2;
  int ixsep1, ixsep2;

  using ControlSignature =
    void(FieldIn points, ExecObject boutppField, WholeCellSetIn<> cells3D, WholeCellSetIn<> cells2D, WholeArrayOut puncIndices, WholeArrayOut result);
  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6);
  using InputDomain = _1;

  template <typename BoutppFieldType, typename CellSet3DType, typename CellSet2DType, typename IndexType, typename ResultType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            const vtkm::Vec3f& pStart,
                            BoutppFieldType& boutppField,
                            const CellSet3DType& cells3D,
                            const CellSet2DType& cells2D,
                            IndexType& puncIndex,
                            ResultType& result) const
  {
    constexpr double twoPi = 2.0 * M_PI;
    vtkm::Id numPunc = 0;
    vtkm::Id numSteps = 0;
    auto p0 = pStart;

    vtkm::Id stepOffset = idx * this->MaxSteps;
    vtkm::Id puncOffset = idx * this->MaxPunctures;
    result.Set(stepOffset, p0);
    //stepOut << idx << ", " << 0 << ", " << p0[0] << ", " << p0[1] << ", " << p0[2] << std::endl;
    int region = 0;

    auto tmp = this->ConvertToCartesian(p0, boutppField, cells2D, 0);

    for (vtkm::Id step = 1; region == 0 && step < this->MaxSteps; step++)
    {
      if (step == 58)
        std::cout << " ***** Issue." << std::endl;
      auto p1 = this->RK4Step(p0, boutppField, cells3D);
      auto xind = this->LinearInterpolate(boutppField.XArray, boutppField.XiArray, p1[0]);
      if (p1[0] > this->grid2DBounds.X.Max)
        region = 12;
      else if (p1[0] < this->grid2DBounds.X.Min)
        region = 11;
      else
      {
        if (xind > static_cast<vtkm::FloatDefault>(this->ixsep1) + 0.5)
          region = 1;
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
      auto tmp = this->ConvertToCartesian(p1, boutppField, cells2D, step);
      //trajOut << "0, " << step << " " << tmp[0] << ", " << tmp[1] << ", " << tmp[2] << ", 0, 0, 0, 0" << std::endl;

      result.Set(stepOffset + step, tmp);

      //point crosses the X=0 plane.
      if (p0[0] * p1[0] < 0)
      {
        puncIndex.Set(puncOffset + numPunc, step - 1);
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
  vtkm::Vec3f ConvertToCartesian(const vtkm::Vec3f& pt, BoutppFieldType& boutppField, const CellSetType& cells, int step) const
  {
    vtkm::Vec3f pCart;
    vtkm::Id yi = static_cast<vtkm::Id>(pt[1]);
    auto xind = this->LinearInterpolate(boutppField.XArray, boutppField.XiArray, pt[0]);
    auto traj2 = xind;
    auto zind = this->LinearInterpolate(boutppField.ZArray, boutppField.ZiArray, pt[2]);
    auto traj4 = zind;
    //vtkm::Vec3f ptXZ(pt[0], pt[2], 0.0);
    vtkm::Vec3f ptXY(pt[0], pt[1], 0.0);
    auto rxyvalue = this->Evaluate(ptXY, boutppField, boutppField.rxy, cells);
    auto zxyvalue = this->Evaluate(ptXY, boutppField, boutppField.zxy, cells);

    auto zsvalue = this->Evaluate(ptXY, boutppField, boutppField.zShift, cells);
    //double zvalue = INTERP(opts.ziarray, opts.zarray, pt.traj4); traj4 = zind
    auto zvalue = this->LinearInterpolate(boutppField.ZiArray, boutppField.ZArray, zind);


    auto rxy = rxyvalue;
    vtkm::Vec3f tmp(pt[0], yi, pt[2]);

    vtkm::FloatDefault x3d_tmp = rxyvalue * std::cos(zsvalue);
    vtkm::FloatDefault y3d_tmp = rxyvalue * std::sin(zsvalue);

    // Compute x3d and y3d
    vtkm::FloatDefault x3d = x3d_tmp * std::cos(zvalue) - y3d_tmp * std::sin(zvalue);
    vtkm::FloatDefault y3d = x3d_tmp * std::sin(zvalue) + y3d_tmp * std::cos(zvalue);
    vtkm::FloatDefault z3d = zxyvalue;
    trajOut << "0, " << step << ", " << x3d << ", " << y3d << ", " << z3d << ", 0, " << yi << ", " << zsvalue << ", " << zvalue << std::endl;

    return vtkm::Vec3f(x3d, y3d, z3d);
  }

  template <typename FieldType1, typename FieldType2>
  vtkm::FloatDefault LinearInterpolate(const FieldType1& x, const FieldType2& y, const vtkm::FloatDefault val) const
  {
    // Find the interval.
    vtkm::Id idx = -1;
    vtkm::Id size = x.GetNumberOfValues();
    for (vtkm::Id i = 0; i < size - 1; i++)
    {
      if (val >= x.Get(i) && val <= x.Get(i + 1))
      {
        idx = i;
        break;
      }
    }
    // Perform linear interpolation
    auto x1 = x.Get(idx);
    auto x2 = x.Get(idx + 1);
    auto y1 = y.Get(idx);
    auto y2 = y.Get(idx + 1);

    auto result = y1 + (y2 - y1) * (val - x1) / (x2 - x1);
    return result;
  }

  vtkm::FloatDefault floatMod(vtkm::FloatDefault val, vtkm::FloatDefault mod_base) const
  {
    vtkm::FloatDefault result = std::fmod(val, mod_base);
    if (result < 0)
      result += mod_base;
    return result;
  }

  template <typename BoutppType, typename CellSetType>
  VTKM_EXEC vtkm::Vec3f RK4Step(const vtkm::Vec3f& pStart, BoutppType& boutppField, const CellSetType& cells) const
  {
    constexpr vtkm::Vec3f yPlus1(0, 1, 0);
    constexpr double twoPi = 2.0 * M_PI;
    constexpr double h = 1.0;
    constexpr double hh = h / 2.0;
    constexpr double h6 = h / 6.0;
    constexpr double direction = 1.0;

    std::cout << "vtkm::RK begin: " << pStart << std::endl;
    //Step 1
    auto res1 = this->Evaluate(pStart, boutppField, cells);
    vtkm::Vec3f p1;
    p1[0] = pStart[0] + direction * hh * res1[0];
    p1[1] = pStart[1];
    p1[2] = pStart[2] + direction * hh * res1[1];
    p1[2] = this->floatMod(p1[2], twoPi);
    std::cout << "step1: " << res1[0] << " " << res1[1] << " :: " << p1[0] << " " << p1[2] << std::endl;

    //Step 2
    auto res2 = this->Evaluate(p1, boutppField, cells);
    res2 += this->Evaluate(p1 + yPlus1, boutppField, cells);
    res2[0] /= 2.0;
    res2[1] /= 2.0;
    vtkm::Vec3f p2;
    p2[0] = pStart[0] + direction * hh * res2[0];
    p2[1] = pStart[1];
    p2[2] = pStart[2] + direction * hh * res2[1];
    p2[2] = this->floatMod(p2[2], twoPi);
    std::cout << "step2: " << res2[0] << " " << res2[1] << " :: " << p2[0] << " " << p2[2] << std::endl;

    //Step 3
    auto res3 = this->Evaluate(p2, boutppField, cells);
    res3 += this->Evaluate(p2 + yPlus1, boutppField, cells);
    res3[0] /= 2.0;
    res3[1] /= 2.0;

    vtkm::Vec3f p3;
    p3[0] = pStart[0] + direction * h * res3[0];
    p3[1] = pStart[1];
    p3[2] = pStart[2] + direction * h * res3[1];
    p3[2] = this->floatMod(p3[2], twoPi);
    std::cout << "step3: " << res3[0] << " " << res3[1] << " :: " << p3[0] << " " << p3[2] << std::endl;

    //Step 4
    auto res4 = this->Evaluate(p3 + yPlus1, boutppField, cells);

    vtkm::Vec3f pEnd;
    pEnd[0] = pStart[0] + direction * h6 * (res1[0] + 2.0 * res2[0] + 2.0 * res3[0] + res4[0]);
    pEnd[1] = pStart[1] + 1;
    pEnd[2] = pStart[2] + direction * h6 * (res1[1] + 2.0 * res2[1] + 2.0 * res3[1] + res4[1]);
    return pEnd;
  }

  template <typename BoutppFieldType, typename CellSetType>
  VTKM_EXEC bool Locate(const vtkm::Vec3f& pt,
                        BoutppFieldType& boutppField,
                        const CellSetType& cells,
                        vtkm::UInt8& cellShape,
                        vtkm::IdComponent& numVerts,
                        vtkm::VecVariable<vtkm::Id, 8>& ptIndices,
                        vtkm::Vec3f& parametric) const
  {
    vtkm::Id cellId;
    boutppField.Locator3D.FindCell(pt, cellId, parametric);
    if (cellId == -1)
    {
      VTKM_ASSERT(false);
      std::cout << "Locator failed: " << pt << std::endl;
      return false;
    }

    cellShape = cells.GetCellShape(cellId).Id;
    numVerts = cells.GetNumberOfIndices(cellId);
    ptIndices = cells.GetIndices(cellId);
    return true;
  }

  template <typename BoutppFieldType, typename ArrayType, typename CellSetType>
  VTKM_EXEC vtkm::FloatDefault Evaluate(const vtkm::Vec3f& pt, BoutppFieldType& boutppField, const ArrayType& field, const CellSetType& cells) const
  {
    vtkm::UInt8 cellShape;
    vtkm::Vec3f parametric;
    vtkm::IdComponent numVerts;
    vtkm::VecVariable<vtkm::Id, 8> ptIndices;

    this->Locate(pt, boutppField, cells, cellShape, numVerts, ptIndices, parametric);
    vtkm::VecVariable<vtkm::FloatDefault, 8> ptVals;
    for (vtkm::IdComponent i = 0; i < numVerts; i++)
      ptVals.Append(field.Get(ptIndices[i]));

    vtkm::FloatDefault val;
    auto status = vtkm::exec::CellInterpolate(ptVals, parametric, cellShape, val);

    if (status != vtkm::ErrorCode::Success)
      this->RaiseError("Evaluate failed.");
    return val;
  }

  template <typename BoutppFieldType, typename CellSetType>
  VTKM_EXEC vtkm::Vec2f Evaluate(const vtkm::Vec3f& pt, BoutppFieldType& boutppField, const CellSetType& cells) const
  {
    vtkm::UInt8 cellShape;
    vtkm::Vec3f parametric;
    vtkm::IdComponent numVerts;
    vtkm::VecVariable<vtkm::Id, 8> ptIndices;

    this->Locate(pt, boutppField, cells, cellShape, numVerts, ptIndices, parametric);

#if 0
      vtkm::Id cellId;
      boutppField.Locator3D.FindCell(pt, cellId, parametric);
      if (cellId == -1)
      {
        std::cout << "Locator failed: " << pt << std::endl;
        VTKM_ASSERT(false);
        this->RaiseError("Locator.FindCell failed.");
        }

        std::cout<<"cellId: "<<cellId<<" nvals= "<<boutppField.dxdy.GetNumberOfValues()<<std::endl;
        auto cellShape = cells.GetCellShape(cellId);
        vtkm::IdComponent numVerts = cells.GetNumberOfIndices(cellId);
        vtkm::VecVariable<vtkm::Id, 8> ptIndices = cells.GetIndices(cellId);
#endif

    vtkm::VecVariable<vtkm::FloatDefault, 8> dxdyVals, dzdyVals;
    for (vtkm::IdComponent i = 0; i < numVerts; i++)
    {
      dxdyVals.Append(boutppField.dxdy.Get(ptIndices[i]));
      dzdyVals.Append(boutppField.dzdy.Get(ptIndices[i]));
    }

    vtkm::Vec2f result;
    auto status1 = vtkm::exec::CellInterpolate(dxdyVals, parametric, cellShape, result[0]);
    auto status2 = vtkm::exec::CellInterpolate(dzdyVals, parametric, cellShape, result[1]);
    if (status1 != vtkm::ErrorCode::Success || status2 != vtkm::ErrorCode::Success)
      this->RaiseError("CellInterpolate failed.");
    return result;
  }

  vtkm::Id MaxPunctures;
  vtkm::Id MaxSteps;
};