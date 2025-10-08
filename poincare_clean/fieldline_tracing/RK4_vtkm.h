#pragma once

#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/CellLocatorRectilinearGrid.h>
#include <viskores/cont/CellLocatorUniformGrid.h>
#include <viskores/cont/DataSet.h>
#include <viskores/exec/CellInterpolate.h>
#include <viskores/worklet/WorkletMapField.h>

#include "tricubicInterpolator.h"
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
//static std::ofstream fout, puncout;

static bool approxLT(viskores::FloatDefault a, viskores::FloatDefault b, viskores::FloatDefault tol = 1e-6)
{
  if (a < (b - tol))
    return true;
  return false;
}
static bool approxGT(viskores::FloatDefault a, viskores::FloatDefault b, viskores::FloatDefault tol = 1e-6)
{
  if (a > (b + tol))
    return true;
  return false;
}

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
  using LocatorType2D = viskores::exec::CellLocatorRectilinearGrid;
  using LocatorType3D = viskores::exec::CellLocatorRectilinearGrid;
  using ArrayType = viskores::cont::ArrayHandle<viskores::FloatDefault>;
  using ArrayPortalType = typename ArrayType::ReadPortalType;

public:
  BoutppFieldExecutionObject(const LocatorType3D& locator3D,
                             const LocatorType2D& locator2D,
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

  LocatorType2D Locator2D;
  LocatorType3D Locator3D;
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
  int ny;
  int nypf1, nypf2;
  int ixsep1, ixsep2;
  int divertor;
  int direction;

  using ControlSignature = void(FieldIn points,
                                FieldIn xind,
                                ExecObject boutppField,
                                WholeCellSetIn<> cells3D,
                                WholeCellSetIn<> cells2D,
                                WholeArrayOut validResult,
                                WholeArrayOut result,
                                WholeArrayOut resultStep,
                                WholeArrayOut puncIndex);
  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9);
  using InputDomain = _1;

  template <typename BoutppFieldType,
            typename CellSet3DType,
            typename CellSet2DType,
            typename ValidPuncType,
            typename ResultType,
            typename StepType,
            typename IndexType>
  VISKORES_EXEC void operator()(const viskores::Id& idx,
                                const viskores::Vec3f& pStart,
                                const viskores::FloatDefault& xindIn,
                                BoutppFieldType& boutppField,
                                const CellSet3DType& cells3D,
                                const CellSet2DType& cells2D,
                                ValidPuncType& validResult,
                                ResultType& result,
                                StepType& resultStep,
                                IndexType& puncIndex) const
  {
    constexpr double twoPi = 2.0 * M_PI;
    viskores::Id numPunc = 0;
    viskores::Id numSteps = 0;
    auto p0 = pStart;
    auto xind = xindIn;
    int region = this->GetRegion(xind, p0[1]);

    viskores::Id stepOffset = idx * this->MaxSteps;
    viskores::Id puncOffset = idx * this->MaxPunctures;
    //resultCart.Set(stepOffset, p0);
    std::ostringstream oss;
    oss << "steps." << std::setw(3) << std::setfill('0') << int(xind) << ".txt";
    auto stepOut = std::ofstream(oss.str());
    stepOut << "ID, STEP, X, Y, Z, REGION" << std::endl;
    stepOut << xind << ", " << 0 << ", " << p0[0] << ", " << p0[1] << ", " << p0[2] << ", " << region << std::endl;
    std::ostringstream oss2;
    oss2 << "punc." << std::setw(3) << std::setfill('0') << int(xind) << ".txt";
    std::cout << " ********** " << oss2.str() << std::endl;
    auto puncOut = std::ofstream(oss2.str());
    puncOut << "ID, STEP, X, Y, Z, Rxy, Zxy, Zvalue" << std::endl;

    auto rootPuncOut = std::ofstream("rootPunc.v.txt");
    rootPuncOut << "ID, STEP, X, Y, Z" << std::endl;

    viskores::Vec3f pCart0, pCart1, vCart, pIndex0, pIndex1, junk;

    //pCart0 = this->ConvertToCartesian(p0, boutppField, cells2D, 0, true);
    pCart0 = this->ToCartesian(p0, boutppField, cells2D);
    //vCart0 = this->EvaluateTangent(p0, boutppField, cells2D);

    pIndex0 = this->ConvertToIndex(p0, boutppField, 0);
    //resultIndex.Set(stepOffset, pIndex0);

    //    tanOut << "0, 0, " << pCart0[0] << ", " << pCart0[1] << ", " << pCart0[2] << ", " << vCart0[0] << ", " << vCart0[1] << ", " << vCart0[2]           << std::endl;
    //std::cout << "Begin: idx= " << idx << " pt= " << p0 << std::endl;

    int cnt = 0;
    for (viskores::Id step = 1; region < 10 && step < this->MaxSteps; step++)
    {
      bool inCFR =
        (region == 0 && p0[1] > static_cast<viskores::FloatDefault>(this->nypf1) && p0[1] < static_cast<viskores::FloatDefault>(this->nypf2));
      bool inSOLPF = (region == 1 || region == 2);

      viskores::Vec3f p1, tan0, tan1;
      this->RK4Step(p0, boutppField, cells3D, p1, tan0, tan1);

      bool inside = true;
      if (approxGT(p1[0], this->grid2DBounds.X.Max))
      {
        region = 12;
        inside = false;
        std::cout << "p1= " << p1 << " *** region 12" << std::endl;
        //VISKORES_ASSERT(false);
      }
      else if (approxLT(p1[0], this->grid2DBounds.X.Min))
      {
        region = 11;
        inside = false;
        std::cout << idx << ": p1= " << p1 << " *** region 11" << " step= " << step << std::endl;
      }

      //Always refresh xind from xEnd when still inside. Needed by twist-shift.
      if (inside)
        xind = this->LinearInterpolate(boutppField.XArray, boutppField.XiArray, p1[0]);

      if (inCFR)
      {
        if (xind > static_cast<viskores::FloatDefault>(this->ixsep1 - 1) + 0.5f)
          region = 1; // enters the SOL.
      }
      else
      {
        //SOL/PFR transitions among CFR/PFR.
        if (xind > static_cast<viskores::FloatDefault>(this->ixsep1 - 1) + 0.5f && p1[1] > this->nypf1 - 1 && p1[1] < this->nypf2)
          region = 0;
        else if ((xind < static_cast<viskores::FloatDefault>(this->ixsep1 - 1) + 0.5f) && (p1[1] > this->nypf2 - 1 - 1 || p1[1] < this->nypf1 - 1))
          region = 2;
      }

      //branch-specific fixes.
      if (inCFR && inside && region == 0)
      {
        //twist-shift at the branc cut (CFR only) -- now using refreshed xind
        if (direction == 1 && p0[1] == this->nypf2 - 1)
        {
          auto shiftangle = this->LinearInterpolate(boutppField.XiArray, boutppField.ShiftAngle, xind);
          p1[2] = p1[2] + shiftangle;
          p1[1] = this->nypf1;
        }
        else if (direction != 1)
        {
          VISKORES_ASSERT(false);
        }
      }
      else if (inSOLPF)
      {
        if (direction == 1 && p0[1] == this->nypf1 - 1 && region == 2)
          p1[1] = this->nypf2;
        else if (direction != 1 && p0[1] == this->nypf2 && region == 2)
        {
          p1[1] = this->nypf1 - 1;
          VISKORES_ASSERT(false);
        }
        if (direction == 1 && p1[1] == this->grid2DBounds.Y.Max - 1)
        {
          // particle reaches the outer diverter.
          region = 14;
        }
        else if (direction == -1 && p1[1] == 0)
        {
          // particle reaches the inner diverter.
          region = 13;
        }
      }

//do zshift stuff.
//Twist-shift at branch cut.
#if 0
      if (p0[1] == this->nypf2 - 1 && region == 0)
      {
        auto sa = this->LinearInterpolate(boutppField.XiArray, boutppField.ShiftAngle, xind);
        //double shiftAngle = 0; //scalarField1DEval(opts.XiArray, opts.ShiftAngle, xind);
        p1[2] = p1[2] + sa;
        p1[1] = this->nypf1;
      }
#endif

      p1[2] = this->floatMod(p1[2], twoPi);
      pIndex1 = this->ConvertToIndex(p0, boutppField, step);
      //resultIndex.Set(stepOffset + step, pIndex1);

      viskores::Vec3f ptXY(p1[0], p1[1], 0.0);
      viskores::FloatDefault rxyVal, zxyVal, zind, zvalue;
      bool res0 = this->EvaluateLinear(ptXY, boutppField, boutppField.rxy, cells2D, rxyVal);
      bool res1 = this->EvaluateLinear(ptXY, boutppField, boutppField.zxy, cells2D, zxyVal);
      if (!res0 || !res1)
        break;

      zind = this->LinearInterpolate(boutppField.ZArray, boutppField.ZiArray, p1[2]);
      zvalue = this->LinearInterpolate(boutppField.ZiArray, boutppField.ZArray, zind);


      stepOut << xind << ", " << step << ", " << p1[0] << ", " << p1[1] << ", " << p1[2] << ", (" << zind << "), " << region << ", " << rxyVal << ", "
              << zxyVal << ", " << zvalue << std::endl;
      //convert to cartesian.
      pCart0 = pCart1;
      viskores::Vec3f vCart;
      pCart1 = this->ToCartesian(p1, boutppField, cells2D);
      /*
      std::cout << step << ": " << xind << " " << p1[1] << " " << zind << " ---> " << rxyVal << " " << zxyVal << " " << zvalue << " " << pCart1
                << std::endl;
      */

      //trajOut << "0, " << step << " " << tmp[0] << ", " << tmp[1] << ", " << tmp[2] << ", 0, 0, 0, 0" << std::endl;

      //resultCart.Set(stepOffset + step, pCart1);
      //validSteps.Set(stepOffset + step, true);
      //pointId.Set(stepOffset + step, idx);

      //point crosses the X=0 plane.
      if (step > 1 && pCart0[0] * pCart1[0] < 0 && pCart0[0] > 0)
      {
        viskores::FloatDefault puncStep;
        auto puncPt = this->FindPuncture(step, p0, p1, boutppField, cells2D, puncStep);
        result.Set(puncOffset + numPunc, puncPt);
        resultStep.Set(puncOffset + numPunc, puncStep);
        rootPuncOut << xindIn << ", " << puncStep << ", " << puncPt[0] << ", " << puncPt[1] << ", " << puncPt[2] << std::endl;
        puncOut << xindIn << ", " << puncStep << ", " << puncPt[0] << ", " << puncPt[1] << ", " << puncPt[2] << ", 0.0, 0.0, 0.0" << std::endl;
        //std::cout << std::setprecision(5);
        //std::cout << "** " << step << " ZC: c= " << pCart0 << " " << pCart1 << "  :: " << pIndex0 << " " << pIndex1 << " p01= " << p0 << " " << p1
        //          << std::endl;
        //std::cout << "SET_PUNC: " << puncOffset << " " << numPunc << " " << (step - 1) << std::endl;
        puncIndex.Set(puncOffset + numPunc, xindIn);
        validResult.Set(puncOffset + numPunc, true);
        numPunc++;
      }

      //fout << step << ", " << xind << ", " << p1[0] << ", " << (int)p1[1] << ", " << p1[2] << ", " << region;

      p0 = p1;
      pIndex0 = pIndex1;
      if (numPunc == this->MaxPunctures)
      {
        std::cout << idx << " Finished" << std::endl;
        break;
      }
      cnt++;
    }

    std::cout << "Done. region= " << region << " #punc= " << numPunc << " #steps= " << cnt << std::endl;
  }

private:
  VISKORES_EXEC int GetRegion(const viskores::FloatDefault& xind, const viskores::FloatDefault& y) const
  {
    int region = 0;

    auto xxx = static_cast<viskores::FloatDefault>(this->ixsep1 - 1) + 0.5f;
    if (xind < xxx)
      std::cout << " Yes" << std::endl;
    if (xind < static_cast<viskores::FloatDefault>(this->ixsep1 - 1) + 0.5f)
    {
      region = 0; // closed flux surface.
      if (y < this->nypf1 + 1 || y > this->nypf2)
        region = 2; // Private flux region (PFR).
    }
    else
      region = 1; //Scrape off layer (SOL).

    if (this->divertor == 1) //single null
    {
      if (y == this->ny - 1 && this->direction == 1)
        region = 14; //starts on divertor.
      else if (y == 0 && direction == -1)
        region = 13; // starts on divertor.
    }


    return region;
  }


  //x/zind: [-0.135435384332,55,0.000486567073233] --> 149.790391597 0.019747086487
  template <typename BoutppFieldType, typename CellSetType>
  VISKORES_EXEC viskores::Vec3f FindPuncture(const viskores::Id& step,
                                             const viskores::Vec3f& p0,
                                             const viskores::Vec3f& p1,
                                             BoutppFieldType& boutppField,
                                             const CellSetType& cells,
                                             viskores::FloatDefault& stepPunc) const
  {
    // Evaluate endpoints in Cartesian
    viskores::Vec3f cart0 = this->ToCartesian(p0, boutppField, cells);
    viskores::Vec3f cart1 = this->ToCartesian(p1, boutppField, cells);
    viskores::FloatDefault f0 = cart0[0];
    viskores::FloatDefault f1 = cart1[0];

    const viskores::FloatDefault tol = static_cast<viskores::FloatDefault>(1e-6);

    // If either endpoint lies on the plane x=0, return it
    if (viskores::Abs(f0) <= tol)
      return cart0;
    if (viskores::Abs(f1) <= tol)
      return cart1;

    // If no sign change, return the endpoint closer to x=0
    if (f0 * f1 > 0)
      return (viskores::Abs(f0) < viskores::Abs(f1)) ? cart0 : cart1;

    // Bracketed: binary search on t in [0,1] where p(t) = p0 + t*(p1-p0)
    viskores::FloatDefault t_lo = 0, t_hi = 1;
    viskores::FloatDefault f_lo = f0, f_hi = f1;

    for (int iter = 0; iter < 100; ++iter)
    {
      viskores::FloatDefault t_mid = static_cast<viskores::FloatDefault>(0.5) * (t_lo + t_hi);

      viskores::Vec3f p_mid = p0 + t_mid * (p1 - p0);
      viskores::Vec3f cart_mid = this->ToCartesian(p_mid, boutppField, cells);
      viskores::FloatDefault f_mid = cart_mid[0];

      if (viskores::Abs(f_mid) <= tol || (t_hi - t_lo) <= static_cast<viskores::FloatDefault>(1e-9))
      {
        stepPunc = static_cast<viskores::FloatDefault>(step) + t_mid;
        return cart_mid; // Best estimate of puncture (Cartesian)
      }

      // Keep sub-interval that brackets the zero
      if ((f_lo < 0 && f_mid > 0) || (f_lo > 0 && f_mid < 0))
      {
        t_hi = t_mid;
        f_hi = f_mid;
      }
      else
      {
        t_lo = t_mid;
        f_lo = f_mid;
      }
    }

    // Max iterations reached: return midpoint estimate
    viskores::FloatDefault t_mid = static_cast<viskores::FloatDefault>(0.5) * (t_lo + t_hi);
    viskores::Vec3f p_mid = p0 + t_mid * (p1 - p0);
    stepPunc = static_cast<viskores::FloatDefault>(step) + t_mid;
    return this->ToCartesian(p_mid, boutppField, cells);
  }

  template <typename BoutppFieldType, typename CellSetType>
  VISKORES_EXEC viskores::Vec3f FindPuncture2(const viskores::Vec3f& p0,
                                              const viskores::Vec3f& p1,
                                              BoutppFieldType& boutppField,
                                              const CellSetType& cells) const
  {
    viskores::Vec3f puncPt;
    auto diffx = p1[0] - p0[0];
    auto diffy = p1[1] - p0[1];
    auto diffz = p1[2] - p0[2];

    int N = 100;
    viskores::FloatDefault _N = static_cast<viskores::FloatDefault>(N);
    viskores::Vec3f dXYZ(diffx / _N, diffy / _N, diffz / _N);
    auto pt = p0;
    for (int i = 0; i < N; i++)
    {
      auto cart = this->ToCartesian(pt, boutppField, cells);
      std::cout << "  iter_" << i << " cart= " << cart[0] << std::endl;
      pt = pt + dXYZ;
    }


    return puncPt;
  }

  template <typename BoutppFieldType, typename CellSetType>
  VISKORES_EXEC viskores::Vec3f EvaluateTangent(const viskores::Vec3f& pt0, BoutppFieldType& boutppField, const CellSetType& cells) const
  {
    constexpr double twoPi = 2.0 * M_PI;
    auto pt0Cart = this->ConvertToCartesian(pt0, boutppField, cells, -1, false);

    auto vec2d = this->EvaluateLinear(pt0, boutppField, cells);
    viskores::Vec3f vec0(vec2d[0], 1, vec2d[1]);
    auto pt1 = pt0 + vec0;
    pt1[2] = this->floatMod(pt1[2], twoPi);
    auto pt1Cart = this->ConvertToCartesian(pt1, boutppField, cells, -1, false);

    auto tangent = pt1Cart - pt0Cart;
    tangent = tangent * 0.75;
    return tangent;
  }

  template <typename BoutppFieldType>
  VISKORES_EXEC viskores::Vec3f ConvertToIndex(const viskores::Vec3f& pt, const BoutppFieldType& boutppField, int step) const
  {
    viskores::Vec3f ptIdx;

    ptIdx[0] = this->LinearInterpolate(boutppField.XArray, boutppField.XiArray, pt[0]);
    ptIdx[1] = static_cast<viskores::Id>(pt[1]);
    ptIdx[2] = this->LinearInterpolate(boutppField.ZArray, boutppField.ZiArray, pt[2]);
    //std::cout << "ConvertToIndex_" << step << ": " << pt << " --> " << ptIdx << std::endl;
    return ptIdx;
  }

  template <typename BoutppFieldType, typename CellSetType>
  VISKORES_EXEC viskores::Vec3f ToCartesian(const viskores::Vec3f& pt, BoutppFieldType& boutppField, const CellSetType& cells) const
  {
    viskores::Vec3f ptXY(pt[0], pt[1], 0.0);
    auto zind = this->LinearInterpolate(boutppField.ZArray, boutppField.ZiArray, pt[2]);

    auto rxy = this->EvaluateCubic1D(ptXY, "rxy");
    auto zxy = this->EvaluateCubic1D(ptXY, "zxy");
    auto zvalue = this->LinearInterpolate(boutppField.ZiArray, boutppField.ZArray, zind);
    //auto zsvalue = this->Evaluate(ptXY, boutppField, boutppField.zShift, cells);
    auto zsvalue = this->EvaluateCubic1D(ptXY, "zShift");

    auto x3d = rxy * std::cos(zsvalue) * std::cos(zvalue) - rxy * std::sin(zsvalue) * std::sin(zvalue);
    auto y3d = rxy * std::cos(zsvalue) * std::sin(zvalue) + rxy * std::sin(zsvalue) * std::cos(zvalue);
    auto z3d = zxy;

    return viskores::Vec3f(x3d, y3d, z3d);
  }

  template <typename BoutppFieldType, typename CellSetType>
  VISKORES_EXEC viskores::Vec3f ConvertToCartesian(const viskores::FloatDefault& _xind,
                                                   const viskores::Vec3f& pt,
                                                   const viskores::Vec3f& vec,
                                                   BoutppFieldType& boutppField,
                                                   const CellSetType& cells,
                                                   int step,
                                                   bool outputIt,
                                                   viskores::Vec3f& ptIndex,
                                                   std::ofstream& puncOut) const
  {
    viskores::Id yi = static_cast<viskores::Id>(pt[1]);
    auto xind = this->LinearInterpolate(boutppField.XArray, boutppField.XiArray, pt[0]);
    auto zind = this->LinearInterpolate(boutppField.ZArray, boutppField.ZiArray, pt[2]);
    ptIndex[0] = xind;
    ptIndex[1] = yi;
    ptIndex[2] = zind;

    viskores::Vec3f ptXY(pt[0], pt[1], 0.0);
    //auto rxy = this->Evaluate(ptXY, boutppField, boutppField.rxy, cells);
    //auto zxy = this->Evaluate(ptXY, boutppField, boutppField.zxy, cells);
    auto rxy = this->EvaluateCubic1D(ptXY, "rxy");
    auto zxy = this->EvaluateCubic1D(ptXY, "zxy");
    auto zvalue = this->LinearInterpolate(boutppField.ZiArray, boutppField.ZArray, zind);
    //auto zsvalue = this->Evaluate(ptXY, boutppField, boutppField.zShift, cells);
    auto zsvalue = this->EvaluateCubic1D(ptXY, "zShift");

    //std::cout << "***** Rxy(" << xind << ", " << ptXY[1] << ")= " << rxy << std::endl;
    //std::cout << "***** Zxy(" << xind << ", " << ptXY[1] << ")= " << zxy << std::endl;
    //std::cout << "***** Zs= " << zsvalue << std::endl;
    //std::cout << "***** Zvalue= " << zvalue << std::endl;
    // Compute x3d and y3d
    viskores::FloatDefault x3d = rxy * std::cos(zsvalue) * std::cos(zvalue) - rxy * std::sin(zsvalue) * std::sin(zvalue);
    viskores::FloatDefault y3d = rxy * std::cos(zsvalue) * std::sin(zvalue) + rxy * std::sin(zsvalue) * std::cos(zvalue);
    viskores::FloatDefault z3d = zxy;
    puncOut << xind << ", " << step << ", " << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << rxy << ", " << zxy << ", " << zvalue << std::endl;

    //std::cout << "***** XYZ: " << x3d << " " << y3d << " " << z3d << std::endl;
#ifndef VISKORES_HIP
    //if (outputIt)
    //  trajOut << "0, " << step << ", " << x3d << ", " << y3d << ", " << z3d << ", 0, " << yi << ", " << zsvalue << ", " << zvalue << std::endl;
#endif

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

      /*
      vecCart[0] = vr * cos_zs * cos_z - vr * sin_zs * sin_z - rxy * vtheta * sin_zs_z;
      vecCart[1] = vr * cos_zs * sin_z + vr * sin_zs * cos_z + rxy * vtheta * cos_zs_z;
      vecCart[2] = vz; // Z component is directly mapped
      */
      //vecCart = vecCart * 0.50;
      //std::cout << " *********** vecCart= " << vecCart << std::endl;
    }

    return viskores::Vec3f(x3d, y3d, z3d);
  }

  template <typename BoutppFieldType, typename CellSetType>
  VISKORES_EXEC viskores::Vec3f ConvertToCartesianWithTangent(const viskores::Vec3f& pt,
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
    auto rxy = this->EvaluateLinear(ptXY, boutppField, boutppField.rxy, cells);
    auto zxy = this->EvaluateLinear(ptXY, boutppField, boutppField.zxy, cells);
    auto zvalue = this->LinearInterpolate(boutppField.ZiArray, boutppField.ZArray, zind);
    auto zsvalue = this->EvaluateLinear(ptXY, boutppField, boutppField.zShift, cells);

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
      //trajOut << "1, " << step << ", " << x3d << ", " << y3d << ", " << z3d << ", " << vx3d << ", " << vy3d << ", " << vz3d << ", " << yi << ", "
      //        << zsvalue << ", " << zvalue << std::endl;
    }

    return viskores::Vec3f(x3d, y3d, z3d);
  }

  template <typename FieldType1, typename FieldType2>
  VISKORES_EXEC viskores::FloatDefault LinearInterpolate(const FieldType1& x, const FieldType2& y, const viskores::FloatDefault val) const
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

  VISKORES_EXEC viskores::FloatDefault floatMod(viskores::FloatDefault val, viskores::FloatDefault mod_base) const
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
    //constexpr double h = 1.0;
    constexpr double h = 1.0;
    constexpr double hh = h / 2.0;
    constexpr double h6 = h / 6.0;
    const viskores::FloatDefault dirVal = static_cast<viskores::FloatDefault>(this->direction);
    constexpr viskores::Vec3f yPlusH(0, h, 0);

    //std::cout << "viskores::RK begin: " << pStart << std::endl;
    //Step 1
    //auto res1 = this->Evaluate(pStart, boutppField, cells);
    auto res1 = this->EvaluateCubic2D(pStart);
    tangent0 = { res1[0], 1, res1[1] };

    viskores::Vec3f p1;
    p1[0] = pStart[0] + dirVal * hh * res1[0];
    p1[1] = pStart[1];
    p1[2] = pStart[2] + dirVal * hh * res1[1];
    p1[2] = this->floatMod(p1[2], twoPi);
    //std::cout << "step1: " << res1[0] << " " << res1[1] << " :: " << p1[0] << " " << p1[2] << std::endl;

    //Step 2
    //auto res2 = this->EvaluateLinear(p1, boutppField, cells);
    //res2 += this->EvaluateLinear(p1 + yPlusH, boutppField, cells);
    auto res2 = this->EvaluateCubic2D(p1);
    res2 += this->EvaluateCubic2D(p1 + yPlusH);

    res2[0] /= 2.0;
    res2[1] /= 2.0;
    viskores::Vec3f p2;
    p2[0] = pStart[0] + dirVal * hh * res2[0];
    p2[1] = pStart[1];
    p2[2] = pStart[2] + dirVal * hh * res2[1];
    p2[2] = this->floatMod(p2[2], twoPi);
    //std::cout << "step2: " << res2[0] << " " << res2[1] << " :: " << p2[0] << " " << p2[2] << std::endl;

    //Step 3
    //auto res3 = this->Evaluate(p2, boutppField, cells);
    //res3 += this->Evaluate(p2 + yPlusH, boutppField, cells);

    auto res3 = this->EvaluateCubic2D(p2);
    res3 += this->EvaluateCubic2D(p2 + yPlusH);
    res3[0] /= 2.0;
    res3[1] /= 2.0;

    viskores::Vec3f p3;
    p3[0] = pStart[0] + dirVal * h * res3[0];
    p3[1] = pStart[1];
    p3[2] = pStart[2] + dirVal * h * res3[1];
    p3[2] = this->floatMod(p3[2], twoPi);
    //std::cout << "step3: " << res3[0] << " " << res3[1] << " :: " << p3[0] << " " << p3[2] << std::endl;

    //Step 4
    //auto res4 = this->Evaluate(p3 + yPlusH, boutppField, cells);
    auto res4 = this->EvaluateCubic2D(p3 + yPlusH);

    //viskores::Vec3f pEnd;
    pEnd[0] = pStart[0] + dirVal * h6 * (res1[0] + 2.0 * res2[0] + 2.0 * res3[0] + res4[0]);
    pEnd[1] = pStart[1] + h;
    pEnd[2] = pStart[2] + dirVal * h6 * (res1[1] + 2.0 * res2[1] + 2.0 * res3[1] + res4[1]);
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
      auto pt1 = pt - viskores::Vec3f(0, 1e-5, 0);
      auto pt2 = pt + viskores::Vec3f(0, 1e-5, 0);
      if (boutppField.Locator3D.IsInside(pt1))
        boutppField.Locator3D.FindCell(pt1, cellId, parametric);
      else if (boutppField.Locator3D.IsInside(pt2))
        boutppField.Locator3D.FindCell(pt2, cellId, parametric);
      else
      {
        return false;
      }
    }

    cellShape = cells.GetCellShape(cellId).Id;
    numVerts = cells.GetNumberOfIndices(cellId);
    ptIndices = cells.GetIndices(cellId);
    return true;
  }

  VISKORES_EXEC viskores::FloatDefault EvaluateCubic1D(const viskores::Vec3f& pt, const std::string& fieldName) const
  {
    viskores::FloatDefault res;
    res = CubicEval(this->ds2D, fieldName, pt);
    return res;
  }

  VISKORES_EXEC viskores::Vec2f EvaluateCubic2D(const viskores::Vec3f& pt) const
  {
    std::string dxdy = "dxdy", dzdy = "dzdy";

    viskores::Vec2f res;
    res[0] = CubicEval(this->ds3D, dxdy, pt);
    res[1] = CubicEval(this->ds3D, dzdy, pt);
    return res;
  }

  template <typename BoutppFieldType, typename ArrayType, typename CellSetType>
  VISKORES_EXEC viskores::FloatDefault EvaluateLinear(const viskores::Vec3f& pt,
                                                      BoutppFieldType& boutppField,
                                                      const ArrayType& field,
                                                      const CellSetType& cells,
                                                      viskores::FloatDefault& val) const
  {
    viskores::UInt8 cellShape;
    viskores::Vec3f parametric;
    viskores::IdComponent numVerts;
    viskores::VecVariable<viskores::Id, 8> ptIndices;

    bool res = this->Locate(pt, boutppField, cells, cellShape, numVerts, ptIndices, parametric);
    if (!res)
      return res;
    viskores::VecVariable<viskores::FloatDefault, 8> ptVals;
    for (viskores::IdComponent i = 0; i < numVerts; i++)
      ptVals.Append(field.Get(ptIndices[i]));

    auto status = viskores::exec::CellInterpolate(ptVals, parametric, cellShape, val);

    if (status != viskores::ErrorCode::Success)
      this->RaiseError("Evaluate failed.");
    return true;
  }

  template <typename BoutppFieldType, typename CellSetType>
  VISKORES_EXEC viskores::Vec2f EvaluateLinear(const viskores::Vec3f& pt, BoutppFieldType& boutppField, const CellSetType& cells) const
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
      std::cout << "dxdy: " << boutppField.dxdy.GetNumberOfValues() << std::endl;
      std::cout << "  dx: " << ptIndices[i] << std::endl;
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

public:
  viskores::cont::DataSet ds3D, ds2D;
};
