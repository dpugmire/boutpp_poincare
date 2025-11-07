#pragma once

#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/CellLocatorRectilinearGrid.h>
#include <viskores/cont/CellLocatorUniformGrid.h>
#include <viskores/cont/DataSet.h>
#include <viskores/exec/CellInterpolate.h>
#include <viskores/worklet/WorkletMapField.h>

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
//static std::ofstream fout, puncout;

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

  using ControlSignature = void(FieldIn IDs,
                                FieldIn points,
                                ExecObject boutppField,
                                WholeCellSetIn<> cells3D,
                                WholeCellSetIn<> cells2D,
                                WholeArrayOut puncIDs,
                                WholeArrayOut puncIndex,
                                WholeArrayOut resultRZThetaPsi);
  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7, _8);
  using InputDomain = _1;

  template <typename BoutppFieldType, typename CellSet3DType, typename CellSet2DType, typename PuncIndexType, typename ResultType>
  VISKORES_EXEC void operator()(const viskores::Id& idx,
                                const viskores::Id& ID,
                                const viskores::Vec3f& pStart,
                                BoutppFieldType& boutppField,
                                const CellSet3DType& cells3D,
                                const CellSet2DType& cells2D,
                                PuncIndexType& puncID,
                                PuncIndexType& puncIndex,
                                ResultType& resultRZThetaPsi) const
  {
    constexpr double twoPi = 2.0 * M_PI;
    viskores::Id numPunc = 0;
    viskores::Id numSteps = 0;
    auto p0 = pStart;
    auto xind = p0[0];
    int region = this->GetRegion(xind, p0[1]);

    viskores::Id stepOffset = idx * this->MaxSteps;
    viskores::Id puncOffset = idx * this->MaxPunctures;

    puncIndex.Set(puncOffset, numPunc);
#if 1
#if TRACE_IO
    std::ostringstream oss;
    oss << "steps." << std::setw(3) << std::setfill('0') << int(idx) << ".txt";
    auto stepOut = std::ofstream(oss.str());
    stepOut << "Xind0, STEP, X, Y, Z, REGION" << std::endl;
    stepOut << xind << ", " << 0 << ", " << p0[0] << ", " << p0[1] << ", " << p0[2] << ", " << region << std::endl;
    std::ostringstream oss2;
    oss2 << "punc." << std::setw(3) << std::setfill('0') << int(idx) << ".txt";
    std::cout << " ********** " << oss2.str() << std::endl;
    auto puncOut = std::ofstream(oss2.str());
    puncOut << "Xind0, STEP, X, Y, Z, Rxy, Zxy, Zvalue" << std::endl;

    auto rootPuncOut = std::ofstream("rootPunc.v.txt");
    rootPuncOut << "ID, STEP, X, Y, Z" << std::endl;
#endif

    viskores::Vec3f pCart0, pCart1, vCart, pIndex0, pIndex1;

    pCart0 = this->ToCartesian(p0, boutppField, cells2D);
    pIndex0 = this->ConvertToIndex(p0, boutppField, 0);

    int cnt = 0;
    for (viskores::Id step = 1; region < 10 && step < this->MaxSteps; step++)
    {
      bool inCFR =
        (region == 0 && p0[1] > static_cast<viskores::FloatDefault>(this->nypf1) && p0[1] < static_cast<viskores::FloatDefault>(this->nypf2));
      bool inSOLPF = (region == 1 || region == 2);

      viskores::Vec3f p1, tan0, tan1;
      this->RK4Step(p0, boutppField, cells3D, p1, tan0, tan1);

      bool inside = true;
      if (this->ApproxGT(p1[0], this->grid2DBounds.X.Max))
      {
        region = 12;
        inside = false;
        std::cout << "p1= " << p1 << " *** region 12" << std::endl;
      }
      else if (this->ApproxLT(p1[0], this->grid2DBounds.X.Min))
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

      viskores::Vec3f ptXY(p1[0], p1[1], 0.0);
      viskores::FloatDefault rxyVal, zxyVal, zind, zvalue;
      bool res0 = this->EvaluateLinear(ptXY, boutppField, boutppField.rxy, cells2D, rxyVal);
      bool res1 = this->EvaluateLinear(ptXY, boutppField, boutppField.zxy, cells2D, zxyVal);
      if (!res0 || !res1)
        break;

      zind = this->LinearInterpolate(boutppField.ZArray, boutppField.ZiArray, p1[2]);
      zvalue = this->LinearInterpolate(boutppField.ZiArray, boutppField.ZArray, zind);

#if TRACE_IO
      stepOut << xindIn << ", " << step << ", " << p1[0] << ", " << p1[1] << ", " << p1[2] << ", (" << zind << "), " << region << ", " << rxyVal
              << ", " << zxyVal << ", " << zvalue << std::endl;
#endif
      //convert to cartesian.
      pCart0 = pCart1;
      viskores::Vec3f vCart;
      pCart1 = this->ToCartesian(p1, boutppField, cells2D);

      //point crosses the X=0 plane.
      if (step > 1 && pCart0[0] * pCart1[0] < 0 && pCart0[0] > 0)
      {
        viskores::FloatDefault puncStep;
        viskores::Vec2f psiThetaPt;
        auto puncPt = this->FindPuncture(step, p0, p1, boutppField, cells2D, puncStep, psiThetaPt);
        viskores::Vec4f rzThetaPsi(puncPt[1], puncPt[2], psiThetaPt[1], psiThetaPt[0]);
        resultRZThetaPsi.Set(puncOffset + numPunc, rzThetaPsi);
        puncID.Set(puncOffset + numPunc, ID);
        puncIndex.Set(puncOffset + numPunc, numPunc);

        //result.Set(puncOffset + numPunc, puncPt);
        //resultStep.Set(puncOffset + numPunc, puncStep);
        //resultPsiTheta.Set(puncOffset + numPunc, psiThetaPt);
#if TRACE_IO
        rootPuncOut << xindIn << ", " << puncStep << ", " << puncPt[0] << ", " << puncPt[1] << ", " << puncPt[2] << std::endl;
        puncOut << xindIn << ", " << puncStep << ", " << puncPt[0] << ", " << puncPt[1] << ", " << puncPt[2] << ", 0.0, 0.0, 0.0" << std::endl;
#endif

        //puncIndex.Set(puncOffset + numPunc, ID);
        //validResult.Set(puncOffset + numPunc, true);
        numPunc++;
      }

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
#endif
  }

private:
  VISKORES_EXEC int GetRegion(const viskores::FloatDefault& xind, const viskores::FloatDefault& y) const
  {
    int region = 0;

    auto xxx = static_cast<viskores::FloatDefault>(this->ixsep1 - 1) + 0.5f;
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
                                             viskores::FloatDefault& stepPunc,
                                             viskores::Vec2f& psiTheta) const
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
      viskores::Vec3f cart_mid = this->ToCartesian(p_mid, boutppField, cells, psiTheta);
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
    auto p = this->ToCartesian(p_mid, boutppField, cells, psiTheta);
    return p;
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
    viskores::Vec2f psiTheta;
    return this->ToCartesian(pt, boutppField, cells, psiTheta);
  }

  template <typename BoutppFieldType, typename CellSetType>
  VISKORES_EXEC viskores::Vec3f ToCartesian(const viskores::Vec3f& pt,
                                            BoutppFieldType& boutppField,
                                            const CellSetType& cells,
                                            viskores::Vec2f& psiTheta) const
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

    psiTheta[0] = this->EvaluateCubic1D(ptXY, "psixy");
    psiTheta[1] = this->CalculateTheta(ptXY);

    return viskores::Vec3f(x3d, y3d, z3d);
  }

  VISKORES_EXEC viskores::FloatDefault interp1_linear(const viskores::FloatDefault* x,
                                                      const viskores::FloatDefault* y,
                                                      viskores::Id n,
                                                      viskores::FloatDefault xq) const
  {
    // Guard cases
    if (n <= 0)
      return viskores::FloatDefault{ 0 };
    if (n == 1)
      return y[0];

    // Clamp to ends
    if (xq <= x[0])
      return y[0];
    if (xq >= x[n - 1])
      return y[n - 1];

    // Binary search for i such that x[i] <= xq <= x[i+1]
    viskores::Id lo = 0;
    viskores::Id hi = n - 1;
    // Invariant: interval [lo, hi] contains xq and hi - lo >= 1
    while (hi - lo > 1)
    {
      const viskores::Id mid = lo + (hi - lo) / 2;
      if (xq >= x[mid])
        lo = mid;
      else
        hi = mid;
    }
    // Now interpolate between lo and lo+1
    const viskores::FloatDefault x0 = x[lo];
    const viskores::FloatDefault x1 = x[lo + 1];
    const viskores::FloatDefault y0 = y[lo];
    const viskores::FloatDefault y1 = y[lo + 1];

    const viskores::FloatDefault dx = x1 - x0;
    if (dx == viskores::FloatDefault{ 0 })
    {
      // Degenerate segment: fallback to left value
      return y0;
    }
    const viskores::FloatDefault t = (xq - x0) / dx;
    return (viskores::FloatDefault(1) - t) * y0 + t * y1;
  }

  VISKORES_EXEC viskores::FloatDefault CalculateTheta(const viskores::Vec3f& pt) const
  {
    auto yind_tmp = pt[1];

    auto theta = this->interp1_linear(this->yiarray_cfr.data(), this->theta_cfr.data(), this->yiarray_cfr.size(), yind_tmp);
    return theta;
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
    res[0] = this->CubicEval(this->ds3D, dxdy, pt);
    res[1] = this->CubicEval(this->ds3D, dzdy, pt);
    return res;
  }

  VISKORES_EXEC viskores::FloatDefault CubicEval(const viskores::cont::DataSet& ds, const std::string& fieldName, const viskores::Vec3f& pt) const
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

    auto val = this->TricubicSampleRectilinear(array, xCoords, yCoords, zCoords, pt);
    return val;
  }

  VISKORES_EXEC viskores::FloatDefault TricubicSampleRectilinear(const viskores::cont::ArrayHandle<viskores::FloatDefault>& _data,
                                                                 const viskores::cont::ArrayHandle<viskores::FloatDefault>& _xCoords,
                                                                 const viskores::cont::ArrayHandle<viskores::FloatDefault>& _yCoords,
                                                                 const viskores::cont::ArrayHandle<viskores::FloatDefault>& _zCoords,
                                                                 const viskores::Vec3f& pt) const
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
    viskores::Id iu = this->FindIndex(_xCoords, x);
    viskores::Id iv = this->FindIndex(_yCoords, y);
    viskores::Id iw = (nz >= 4 ? this->FindIndex(_zCoords, z) : 0);

    // If we don’t have 4 Z‐levels, do 2D bicubic in X–Y only
    if (nz < 4)
    {
      // --- bicubic: gather a 4×4 patch in X–Y at single k = clamp(iw,0,ny−1) ---
      double P2d[4 * 4];
      for (viskores::Id jj = 0; jj < 4; ++jj)
      {
        viskores::Id j = this->Clamp(iv - 1 + jj, ny);
        for (viskores::Id ii = 0; ii < 4; ++ii)
        {
          viskores::Id i = this->Clamp(iu - 1 + ii, nx);
          // flatten (i,j, 0) → data index
          P2d[jj * 4 + ii] = data.Get((j * nx) + i);
        }
      }
      // 3) bicubic along X → C2[4]
      double C2[4];
      for (int jj = 0; jj < 4; ++jj)
      {
        auto x0 = xCoords.Get(iu - 1), x1 = xCoords.Get(iu), x2 = xCoords.Get(iu + 1), x3 = xCoords.Get(iu + 2);
        C2[jj] = this->CubicInterpolateNonUniform(x0, x1, x2, x3, P2d[jj * 4 + 0], P2d[jj * 4 + 1], P2d[jj * 4 + 2], P2d[jj * 4 + 3], x);
      }
      // 4) bicubic along Y on C2 → final
      auto y0 = yCoords.Get(iv - 1), y1 = yCoords.Get(iv), y2 = yCoords.Get(iv + 1), y3 = yCoords.Get(iv + 2);
      return this->CubicInterpolateNonUniform(y0, y1, y2, y3, C2[0], C2[1], C2[2], C2[3], y);
    }

    // --- otherwise do full tricubic in X–Y–Z ---
    double P[4 * 4 * 4];
    for (viskores::Id kk = 0; kk < 4; ++kk)
    {
      viskores::Id k = this->Clamp(iw - 1 + kk, nz);
      for (viskores::Id jj = 0; jj < 4; ++jj)
      {
        viskores::Id j = this->Clamp(iv - 1 + jj, ny);
        for (viskores::Id ii = 0; ii < 4; ++ii)
        {
          viskores::Id i = this->Clamp(iu - 1 + ii, nx);
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
        Cbuf[kk * 4 + jj] = this->CubicInterpolateNonUniform(x0, x1, x2, x3, p0, p1, p2, p3, x);
      }

    // interpolate in Y → D[4]
    viskores::FloatDefault D[4];
    for (int kk = 0; kk < 4; ++kk)
    {
      auto y0 = yCoords.Get(iv - 1), y1 = yCoords.Get(iv), y2 = yCoords.Get(iv + 1), y3 = yCoords.Get(iv + 2);
      D[kk] = this->CubicInterpolateNonUniform(y0, y1, y2, y3, Cbuf[kk * 4 + 0], Cbuf[kk * 4 + 1], Cbuf[kk * 4 + 2], Cbuf[kk * 4 + 3], y);
    }

    // interpolate in Z
    auto z0 = zCoords.Get(iw - 1), z1 = zCoords.Get(iw), z2 = zCoords.Get(iw + 1), z3 = zCoords.Get(iw + 2);
    return this->CubicInterpolateNonUniform(z0, z1, z2, z3, D[0], D[1], D[2], D[3], z);
  }

  VISKORES_EXEC viskores::FloatDefault CubicInterpolateNonUniform(viskores::FloatDefault x0,
                                                                  viskores::FloatDefault x1,
                                                                  viskores::FloatDefault x2,
                                                                  viskores::FloatDefault x3,
                                                                  viskores::FloatDefault p0,
                                                                  viskores::FloatDefault p1,
                                                                  viskores::FloatDefault p2,
                                                                  viskores::FloatDefault p3,
                                                                  viskores::FloatDefault x) const
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

  VISKORES_EXEC viskores::Id Clamp(viskores::Id i, viskores::Id n) const
  {
    if (i < 0)
      return 0;
    else if (i >= n)
      return n - 1;
    else
      return i;
  }

  VISKORES_EXEC viskores::Id FindIndex(const viskores::cont::ArrayHandle<viskores::FloatDefault>& _coords, viskores::FloatDefault val) const
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

  VISKORES_EXEC bool ApproxLT(viskores::FloatDefault a, viskores::FloatDefault b, viskores::FloatDefault tol = 1e-6) const
  {
    if (a < (b - tol))
      return true;
    return false;
  }
  VISKORES_EXEC bool ApproxGT(viskores::FloatDefault a, viskores::FloatDefault b, viskores::FloatDefault tol = 1e-6) const
  {
    if (a > (b + tol))
      return true;
    return false;
  }

  viskores::Id MaxPunctures;
  viskores::Id MaxSteps;
  viskores::FloatDefault xMin, xMax;

public:
  viskores::cont::DataSet ds3D, ds2D;
  viskores::Bounds grid3DBounds, grid2DBounds;
  viskores::Vec2f Center;
  int ny;
  int nypf1, nypf2;
  int ixsep1, ixsep2;
  int divertor;
  int direction;

  std::vector<viskores::FloatDefault> yiarray_cfr, theta_cfr;
};
