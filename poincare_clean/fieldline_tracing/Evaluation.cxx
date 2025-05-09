#include "Evaluation.h"
#include <viskores/cont/ArrayHandleConstant.h>
#include <viskores/cont/ArrayHandleCounting.h>
#include <viskores/cont/CellLocatorUniformGrid.h>
#include <viskores/cont/DataSetBuilderExplicit.h>
#include <viskores/cont/DataSetBuilderUniform.h>
#include <viskores/exec/CellInterpolate.h>
#include <viskores/filter/flow/worklet/CellInterpolationHelper.h>
#include <viskores/filter/resampling/Probe.h>
#include <viskores/io/VTKDataSetWriter.h>
#include <viskores/worklet/WorkletMapField.h>


namespace
{
class LocatorWorklet : public viskores::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn, WholeArrayIn, ExecObject, ExecObject, FieldOut);
  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5);

  template <typename T, typename FieldType, typename LocatorType, typename InterpolationHelperType, typename U>
  VISKORES_EXEC void operator()(const viskores::Id& idx,
                                const T& pt,
                                const FieldType& field,
                                const LocatorType& locator,
                                const InterpolationHelperType& ipHelper,
                                U& result) const
  {
    viskores::Id cellId;
    viskores::Vec3f parametric;
    locator.FindCell(pt, cellId, parametric);

    viskores::UInt8 cellShape;
    viskores::IdComponent nVerts;
    viskores::VecVariable<viskores::Id, 8> ptIndices;
    viskores::VecVariable<viskores::FloatDefault, 8> fieldValues;
    ipHelper.GetCellInfo(cellId, cellShape, nVerts, ptIndices);
    //field.GetValue(0);
    //viskores::exec::CellInterpolate(fieldValues, parametric, cellShape, nVerts, ptIndices, result);
    //viskores::exec::CellInterpolate(ptIndices, fieldValues, parametric, cellShape, result);
  }
};

class LinearInterpolationWorklet : public viskores::worklet::WorkletMapField
{
public:
  LinearInterpolationWorklet(viskores::Id size)
    : Size(size)
  {
  }

  using ControlSignature = void(FieldIn val, WholeArrayIn x, WholeArrayIn y, FieldOut result);
  using ExecutionSignature = void(_1, _2, _3, _4);

  template <typename XArrayType, typename YArrayType, typename ResultType>
  VISKORES_EXEC void operator()(const viskores::FloatDefault& val, const XArrayType& x, const YArrayType& y, ResultType& result) const
  {
    // Find the interval.
    viskores::Id idx = -1;
    for (viskores::Id i = 0; i < this->Size - 1; i++)
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

    result = y1 + (y2 - y1) * (val - x1) / (x2 - x1);
  }

private:
  viskores::Id Size;
};
} //namespace

viskores::FloatDefault scalarField1DEval(const viskores::cont::ArrayHandle<viskores::FloatDefault>& x,
                                         const viskores::cont::ArrayHandle<viskores::FloatDefault>& y,
                                         const viskores::FloatDefault& val)
{
  std::vector<viskores::FloatDefault> vals;
  vals.push_back(val);
  return scalarField1DEval(x, y, vals);
}

viskores::FloatDefault scalarField1DEval(const viskores::cont::ArrayHandle<viskores::FloatDefault>& x,
                                         const viskores::cont::ArrayHandle<viskores::FloatDefault>& y,
                                         const std::vector<viskores::FloatDefault>& vals)
{
  viskores::cont::Invoker invoker;

  auto worklet = LinearInterpolationWorklet(x.GetNumberOfValues());

  viskores::cont::ArrayHandle<viskores::FloatDefault> result;
  auto input = viskores::cont::make_ArrayHandle<viskores::FloatDefault>(vals, viskores::CopyFlag::On);
  invoker(worklet, input, x, y, result);

  return result.ReadPortal().Get(0);
}

viskores::FloatDefault scalarField2DEval(const viskores::cont::DataSet& dataset, const std::string& fieldname, const viskores::Vec3f& pt)
{
  return scalarField3DEval(dataset, fieldname, pt);
}



#if 0
viskores::FloatDefault scalarField2DEval(const viskores::cont::DataSet& dataset,
                                     const std::string& fieldname,
                                     const std::vector<viskores::Vec3f>& pts)
{
  return scalarField3DEval(dataset, fieldname, pts);
  /*
  viskores::filter::resampling::Probe probe;

  viskores::cont::DataSetBuilderExplicit builder;
  auto inPts = builder.Create(pts, { viskores::CELL_SHAPE_VERTEX }, { 1 }, { 0 });

  probe.SetGeometry(inPts);
  auto output = probe.Execute(dataset);

  viskores::cont::ArrayHandle<viskores::FloatDefault> fieldArray;
  output.GetField(fieldname).GetData().AsArrayHandle(fieldArray);

  return fieldArray.ReadPortal().Get(0);
  */

#if 0
    viskores::cont::CellLocatorUniformGrid locator;
    locator.SetCoordinates(dataset.GetCoordinateSystem());
    locator.SetCellSet(dataset.GetCellSet());
    locator.Update();

    viskores::cont::CellInterpolationHelper interpolationHelper(dataset.GetCellSet());
    viskores::cont::Invoker invoker;
    viskores::cont::ArrayHandle<viskores::Vec3f> input = viskores::cont::make_ArrayHandle<viskores::Vec3f>(pts, viskores::CopyFlag::Off);
    viskores::cont::ArrayHandle<viskores::FloatDefault> output;
    invoker(LocatorWorklet{}, input, dataset.GetField(fieldname), locator, interpolationHelper, output);
    //  interpolationHelper.
#endif
}
#endif

viskores::FloatDefault scalarField3DEval(const viskores::cont::DataSet& dataset, const std::string& fieldname, const viskores::Vec3f& pt)
{
  viskores::cont::DataSetBuilderExplicit builder;
  viskores::Id num = 1;

  if (num != 1)
    throw std::runtime_error("Only 1 point is supported");

  std::vector<viskores::Vec3f> pts = { pt };
  std::vector<viskores::Id> ptIds;
  std::vector<viskores::IdComponent> numPts;
  std::vector<viskores::UInt8> cellTypes;
  for (viskores::Id i = 0; i < num; i++)
  {
    ptIds.push_back(i);
    numPts.push_back(1);
    cellTypes.push_back(viskores::CELL_SHAPE_VERTEX);
  }

  auto samplePts = builder.Create(pts, cellTypes, numPts, ptIds);
  //samplePts.PrintSummary(std::cout);

  viskores::filter::resampling::Probe probe;
  probe.SetGeometry(samplePts);
  auto output = probe.Execute(dataset);
  /*
    viskores::io::VTKDataSetWriter writer("grid3D.vtk");
    writer.WriteDataSet(dataset);
    writer = viskores::io::VTKDataSetWriter("samplePts.vtk");
    writer.WriteDataSet(samplePts);
    output.PrintSummary(std::cout);
    */

  viskores::cont::ArrayHandle<viskores::FloatDefault> fieldArray;
  output.GetField(fieldname).GetData().AsArrayHandle(fieldArray);

  return fieldArray.ReadPortal().Get(0);
}

viskores::FloatDefault scalarField3DEvalFwdAvg(const viskores::cont::DataSet& dataset, const std::string& fieldname, const viskores::Vec3f& pt)
{
  auto v0 = scalarField3DEval(dataset, fieldname, pt);
  auto v1 = scalarField3DEval(dataset, fieldname, { pt[0], pt[1] + 1, pt[2] });

  return (v0 + v1) / 2.0;
}