#include "Evaluation.h"
#include <vtkm/cont/CellLocatorUniformGrid.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/exec/CellInterpolate.h>
#include <vtkm/filter/flow/worklet/CellInterpolationHelper.h>
#include <vtkm/filter/resampling/Probe.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/ArrayHandleConstant.h>
#include <vtkm/cont/ArrayHandleCounting.h>


namespace
{
class LocatorWorklet : public vtkm::worklet::WorkletMapField
{
  public:
    using ControlSignature = void(FieldIn, WholeArrayIn, ExecObject, ExecObject, FieldOut);
    using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5);

    template <typename T, typename FieldType, typename LocatorType, typename InterpolationHelperType, typename U>
    VTKM_EXEC void operator()(const vtkm::Id& idx,
                              const T& pt,
                              const FieldType& field,
                              const LocatorType& locator,
                              const InterpolationHelperType& ipHelper,
                              U &result) const
    {
        vtkm::Id cellId;
        vtkm::Vec3f parametric;
        locator.FindCell(pt, cellId, parametric);

        vtkm::UInt8 cellShape;
        vtkm::IdComponent nVerts;
        vtkm::VecVariable<vtkm::Id, 8> ptIndices;
        vtkm::VecVariable<vtkm::FloatDefault, 8> fieldValues;
        ipHelper.GetCellInfo(cellId, cellShape, nVerts, ptIndices);
        //field.GetValue(0);
        //vtkm::exec::CellInterpolate(fieldValues, parametric, cellShape, nVerts, ptIndices, result);
        //vtkm::exec::CellInterpolate(ptIndices, fieldValues, parametric, cellShape, result);
    }
};

class LinearInterpolationWorklet : public vtkm::worklet::WorkletMapField
{
public:
    LinearInterpolationWorklet(vtkm::Id size) : Size(size) {}

    using ControlSignature = void(FieldIn val, WholeArrayIn x, WholeArrayIn y,FieldOut result);
    using ExecutionSignature = void(_1, _2, _3, _4);

    template <typename XArrayType, typename YArrayType, typename ResultType>
    VTKM_EXEC void operator()(const vtkm::FloatDefault& val,
                              const XArrayType& x,
                              const YArrayType& y,
                              ResultType& result) const
    {
      // Find the interval.
      vtkm::Id idx = -1;
      for (vtkm::Id i = 0; i < this->Size-1; i++)
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
  vtkm::Id Size;
};
} //namespace

vtkm::FloatDefault scalarField1DEval(const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& x,
    const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& y,
    const vtkm::FloatDefault& val)
{
  std::vector<vtkm::FloatDefault> vals;
  vals.push_back(val);
  return scalarField1DEval(x, y, vals);
}

vtkm::FloatDefault scalarField1DEval(const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& x,
                                     const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& y,
                                     const std::vector<vtkm::FloatDefault>& vals)
{
    vtkm::cont::Invoker invoker;

    auto worklet = LinearInterpolationWorklet( x.GetNumberOfValues());

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> result;
    auto input = vtkm::cont::make_ArrayHandle<vtkm::FloatDefault>(vals, vtkm::CopyFlag::On);
    invoker(worklet, input, x, y, result);

    return result.ReadPortal().Get(0);
}

vtkm::FloatDefault scalarField2DEval(const vtkm::cont::DataSet& dataset,
                                     const std::string& fieldname,
                                     const vtkm::Vec3f& pt)
{
    return scalarField3DEval(dataset, fieldname, pt);
}



#if 0
vtkm::FloatDefault scalarField2DEval(const vtkm::cont::DataSet& dataset,
                                     const std::string& fieldname,
                                     const std::vector<vtkm::Vec3f>& pts)
{
  return scalarField3DEval(dataset, fieldname, pts);
  /*
  vtkm::filter::resampling::Probe probe;

  vtkm::cont::DataSetBuilderExplicit builder;
  auto inPts = builder.Create(pts, { vtkm::CELL_SHAPE_VERTEX }, { 1 }, { 0 });

  probe.SetGeometry(inPts);
  auto output = probe.Execute(dataset);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> fieldArray;
  output.GetField(fieldname).GetData().AsArrayHandle(fieldArray);

  return fieldArray.ReadPortal().Get(0);
  */

#if 0
    vtkm::cont::CellLocatorUniformGrid locator;
    locator.SetCoordinates(dataset.GetCoordinateSystem());
    locator.SetCellSet(dataset.GetCellSet());
    locator.Update();

    vtkm::cont::CellInterpolationHelper interpolationHelper(dataset.GetCellSet());
    vtkm::cont::Invoker invoker;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> input = vtkm::cont::make_ArrayHandle<vtkm::Vec3f>(pts, vtkm::CopyFlag::Off);
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> output;
    invoker(LocatorWorklet{}, input, dataset.GetField(fieldname), locator, interpolationHelper, output);
    //  interpolationHelper.
#endif
}
#endif

vtkm::FloatDefault
scalarField3DEval(const vtkm::cont::DataSet& dataset,
                  const std::string& fieldname,
                  const vtkm::Vec3f& pt)
{
    vtkm::cont::DataSetBuilderExplicit builder;
    vtkm::Id num = 1;

    if (num != 1)
      throw std::runtime_error("Only 1 point is supported");

    std::vector<vtkm::Vec3f> pts = {pt};
    std::vector<vtkm::Id> ptIds;
    std::vector<vtkm::IdComponent> numPts;
    std::vector<vtkm::UInt8> cellTypes;
    for (vtkm::Id i = 0; i < num; i++)
    {
        ptIds.push_back(i);
        numPts.push_back(1);
        cellTypes.push_back(vtkm::CELL_SHAPE_VERTEX);
    }

    auto samplePts = builder.Create(pts, cellTypes, numPts, ptIds);
    //samplePts.PrintSummary(std::cout);

    vtkm::filter::resampling::Probe probe;
    probe.SetGeometry(samplePts);
    auto output = probe.Execute(dataset);
    /*
    vtkm::io::VTKDataSetWriter writer("grid3D.vtk");
    writer.WriteDataSet(dataset);
    writer = vtkm::io::VTKDataSetWriter("samplePts.vtk");
    writer.WriteDataSet(samplePts);
    output.PrintSummary(std::cout);
    */

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> fieldArray;
    output.GetField(fieldname).GetData().AsArrayHandle(fieldArray);

    return fieldArray.ReadPortal().Get(0);
}

vtkm::FloatDefault
scalarField3DEvalFwdAvg(const vtkm::cont::DataSet& dataset,
                  const std::string& fieldname,
                  const vtkm::Vec3f& pt)
{
  auto v0 = scalarField3DEval(dataset, fieldname, pt);
  auto v1 = scalarField3DEval(dataset, fieldname, {pt[0], pt[1]+1, pt[2]});

  return (v0+v1) / 2.0;
}