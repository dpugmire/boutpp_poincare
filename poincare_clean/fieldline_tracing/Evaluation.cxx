#include "Evaluation.h"


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
vtkm::filter::resampling::Probe probe;

vtkm::cont::DataSetBuilderExplicit builder;
std::vector<vtkm::Vec3f> pts = { pt };
auto inPts = builder.Create(pts, {vtkm::CELL_SHAPE_VERTEX}, {1}, {0});

probe.SetGeometry(inPts);
auto output = probe.Execute(dataset);

vtkm::cont::ArrayHandle<vtkm::FloatDefault> fieldArray;
output.GetField(fieldname).GetData().AsArrayHandle(fieldArray);

return fieldArray.ReadPortal().Get(0);

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