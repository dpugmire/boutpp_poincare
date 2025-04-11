#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/CellLocatorUniformGrid.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/filter/flow/worklet/CellInterpolationHelper.h>
#include <vtkm/filter/resampling/Probe.h>
#include <vtkm/exec/CellInterpolate.h>

#ifndef _Evauluation_h_
#define _Evaluation_h_

vtkm::FloatDefault scalarField1DEval(const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& x,
    const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& y,
    const std::vector<vtkm::FloatDefault>& vals);

vtkm::FloatDefault scalarField2DEval(const vtkm::cont::DataSet& dataset,
    const std::string& fieldname,
    const vtkm::Vec3f& pt);


#endif // _Evauluation_h_