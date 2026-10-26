/*
 * (C) Copyright 2026 NOAA/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_MODELLEVELSBELOWSURFACE_H_
#define UFO_FILTERS_OBSFUNCTIONS_MODELLEVELSBELOWSURFACE_H_

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling the ModelLevelsBelowSurface ObsFunction
class ModelLevelsBelowSurfaceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelLevelsBelowSurfaceParameters, Parameters)

 public:
  oops::Parameter<Variable> obsPressure{"observation pressure",
    "Observed pressure (Pa)",
    Variable("MetaData/pressure"), this};
};

// -----------------------------------------------------------------------------

/// \brief Number of model levels an observation is below the model surface.
///
/// \details Returns ln(p_obs/p_sfc) / ln(p_1/p_2) where the observed pressure p_obs is greater
/// than the model surface pressure p_sfc (GeoVaLs air_pressure_at_surface), and zero otherwise.
/// p_1 and p_2 are the two lowest levels of GeoVaLs air_pressure. This is the GSI rlow term,
/// which the GSI uses to inflate the errors of observations below the model surface: the GSI
/// locates the observation and the model surface in units of model levels, extrapolating
/// linearly in ln(pressure) from the two lowest levels.
///
/// The output is missing where the observed pressure or the GeoVaLs are missing.
///
/// Example: add 8 K per model level below the model surface to the temperature error, as in the
/// GSI setupt.
///
///  - filter: Variable Assignment
///    assignments:
///    - name: DerivedMetaData/modelLevelsBelowSurface
///      type: float
///      function:
///        name: ObsFunction/ModelLevelsBelowSurface
///  - filter: Perform Action
///    filter variables:
///    - name: airTemperature
///    action:
///      name: assign error
///      error function:
///        name: ObsFunction/Arithmetic
///        options:
///          variables:
///          - name: ObsErrorData/airTemperature
///          - name: DerivedMetaData/modelLevelsBelowSurface
///          coefs: [1.0, 8.0]
///
class ModelLevelsBelowSurface : public ObsFunctionBase<float> {
 public:
  explicit ModelLevelsBelowSurface(const eckit::LocalConfiguration &);
  ~ModelLevelsBelowSurface();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ModelLevelsBelowSurfaceParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_MODELLEVELSBELOWSURFACE_H_
