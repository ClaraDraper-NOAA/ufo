/*
 * (C) Copyright 2026 NOAA/NCEP/OMD
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ModelLevelsBelowSurface.h"

#include <cmath>
#include <cstddef>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

static ObsFunctionMaker<ModelLevelsBelowSurface>
  makerModelLevelsBelowSurface_("ModelLevelsBelowSurface");

// -----------------------------------------------------------------------------

ModelLevelsBelowSurface::ModelLevelsBelowSurface(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::trace() << "ModelLevelsBelowSurface constructor" << std::endl;
  options_.validateAndDeserialize(conf);

  invars_ += options_.obsPressure.value();
  invars_ += Variable("GeoVaLs/air_pressure_at_surface");
  invars_ += Variable("GeoVaLs/air_pressure");
}

// -----------------------------------------------------------------------------

ModelLevelsBelowSurface::~ModelLevelsBelowSurface() {
  oops::Log::trace() << "ModelLevelsBelowSurface destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void ModelLevelsBelowSurface::compute(const ObsFilterData & in,
                                      ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "ModelLevelsBelowSurface compute start" << std::endl;

  const float missing = util::missingValue<float>();

  const size_t nlocs = in.nlocs();
  const Variable modPrsVar("GeoVaLs/air_pressure");
  const size_t nlevs = in.nlevs(modPrsVar);

  std::vector<float> obsPressure(nlocs);
  in.get(options_.obsPressure.value(), obsPressure);
  std::vector<float> modSfcPressure(nlocs);
  in.get(Variable("GeoVaLs/air_pressure_at_surface"), modSfcPressure);

  // get depth of the lowest model layer.
  std::vector<float> modPrs1(nlocs), modPrs2(nlocs);
  in.get(modPrsVar, nlevs - 1, modPrs1);
  in.get(modPrsVar, nlevs - 2, modPrs2);

  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (obsPressure[jj] == missing || modSfcPressure[jj] == missing ||
        modPrs1[jj] == missing || modPrs2[jj] == missing) {
      out[0][jj] = missing;
    } else if (obsPressure[jj] > modSfcPressure[jj]) {
      out[0][jj] = std::log(obsPressure[jj] / modSfcPressure[jj]) /
                     std::log(modPrs1[jj] / modPrs2[jj]);
    } else {
      out[0][jj] = 0.0f;
    }
  }
  oops::Log::trace() << "ModelLevelsBelowSurface compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ModelLevelsBelowSurface::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
