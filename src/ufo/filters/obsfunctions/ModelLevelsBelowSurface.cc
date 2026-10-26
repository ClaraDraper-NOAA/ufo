/*
 * (C) Copyright 2026 NOAA/NCEP/EMC
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

  const std::size_t nlocs = in.nlocs();
  const Variable prsVar("GeoVaLs/air_pressure");
  const std::size_t nlevs = in.nlevs(prsVar);
  ASSERT(nlevs > 1);

  std::vector<float> obsPressure(nlocs);
  in.get(options_.obsPressure.value(), obsPressure);
  std::vector<float> sfcPressure(nlocs);
  in.get(Variable("GeoVaLs/air_pressure_at_surface"), sfcPressure);

  // The two lowest model levels (GeoVaLs are ordered from the top of the atmosphere down).
  std::vector<float> prs1(nlocs), prs2(nlocs);
  in.get(prsVar, nlevs - 1, prs1);
  in.get(prsVar, nlevs - 2, prs2);

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    out[0][jloc] = missing;
    if (obsPressure[jloc] == missing || sfcPressure[jloc] == missing ||
        prs1[jloc] == missing || prs2[jloc] == missing) continue;
    // The model surface is below the lowest model level.
    ASSERT(sfcPressure[jloc] >= prs1[jloc] && prs1[jloc] > prs2[jloc]);
    if (obsPressure[jloc] > sfcPressure[jloc]) {
      out[0][jloc] = std::log(obsPressure[jloc] / sfcPressure[jloc]) /
                     std::log(prs1[jloc] / prs2[jloc]);
    } else {
      out[0][jloc] = 0.0f;
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
