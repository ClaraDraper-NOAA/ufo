/*
 * (C) Copyright 2026 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorDensity.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "ioda/distribution/Distribution.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/GeometryCalculations.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorDensity> makerObsErrorFactorDensity_("ObsErrorFactorDensity");

namespace {

/// Search boxes: latitude bands of the horizontal mesh, each split into longitude boxes of about
/// the same size, and vertical layers of the vertical mesh, within a time interval.
struct SearchGrid {
  double dlat;    // band width (degrees)
  int nlat;       // number of bands
  int nlonmax;    // number of longitude boxes at the equator

  int band(double lat) const {
    return std::min(nlat - 1, std::max(0, static_cast<int>((lat + 90.0) / dlat)));
  }
  int nlon(int ilat) const {
    const double latc = -90.0 + (ilat + 0.5) * dlat;
    return std::max(1, static_cast<int>(std::round(nlonmax * std::cos(latc * Constants::deg2rad))));
  }
  int lonbox(double lon, int ilat) const {
    double lon360 = std::fmod(lon, 360.0);
    if (lon360 < 0.0) lon360 += 360.0;
    const int n = nlon(ilat);
    return std::min(n - 1, static_cast<int>(lon360 / 360.0 * n));
  }
  int64_t key(int64_t itime, int64_t ivert, int64_t ilat, int64_t ilon) const {
    // vertical index offset so that negative values are allowed
    return ((itime * 2000001 + (ivert + 1000000)) * nlat + ilat) * nlonmax + ilon;
  }
};

}  // namespace

// -----------------------------------------------------------------------------

ObsErrorFactorDensity::ObsErrorFactorDensity(const eckit::LocalConfiguration & config)
  : invars_() {
  oops::Log::trace() << "ObsErrorFactorDensity constructor" << std::endl;
  options_.validateAndDeserialize(config);
  invars_ += Variable(options_.testQCflag.value() + "/" + options_.variable.value());
  invars_ += Variable("MetaData/latitude");
  invars_ += Variable("MetaData/longitude");
  invars_ += Variable("MetaData/dateTime");
  invars_ += options_.verticalCoordinate.value();
}

// -----------------------------------------------------------------------------

ObsErrorFactorDensity::~ObsErrorFactorDensity() {
  oops::Log::trace() << "ObsErrorFactorDensity destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorFactorDensity::compute(const ObsFilterData & data,
                                    ioda::ObsDataVector<float> & factor) const {
  oops::Log::trace() << "ObsErrorFactorDensity compute start" << std::endl;
  const ioda::ObsSpace & obsdb = data.obsspace();
  const size_t nlocs = obsdb.nlocs();
  const float missing = util::missingValue<float>();

  std::vector<int> qcflag(nlocs);
  data.get(Variable(options_.testQCflag.value() + "/" + options_.variable.value()), qcflag);
  std::vector<float> lat(nlocs), lon(nlocs), vcoord(nlocs);
  data.get(Variable("MetaData/latitude"), lat);
  data.get(Variable("MetaData/longitude"), lon);
  data.get(options_.verticalCoordinate.value(), vcoord);
  std::vector<util::DateTime> datetime(nlocs);
  data.get(Variable("MetaData/dateTime"), datetime);

  const double hmesh = 1000.0 * options_.horizontalMesh.value();  // m
  const double vmesh = options_.verticalMesh.value();
  const int64_t tmesh = options_.timeMesh.value().toSeconds();
  const float threshold = options_.densityThreshold.value();
  const util::DateTime tstart = obsdb.windowStart();

  SearchGrid grid;
  grid.dlat = hmesh / (Constants::mean_earth_rad_m * Constants::deg2rad);
  grid.nlat = std::max(1, static_cast<int>(std::ceil(180.0 / grid.dlat)));
  grid.nlonmax = std::max(1, static_cast<int>(std::ceil(360.0 / grid.dlat)));

  // Observations that are counted: passing QC, with a location, time and vertical coordinate.
  std::vector<int> used(nlocs, 0), itime(nlocs, 0);
  std::vector<double> dlat(nlocs), dlon(nlocs), dvert(nlocs);
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    dlat[iloc] = lat[iloc];
    dlon[iloc] = lon[iloc];
    dvert[iloc] = vcoord[iloc];
    if (qcflag[iloc] != 0 || lat[iloc] == missing || lon[iloc] == missing ||
        vcoord[iloc] == missing || hmesh <= 0.0 || vmesh <= 0.0 || tmesh <= 0) continue;
    used[iloc] = 1;
    const int64_t dt = (datetime[iloc] - tstart).toSeconds();
    itime[iloc] = static_cast<int>(std::max<int64_t>(0, dt / tmesh));
  }

  // Gather the counted observations from all tasks (each observation appears once).
  std::vector<int> gused(used), gtime(itime);
  std::vector<double> glat(dlat), glon(dlon), gvert(dvert);
  obsdb.distribution()->allGatherv(gused);
  obsdb.distribution()->allGatherv(gtime);
  obsdb.distribution()->allGatherv(glat);
  obsdb.distribution()->allGatherv(glon);
  obsdb.distribution()->allGatherv(gvert);

  // Index the counted observations by search box.
  std::unordered_map<int64_t, std::vector<size_t>> boxes;
  for (size_t j = 0; j < gused.size(); ++j) {
    if (!gused[j]) continue;
    const int ilat = grid.band(glat[j]);
    const int64_t ivert = static_cast<int64_t>(std::floor(gvert[j] / vmesh));
    boxes[grid.key(gtime[j], ivert, ilat, grid.lonbox(glon[j], ilat))].push_back(j);
  }

  // Density: the number of observations within a cylinder centred on each observation, with the
  // horizontal area of a square of the horizontal mesh size (radius hmesh/sqrt(pi)) and the height
  // of the vertical mesh size (+/- vmesh/2). Its volume is that of GSI's counting volume, and,
  // being centred on each observation, it has no box edges. For observations along lines (e.g.
  // aircraft tracks) it also counts like a box of the mesh sizes, on average.
  const double radius = hmesh / std::sqrt(M_PI);
  for (size_t ivar = 0; ivar < factor.nvars(); ++ivar) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) factor[ivar][iloc] = 1.0f;
  }
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (!used[iloc] || threshold <= 0.0f) continue;
    const int ilat0 = grid.band(dlat[iloc]);
    const int64_t ivert0 = static_cast<int64_t>(std::floor(dvert[iloc] / vmesh));
    double count = 0.0;
    for (int ilat = std::max(0, ilat0 - 1); ilat <= std::min(grid.nlat - 1, ilat0 + 1); ++ilat) {
      // longitude boxes within one mesh size of the observation in this band
      const int n = grid.nlon(ilat);
      const double latedge = std::min(89.999, std::max(std::abs(-90.0 + ilat * grid.dlat),
                                                       std::abs(-90.0 + (ilat + 1) * grid.dlat)));
      const double dlonmax = grid.dlat / std::max(1.0e-6, std::cos(latedge * Constants::deg2rad));
      const int nbox = (2.0 * dlonmax >= 360.0) ? n :
                       std::min(n, 2 * static_cast<int>(std::ceil(dlonmax / 360.0 * n)) + 1);
      const int ilonc = grid.lonbox(dlon[iloc], ilat);
      for (int k = 0; k < nbox; ++k) {
        const int ilon = (nbox == n) ? k : ((ilonc - nbox / 2 + k) % n + n) % n;
        for (int64_t ivert = ivert0 - 1; ivert <= ivert0 + 1; ++ivert) {
          const auto it = boxes.find(grid.key(itime[iloc], ivert, ilat, ilon));
          if (it == boxes.end()) continue;
          for (const size_t j : it->second) {
            if (std::abs(gvert[j] - dvert[iloc]) > 0.5 * vmesh) continue;
            if (haversine(dlat[iloc], dlon[iloc], glat[j], glon[j]) > radius) continue;
            count += 1.0;
          }
        }
      }
    }
    const float ratio = static_cast<float>(count) / threshold;
    for (size_t ivar = 0; ivar < factor.nvars(); ++ivar) {
      factor[ivar][iloc] = std::sqrt(std::max(1.0f, ratio));
    }
  }
  oops::Log::trace() << "ObsErrorFactorDensity compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorDensity::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
