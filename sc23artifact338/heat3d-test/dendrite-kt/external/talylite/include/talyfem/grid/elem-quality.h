
#ifndef GRID_ELEM_QUALITY_H_
#define GRID_ELEM_QUALITY_H_

#include <stdint.h>

namespace TALYFEMLIB {

/**
 * List of quality metrics known to TalyFEM
 */
enum QualityMetric {
    kVolume,
    kAngle,
    kFaceArea,
    kAspectRatio
};

/**
 * Often, metrics have multiple values per element. These metrics take an
 * additional "type" parameter to pick the min, max, avg.
 */
enum QualityMetricType {
    kMin,
    kMax,
    kAvg
};

};

#endif