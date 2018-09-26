#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

// Suppress warnings
#pragma GCC diagnostic push
#include "Eigen/Dense"
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic pop

class MeasurementPackage {
public:
  long long timestamp_;

  enum SensorType{
    LASER,
    RADAR
  } sensor_type_;

  Eigen::VectorXd raw_measurements_;
};

#endif /* MEASUREMENT_PACKAGE_H_ */
