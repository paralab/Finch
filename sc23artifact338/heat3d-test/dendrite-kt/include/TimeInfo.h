//
// Modified by boshun on 06/05/19.
//

#ifndef DENDRITE2_0_TIMEINFO_H
#define DENDRITE2_0_TIMEINFO_H
#include <vector>
#include <cassert>
#include <talyfem/utils/utils.h>

/**
 * Variable time step support, user will provide a std::vector<double> of dt and a same size of totalT
 * For example, dt = [0.5, 1], totalT = [5, 10] means that the first 5 seconds is using dt = 0.5 and
 * 5-10 is using dt = 1.
 */
class TimeInfo {
 public:
  TimeInfo(double t, std::vector<double> &dt, std::vector<double> &totalT)
      : currentTime_(t), dt_(dt), totalT_(totalT) {
    assert(dt.size() == totalT.size());
    currentStepNumber_ = 0;
  }

  /**
   * find out the current timestep based on the current time
   * @return
   */
  inline double getCurrentStep() const {
    for (int idx = totalT_.size() - 1; idx > 0; idx--) {
      if (currentTime_ >= totalT_.at(idx - 1)) {
        return dt_.at(idx);
      }
    }
    return *dt_.begin();
  }

  inline void print() const {
    TALYFEMLIB::PrintInfo("Timestep ", getTimeStepNumber(), " - ", getCurrentTime());
  }

  inline void increment() {
    currentStepNumber_++;
    currentTime_ += getCurrentStep();
  }

  inline void setTimeStepNumber(unsigned int currentstep) {
    currentStepNumber_ = currentstep;
  }

  inline void setCurrentTime(double current) {
    currentTime_ = current;
  }

  unsigned int getTimeStepNumber() const {
    return currentStepNumber_;
  }

  double getCurrentTime() const {
    return currentTime_;
  }

  double getEndTime() const {
    return *(totalT_.end() - 1);
  }

 private:
  // timestep
  unsigned int currentStepNumber_;
  // real time
  double currentTime_;
  // step control
  std::vector<double>& dt_;
  std::vector<double>& totalT_;
};

#endif //DENDRITE2_0_TIMEINFO_H
