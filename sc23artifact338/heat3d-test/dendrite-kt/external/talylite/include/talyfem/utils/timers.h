/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#ifndef UTILS_TIMERS_H_
#define UTILS_TIMERS_H_

// Defines Timer classes to measure time used by a block of code, including
// support for averaging time over repeated calls.
//
// These are not intended to be lightweight timers and should not be used in
// performance critical sections of code except for testing.
//
// Two classes are defined:
//   MPITimer - gives a portable time based on MPI_wtime() function
//   HiResTimer - high resolution timer that is less portable but should have
//                nanosecond precision on systems that support it.
//
//
// Example usage for the MPITimer class is shown below. Other timers follow
// the same structure.
//
// MPITimer timer("timer description");
// for (int i = 0; i < 100; i++) {
//   // do some boring stuff we don't want to time
//   timer.Start();
//   // do interesting stuff we want to time
//   timer.Stop();
// }
// timer.PrintAvgSeconds();
//

#include <time.h>

#include <string>
#include <vector>


/**
 * Base class for the timers.
 *
 * This is an interface and not intended to be used for actual objects.
 * Timers should use this as a base class.
 */
class TimerBase {
 public:
  TimerBase()
      : label_("Timer"),
        n_timings_(0),
        timer_running_(false) {
  }

  /**
   * Construct the timer with the given descriptive label
   *
   * @param label description of the timer
   */
  explicit TimerBase(std::string label)
      : label_(label),
        n_timings_(0),
        timer_running_(false) {
  }

  virtual ~TimerBase() { }

  /**
   * Starts a cycle of the timer.
   * If the timer is currently running, the elapsed time is added to the total.
   */
  virtual void Start() = 0;

  /**
   * ReStarts a cycle of the timer.
   * If the timer is currently running, the elapsed time is discarded and the
   * timer is started again. If the timer was not running, it is started.
   */
  virtual void ReStart() = 0;

  /**
   * Stops the current running of the timer and adds elasped time to total time.
   * If the timer is not currently running, this has no effect.
   */
  virtual void Stop() = 0;

  /**
   * Stops the current running of the timer and ignores the elasped time.
   * If the timer is not currently running, this has no effect.
   */
  inline void Discard() {
    timer_running_ = false;
  }

  /**
   * Resets all timer data except the timer label.
   */
  virtual void Reset() = 0;

  /**
   * Return a bool of whether the timer is running.
   *
   * @return true if the timer is currently running
   */
  inline bool timer_running() {
    return timer_running_;
  }

  /**
   * Set the label for the timer
   *
   * @param label The new label to assign to the timer.
   */
  inline void set_label(std::string label) {
    label_ = label;
  }

  /**
   * Returns the average time in seconds.
   * If the timer is currently running, the current timing is not included.
   *
   * @return the average time per call, in seconds
   */
  virtual double GetAvgTimeSeconds() const {
    return 0.0;
  }

  /**
   * Returns the total time in seconds.
   * If the timer is currently running, the current timing is not included.
   *
   * @return the total time spent over all calls, in seconds
   */
  virtual double GetTotalTimeSeconds() const {
    return 0.0;
  }

  /**
   * Returns the elapsed time since the start of the current timing iteration
   * (i.e. time since last call to Start() or ReStart()).
   * If the timer is not currently running, 0.0 is returned.
   *
   * @return time since the start of the current timing iteration.
   */
  virtual double GetElapsedTimeSeconds() const = 0;

  /**
   * Returns the time, in seconds, taken by the last timing iteration
   * (i.e. time between last calls of Start() and either Stop() or ReStart()).
   * If there was no previous timing iteration, this returns 0.0
   * This has no effect on any current timing.
   *
   * @return time taken by last measure cycle
   */
  virtual double GetLastTimeSeconds() const = 0;

  /**
   * Prints the average time in seconds.
   */
  void PrintAvgSeconds() const;

  /**
   * Prints the total time in seconds.
   * The total time is the sum of the time of all measurement cycles.
   */
  void PrintTotalTimeSeconds() const;

  /**
   * Prints the time, in seconds, taken by the last timing iteration
   * (i.e. time between last calls of Start() and either Stop() or ReStart()).
   * If there was no previous timing iteration, this prints 0.0
   */
  void PrintLastTimeSeconds() const;

  /**
   * Prints MPI global values for the average times in seconds.
   * Optionally, can also print min and max values.
   * The values printed are based on each process's average value.
   *
   * @param include_min_max whether to also print the min and max values
   */
  void PrintGlobalAverageSeconds(bool include_min_max = false) const;

  /**
   * Prints MPI global values for the average in seconds.
   * Optionally, can also print min and max values.
   * The values printed are based on each process's total value.
   *
   * @param include_min_max whether to also print the min and max values
   */
  void PrintGlobalTotalSeconds(bool include_min_max = false) const;

  /**
   * Prints an ASCII bar chart of how long each process took in this timer (in seconds).
   */
  void PrintGlobalChart() const;

  /**
   * Adds the given value to the total time.
   *
   * @param adjustment value in seconds to add to total time.
   */
  virtual void AddToTotalTime(double adjustment) = 0;

 protected:
  std::string label_;  ///< string describing the timer
  int n_timings_;  ///< number of timing sets that have been recorded
  bool timer_running_;  ///< whether the timer is currently running
};

/**
 * Timer class using MPI_Wtime
 *
 * See file comment block at top of file for example usage.
 */
class MPITimer : public TimerBase {
 public:
  MPITimer()
      : TimerBase(),
        start_time_(0.0),
        total_time_(0.0),
        last_time_(0.0) {
  }

  /**
   * Construct the timer with the given descriptive label
   *
   * @param label description of the timer
   */
  explicit MPITimer(std::string label)
      : TimerBase(label),
        start_time_(0.0),
        total_time_(0.0),
        last_time_(0.0) {
  }

  virtual ~MPITimer() {
  }

  /**
   * Starts a cycle of the timer.
   * If the timer is currently running, the elapsed time is added to the total.
   */
  virtual void Start();

  /**
   * ReStarts a cycle of the timer.
   * If the timer is currently running, the elapsed time is discarded and the
   * timer is started again. If the timer was not running, it is started.
   */
  virtual void ReStart();

  /**
   * Stops the current running of the timer and adds elasped time to total time.
   * If the timer is not currently running, this has no effect.
   */
  virtual void Stop();

  /**
   * Resets all timer data except the timer label.
   */
  virtual void Reset();

  /**
   * Returns the average time in seconds.
   * If the timer is currently running, the current timing is not included.
   *
   * @return the average time per call, in seconds
   */
  virtual double GetAvgTimeSeconds() const;

  /**
   * Returns the total time in seconds.
   * If the timer is currently running, the current timing is not included.
   *
   * @return the total time spent over all calls, in seconds
   */
  virtual double GetTotalTimeSeconds() const;

  /**
   * Returns the elapsed time since the start of the current timing iteration
   * (i.e. time since last call to Start() or ReStart()).
   * If the timer is not currently running, 0.0 is returned.
   *
   * @return time since the start of the current timing iteration.
   */
  virtual double GetElapsedTimeSeconds() const;

  /**
   * Returns the time, in seconds, taken by the last timing iteration
   * (i.e. time between last calls of Start() and either Stop() or ReStart()).
   * If there was no previous timing iteration, this returns 0.0
   * This has no effect on any current timing.
   *
   * @return time taken by last measure cycle
   */
  virtual double GetLastTimeSeconds() const;

  /**
   * Adds the given value to the total time.
   *
   * @param adjustment value in seconds to add to total time.
   */
  void AddToTotalTime(double adjustment);

 protected:
  /**
   * Returns the time in seconds since the start of the epoch.
   *
   * @return time since beginning of the epoch in seconds
   */
  virtual double GetCurrentTimeSeconds() const;

  double start_time_;  ///< time when the last timing began
  double total_time_;  ///< total time that has been measured
  double last_time_;  ///< time elasped in last cycle
};

#if !defined(__MACH__) && !defined(__APPLE__)
/**
 * High resolution timer class.
 *
 * This class implements a timer measuring the time elapsed for multiple calls
 * of a code block. This has nanosecond resolution on systems that support it.
 * Currently, this does not work on Mac because Macs don't support
 * clock_gettime.
 *
 * See file comment block at top of file for example usage.
 */
class HiResTimer : public TimerBase {
 public:
  HiResTimer()
      : TimerBase() {
    start_time_.tv_sec = 0.0;
    start_time_.tv_nsec = 0.0;
    total_time_.tv_sec = 0.0;
    total_time_.tv_nsec = 0.0;
    last_time_.tv_sec = 0.0;
    last_time_.tv_nsec = 0.0;
  }

  /**
   * Construct the timer with the given descriptive label
   *
   * @param label description of the timer
   */
  explicit HiResTimer(std::string label)
      : TimerBase(label) {
    start_time_.tv_sec = 0.0;
    start_time_.tv_nsec = 0.0;
    total_time_.tv_sec = 0.0;
    total_time_.tv_nsec = 0.0;
    last_time_.tv_sec = 0.0;
    last_time_.tv_nsec = 0.0;
  }

  ~HiResTimer() {
  }

  /**
   * Starts a cycle of the timer.
   * If the timer is currently running, the elapsed time is added to the total.
   */
  void Start();

  /**
   * ReStarts a cycle of the timer.
   * If the timer is currently running, the elapsed time is discarded and the
   * timer is started again. If the timer was not running, it is started.
   */
  void ReStart();

  /**
   * Stops the current running of the timer and adds elasped time to total time.
   * If the timer is not currently running, this has no effect.
   */
  void Stop();

  /**
   * Resets all timer data except the timer label.
   */
  void Reset();

  /**
   * Returns the average time in seconds.
   * If the timer is currently running, the current timing is not included.
   *
   * @return the average time per call, in seconds
   */
  double GetAvgTimeSeconds() const;

  /**
   * Returns the average time in nanoseconds.
   * If the timer is currently running, the current timing is not included.
   *
   * @return the average time per call, in nanoseconds
   */
  long GetAvgTimeNSeconds() const;  // NOLINT(runtime/int) - is long in time.h

  /**
   * Returns the total time in seconds.
   * If the timer is currently running, the current timing is not included.
   *
   * @return the total time spent over all calls, in seconds
   */
  double GetTotalTimeSeconds() const;

  /**
   * Returns the elapsed time since the start of the current timing iteration
   * (i.e. time since last call to Start() or ReStart()).
   * If the timer is not currently running, 0.0 is returned.
   *
   * @return time since the start of the current timing iteration.
   */
  double GetElapsedTimeSeconds() const;

  /**
   * Returns the time taken by the last timing iteration
   * (i.e. time between last calls of Start() and either Stop() or ReStart()).
   * If there was no previous timing iteration, this returns 0.0
   * This has no effect on any current timing.
   *
   * @return time taken by last measure cycle
   */
  timespec GetLastTime() const;

  /**
   * Returns the time, in seconds, taken by the last timing iteration
   * (i.e. time between last calls of Start() and either Stop() or ReStart()).
   * If there was no previous timing iteration, this returns 0.0
   * This has no effect on any current timing.
   *
   * @return time taken by last measure cycle
   */
  virtual double GetLastTimeSeconds() const;

  /**
   * Adds the given value to the total time.
   *
   * @param adjustment value in seconds to add to total time.
   */
  void AddToTotalTime(double adjustment);

  /**
   * Prints MPI global values for the average times in nanoseconds.
   * Optionally, can also print min and max values.
   * The values printed are based on each process's average value.
   *
   * @param include_min_max whether to also print the min and max values
   */
  void PrintGlobalAverageNSeconds(bool include_min_max = false) const;

 private:
  /**
   * Returns the current time from the clock
   *
   * @return time as given by clock
   */
  timespec GetCurrentTime() const;

  timespec start_time_;  ///< time when the last timing began
  timespec total_time_;  ///< total time that has been measured
  timespec last_time_;  ///< time elasped in last cycle

  ///< Clock to use for timing
  static const clockid_t kClockType = CLOCK_MONOTONIC;

  ///< precision of nsec value. (number per second)
  static const time_t kPrecision = 100000000L;
};
#endif  // if !defined(__MACH__) && !defined(__APPLE__)


/**
 * Group object for timers.
 *
 * TODO: add example
 */
template<class TimerClass>
class TimerGroup {
 public:
  TimerGroup() : timers_() { }

  ~TimerGroup() { }

  /**
   * Adds a Timer to the group with the given label.
   *
   * The label is passed to the Time and is not used by this object. The labels
   * need not be unique. The returned value is intended to be a reference to
   * access the timer later and should likely be saved.
   *
   * @param label the label for the new timer.
   * @return an integer index of the timer that is being added
   */
  int AddTimer(std::string label) {
    timers_.push_back(TimerClass(label));
    return timers_.size() - 1;
  }

  /**
   * Returns a pointer to the Timer specified by the given index.
   *
   * The index values are the values that were returned when the timer was
   * added to the group. Index values outside of the range of existing timers
   * return NULL.
   *
   * @param index the index of the timer to return
   * @return timer specified by the index
   */
  TimerClass* GetTimer(unsigned int index) {
    // if out of range, return NULL
    if (index >= timers_.size()) {
      return NULL;
    }
    return &(timers_[index]);
  }

  /**
   * Start the timer specified by the given index.
   *
   * The index values are the values that were returned when the timer was
   * added to the group. Index values outside of the range of existing timers
   * do nothing.
   * TODO: should this throw an exception??
   *
   * @param index the index of the timer to start
   * @return timer specified by the index
   */
  void Start(unsigned int index) {
    // if out of range, do nothing
    if (index >= timers_.size()) {
      return;
    }
    timers_[index].Start();
  }

  /**
   * Stop the timer specified by the given index.
   *
   * The index values are the values that were returned when the timer was
   * added to the group. Index values outside of the range of existing timers
   * do nothing.
   * TODO: should this throw an exception??
   *
   * @param index the index of the timer to Stop
   * @return timer specified by the index
   */
  void Stop(unsigned int index) {
    // if out of range, do nothing
    if (index >= timers_.size()) {
      return;
    }
    timers_[index].Stop();
  }

  /**
   * Resets data from all timers.
   */
  void ResetAll() {
    typename std::vector<TimerClass>::iterator iter;
    for (iter = timers_.begin(); iter != timers_.end(); ++iter) {
      iter->Reset();
    }
  }

  /**
   * Prints the average time in seconds for each timers.
   */
  void PrintAvgSeconds() {
    typename std::vector<TimerClass>::iterator iter;
    for (iter = timers_.begin(); iter != timers_.end(); ++iter) {
      iter->PrintAvgSeconds();
    }
  }

  /**
   * Prints the total time in seconds for each timeres.
   * The total time is the sum of the time of all measurement cycles.
   */
  void PrintTotalTimeSeconds() {
    typename std::vector<TimerClass>::iterator iter;
    for (iter = timers_.begin(); iter != timers_.end(); ++iter) {
      iter->PrintTotalTimeSeconds();
    }
  }

 protected:
  std::vector<TimerClass> timers_;  ///< timers in group
};

#endif  // UTILS_TIMERS_H_
