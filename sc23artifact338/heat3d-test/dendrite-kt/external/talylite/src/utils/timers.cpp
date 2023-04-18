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
#include <talyfem/utils/timers.h>

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/mpi.h>
#else
#include <mpi.h>
#endif

#include <vector>

#include <talyfem/utils/utils.h>

using TALYFEMLIB::PrintTime;

void TimerBase::PrintAvgSeconds() const {
  PrintTime(label_, ": ", GetAvgTimeSeconds());
}

void TimerBase::PrintTotalTimeSeconds() const {
  PrintTime(label_, ": ", GetTotalTimeSeconds());
}

void TimerBase::PrintLastTimeSeconds() const {
  PrintTime(label_, ": ", GetLastTimeSeconds());
}

void TimerBase::PrintGlobalAverageSeconds(bool include_min_max) const {
  // set the values for this process
  double loc_time = GetAvgTimeSeconds();
  double avg_time = 0;
  double min_time = 0;
  double max_time = 0;

  MPI_Reduce(&loc_time, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (include_min_max) {
    MPI_Reduce(&loc_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&loc_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  }

  int mpi_size;  // used to calculate average
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  avg_time /= mpi_size;

  if (include_min_max) {
    PrintTime(label_, " (global_average_sec, min, max): ", avg_time, "  ",
              min_time, "  ", max_time);
  } else {
    PrintTime(label_, " (global_average_sec): ", avg_time);
  }
}

void TimerBase::PrintGlobalTotalSeconds(bool include_min_max) const {
  // set the values for this process
  double loc_time = GetTotalTimeSeconds();
  double avg_time = 0;
  double min_time = 0;
  double max_time = 0;

  MPI_Reduce(&loc_time, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (include_min_max) {
    MPI_Reduce(&loc_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&loc_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  }

  int mpi_size;  // used to calculate average
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  avg_time /= mpi_size;

  if (include_min_max) {
    PrintTime(label_, " (global_total_sec avg, min, max): ", avg_time, "  ",
              min_time, "  ", max_time);
  } else {
    PrintTime(label_, " (global_total_sec avg): ", avg_time);
  }
}

void TimerBase::PrintGlobalChart() const {
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  double loc_time = GetTotalTimeSeconds();

  std::vector<double> all_times;
  all_times.resize(mpi_size);
  MPI_Gather(&loc_time, 1, MPI_DOUBLE, all_times.data(), 1, MPI_DOUBLE, 0,
             PETSC_COMM_WORLD);
  TALYFEMLIB::PrintBarChart(all_times.data(), mpi_size, label_);
}


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void MPITimer::Start() {
  if (timer_running_)
    Stop();
  timer_running_ = true;
  start_time_ = GetCurrentTimeSeconds();
}

void MPITimer::ReStart() {
  timer_running_ = true;
  start_time_ = GetCurrentTimeSeconds();
}

void MPITimer::Stop() {
  double end_time = GetCurrentTimeSeconds();
  if (!timer_running_)
    return;
  timer_running_ = false;
  n_timings_++;
  last_time_ = (end_time - start_time_);
  total_time_ += last_time_;
}

void MPITimer::Reset() {
  timer_running_ = false;
  last_time_ = 0.0;
  start_time_ = 0.0;
  total_time_ = 0.0;
  n_timings_ = 0;
}

double MPITimer::GetAvgTimeSeconds() const {
  return total_time_ / n_timings_;
}

double MPITimer::GetTotalTimeSeconds() const {
  return total_time_;
}

double MPITimer::GetElapsedTimeSeconds() const {
  double current_time = GetCurrentTimeSeconds();
  if (!timer_running_)
    return 0.0;
  return current_time - start_time_;
}

double MPITimer::GetCurrentTimeSeconds() const {
  return MPI_Wtime();
}

double MPITimer::GetLastTimeSeconds() const {
  return last_time_;
}

void MPITimer::AddToTotalTime(double adjustment) {
  total_time_ += adjustment;
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#if !defined(__MACH__) && !defined(__APPLE__)
void HiResTimer::Start() {
  if (timer_running_)
    Stop();
  timer_running_ = true;
  start_time_ = GetCurrentTime();
}

void HiResTimer::ReStart() {
  timer_running_ = true;
  start_time_ = GetCurrentTime();
}

void HiResTimer::Stop() {
  timespec end_time = GetCurrentTime();
  if (!timer_running_)
    return;
  timer_running_ = false;
  n_timings_++;
  last_time_.tv_sec = (end_time.tv_sec - start_time_.tv_sec);
  last_time_.tv_nsec = (end_time.tv_nsec - start_time_.tv_nsec);
  total_time_.tv_sec += last_time_.tv_sec;
  total_time_.tv_nsec += last_time_.tv_nsec;
}

void HiResTimer::Reset() {
  timer_running_ = false;
  n_timings_ = 0;
  start_time_.tv_sec = 0.0;
  start_time_.tv_nsec = 0.0;
  total_time_.tv_sec = 0.0;
  total_time_.tv_nsec = 0.0;
  last_time_.tv_sec = 0.0;
  last_time_.tv_nsec = 0.0;
}

double HiResTimer::GetAvgTimeSeconds() const {
  time_t sec = total_time_.tv_sec;
  long nsec = total_time_.tv_nsec;  // NOLINT(runtime/int)
  double rval = sec + static_cast<double>(nsec) / kPrecision;
  rval /= n_timings_;
  return rval;
}

long HiResTimer::GetAvgTimeNSeconds() const {  // NOLINT(runtime/int)
  time_t sec = total_time_.tv_sec;
  long nsec = total_time_.tv_nsec;  // NOLINT(runtime/int)
  double rval = kPrecision * sec + nsec;
  rval /= n_timings_;
  return rval;
}

double HiResTimer::GetTotalTimeSeconds() const {
  time_t sec = total_time_.tv_sec;
  long nsec = total_time_.tv_nsec;  // NOLINT(runtime/int)
  double rval = sec + static_cast<double>(nsec) / kPrecision;
  return rval;
}

double HiResTimer::GetElapsedTimeSeconds() const {
  timespec time_now = GetCurrentTime();
  time_t sec = time_now.tv_sec - start_time_.tv_sec;
  long nsec = time_now.tv_nsec - start_time_.tv_nsec;  // NOLINT(runtime/int)
  double rval = sec + static_cast<double>(nsec) / kPrecision;
  return rval;
}

timespec HiResTimer::GetLastTime() const {
  return last_time_;
}

double HiResTimer::GetLastTimeSeconds() const {
  time_t sec = last_time_.tv_sec;
  long nsec = last_time_.tv_nsec;  // NOLINT(runtime/int)
  double rval = sec + static_cast<double>(nsec) / kPrecision;
  return rval;
}

void HiResTimer::AddToTotalTime(double adjustment) {
  // split double into seconds and nanoseconds
  int sec = static_cast<int>(adjustment);
  int nsec = kPrecision * (adjustment - sec);
  total_time_.tv_sec += static_cast<time_t>(sec);
  total_time_.tv_nsec += static_cast<time_t>(nsec);
}

void HiResTimer::PrintGlobalAverageNSeconds(bool include_min_max) const {
  // set the values for this process
  long loc_time = GetAvgTimeNSeconds();  // NOLINT(runtime/int)
  long avg_time = 0;  // NOLINT(runtime/int)
  long min_time = 0;  // NOLINT(runtime/int)
  long max_time = 0;  // NOLINT(runtime/int)

  MPI_Reduce(&loc_time, &avg_time, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if (include_min_max) {
    MPI_Reduce(&loc_time, &min_time, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&loc_time, &max_time, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
  }

  int mpi_size;  // used to calculate average
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  avg_time /= mpi_size;

  if (include_min_max) {
    PrintTime(label_, " (global_average_nsecs, min, max): ", avg_time, "  ",
              min_time, "  ", max_time);
  } else {
    PrintTime(label_, " (global_average_nsecs): ", avg_time);
  }
}

timespec HiResTimer::GetCurrentTime() const {
  timespec rval;
  clock_gettime(kClockType, &rval);
  return rval;
}

#endif  // if !defined(__MACH__) && !defined(__APPLE__)
