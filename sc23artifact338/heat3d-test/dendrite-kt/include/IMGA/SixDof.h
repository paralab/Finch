#pragma once

#include "IMGA/MathOperation.h"
#include "talyfem/talyfem.h"

using namespace MathOp;
using namespace std;
using namespace TALYFEMLIB;

class SixDof {
 public:

  /// Kinematic property
  struct KinematicProperty {
    double mass = 0.0;
    std::vector<ZEROPTV> MOI;
    std::vector<ZEROPTV> InverseMOI;
  };

  struct ForceInfo {
    ZEROPTV Force;
    ZEROPTV Torque;
  };

  /// Kinematic state, could be changing
  struct KinematicState {
    // center of mass
    ZEROPTV x;

    // rotation position (body fitted axis compared to global)
    std::vector<ZEROPTV> R;

    // linear momentum
    ZEROPTV P;

    // rotational momentum
    ZEROPTV L;

    // body moment inertia (transformed)
    std::vector<ZEROPTV> InverseMOIbody;

    // translational velocity
    ZEROPTV v;
    // rotation velocity
    ZEROPTV omega;

    // linear displacement to update since x is a global one
    ZEROPTV LinearDisplacementUpdate;
    // rotation matrix used to update the location since R is a global one.
    std::vector<ZEROPTV> RotationMatrixUpdate;
    KinematicState() {
      R.resize(3);
      for (int dim = 0; dim < 3; dim++) {
        R[dim](dim) = 1.0;
      }
      RotationMatrixUpdate.resize(3);
      InverseMOIbody.resize(3);
    };
  };

  KinematicProperty kp_;
  // Maximum of 7 to support bdf6
  int validKinamatic = 0;
  std::vector<KinematicState> ks_{7};
  std::vector<ForceInfo> fi_{7};

  SixDof() {
    kp_.MOI.resize(3);
    kp_.InverseMOI.resize(3);
  }

  ~SixDof() = default;

  /// ode45 solver (should be same as backwardEuler, haven't verified yet)
  void ode45(const ZEROPTV &force, const ZEROPTV &torque, const double dt) {
    throw TALYFEMLIB::TALYException() << "shouldn't use ode45";
    // does not support time yet, basically force and torque is constant and not dependent on the time, this is OK
    // because we are inside a timestep and we do not store the force and torque over time. Surely not estimate on the
    // forward time
    std::vector<double> array_in = StateToArray(ks_.back());
    std::vector<double> k1 = dydt(array_in, force, torque);
    std::vector<double> k2 = dydt(stdVectorAdd(array_in, stdVectorMultiply(k1, 0.5 * dt)), force, torque);
    std::vector<double> k3 = dydt(stdVectorAdd(array_in, stdVectorMultiply(k2, 0.5 * dt)), force, torque);
    std::vector<double> k4 = dydt(stdVectorAdd(array_in, stdVectorMultiply(k3, dt)), force, torque);

    // sum of all the arrays
    std::vector<double> k2_prime = stdVectorMultiply(k2, 2.0);
    std::vector<double> k3_prime = stdVectorMultiply(k3, 2.0);
    std::vector<double> temp_1 = stdVectorAdd(k1, k4);
    std::vector<double> temp_2 = stdVectorAdd(temp_1, k2_prime);
    std::vector<double> temp_3 = stdVectorAdd(temp_2, k3_prime);
    std::vector<double> array_out = stdVectorAdd(array_in, stdVectorMultiply(temp_3, (dt / 6.0)));

    // use the final state to calculate the new state
    KinematicState ks_out{};
    ArrayToState(array_out, ks_out);
    // update state, calculate the displacement along the way
    setArrayToState(ks_out);

  }

  /// set ks_in while finished, each of this function call set up a
  void forwardEuler(const ZEROPTV &force, const ZEROPTV &torque, const double dt) {
    // simple forward euler, validated with matlab code
    std::vector<double> array_in = StateToArray(ks_.back());  // y_n
    std::vector<double> k1 = dydt(array_in, force, torque); // f(t_n1, y_n1)

    // sum of all the arrays
    std::vector<double> array_out = stdVectorAdd(array_in, stdVectorMultiply(k1, dt));

    // use the final state to calculate the new state
    KinematicState ks_out{};
    ArrayToState(array_out, ks_out);
    // update state, calculate the displacement along the way
    setArrayToState(ks_out);

  }

  /// backwardEuler (BDF1)
  void backwardEuler(const ZEROPTV &force, const ZEROPTV &torque, const double dt) {
    assert(validKinamatic >= 1);
    // use a fixed point iteration
    std::vector<double> array_in = StateToArray(ks_.back());  // y_n
    // guess
    std::vector<double> y_temp = stdVectorAdd(array_in, stdVectorMultiply(dydt(array_in, force, torque), dt));

    double error = 1.0;
    const double tol = 1e-10;
    while (error > tol) {
      std::vector<double> array_out = stdVectorAdd(array_in, stdVectorMultiply(dydt(y_temp, force, torque), dt));
      error = 0.0;
      for (int i = 0; i < 18; i++) {
        error += (y_temp[i] - array_out[i]) * (y_temp[i] - array_out[i]);
      }
      error = sqrt(error);
      y_temp = array_out;
    }

    // use the final state to calculate the new state
    KinematicState ks_out{};
    ArrayToState(y_temp, ks_out);
    // update state, calculate the displacement along the way
    setArrayToState(ks_out);

  }

  /// BDF2, fall back to BDF1 (backward euler) when there is not enough data
  void BDF2(const ZEROPTV &force, const ZEROPTV &torque, const double dt) {
    if (validKinamatic <= 1) {
      backwardEuler(force, torque, dt);
      return;
    }
    assert(validKinamatic >= 2);
    // use a fixed point iteration
    std::vector<double> y_n = StateToArray(ks_.rbegin()[1]);  // y_n
    std::vector<double> y_n2_term = stdVectorMultiply(StateToArray(ks_.rbegin()[2]), -1.0 / 3.0);
    // guess
    std::vector<double> y_temp = stdVectorAdd(y_n, stdVectorMultiply(dydt(y_n, force, torque), dt));

    double error = 1.0;
    const double tol = 1e-10;
    while (error > tol) {
      std::vector<double> t1 = stdVectorMultiply(dydt(y_temp, force, torque), dt * 2.0 / 3.0);
      std::vector<double> t2 = stdVectorAdd(stdVectorMultiply(y_n, 4.0 / 3.0), t1);
      std::vector<double> array_out = stdVectorAdd(t2, y_n2_term);

      error = 0.0;
      for (int i = 0; i < 18; i++) {
        error += (y_temp[i] - array_out[i]) * (y_temp[i] - array_out[i]);
      }
      error = sqrt(error);
      y_temp = array_out;
    }

    // use the final state to calculate the new state
    KinematicState ks_out{};
    ArrayToState(y_temp, ks_out);
    // update state, calculate the displacement along the way
    setArrayToState(ks_out);
  }

  /// BDF3, fall back to BDF2 when there is not enough data
  void BDF3(const ZEROPTV &force, const ZEROPTV &torque, const double dt) {
    if (validKinamatic <= 2) {
      BDF2(force, torque, dt);
      return;
    }
    assert(validKinamatic >= 3);
    // use a fixed point iteration
    std::vector<double> y_n = StateToArray(ks_.rbegin()[1]);  // y_n
    std::vector<double> y_n2_term = stdVectorMultiply(StateToArray(ks_.rbegin()[2]), -9.0 / 11.0);
    std::vector<double> y_n3_term = stdVectorMultiply(StateToArray(ks_.rbegin()[3]), 2.0 / 11.0);
    // guess
    std::vector<double> y_temp = stdVectorAdd(y_n, stdVectorMultiply(dydt(y_n, force, torque), dt));

    double error = 1.0;
    const double tol = 1e-10;
    while (error > tol) {
      std::vector<double> t1 = stdVectorMultiply(dydt(y_temp, force, torque), dt * 6.0 / 11.0);
      std::vector<double> t2 = stdVectorAdd(stdVectorMultiply(y_n, 18.0 / 11.0), t1);
      std::vector<double> t3 = stdVectorAdd(t2, y_n2_term);
      std::vector<double> array_out = stdVectorAdd(t3, y_n3_term);
      error = 0.0;
      for (int i = 0; i < 18; i++) {
        error += (y_temp[i] - array_out[i]) * (y_temp[i] - array_out[i]);
      }
      error = sqrt(error);
      y_temp = array_out;
    }

    // use the final state to calculate the new state
    KinematicState ks_out{};
    ArrayToState(y_temp, ks_out);
    // update state, calculate the displacement along the way
    setArrayToState(ks_out);
  }

  /// BDF4, fall back to BDF3 when there is not enough data
  void BDF4(const ZEROPTV &force, const ZEROPTV &torque, const double dt) {
    if (validKinamatic <= 3) {
      BDF3(force, torque, dt);
      return;
    }
    assert(validKinamatic >= 4);
    // use a fixed point iteration
    std::vector<double> y_n = StateToArray(ks_.rbegin()[1]);  // y_n
    std::vector<double> y_n2_term = stdVectorMultiply(StateToArray(ks_.rbegin()[2]), -36.0 / 25.0);
    std::vector<double> y_n3_term = stdVectorMultiply(StateToArray(ks_.rbegin()[3]), 16.0 / 25.0);
    std::vector<double> y_n4_term = stdVectorMultiply(StateToArray(ks_.rbegin()[4]), -3.0 / 25.0);
    // guess
    std::vector<double> y_temp = stdVectorAdd(y_n, stdVectorMultiply(dydt(y_n, force, torque), dt));

    double error = 1.0;
    const double tol = 1e-10;
    while (error > tol) {
      std::vector<double> t1 = stdVectorMultiply(dydt(y_temp, force, torque), dt * 12.0 / 25.0);
      std::vector<double> t2 = stdVectorAdd(stdVectorMultiply(y_n, 48.0 / 25.0), t1);
      std::vector<double> t3 = stdVectorAdd(t2, y_n2_term);
      std::vector<double> t4 = stdVectorAdd(t3, y_n3_term);
      std::vector<double> array_out = stdVectorAdd(t4, y_n4_term);
      error = 0.0;
      for (int i = 0; i < 18; i++) {
        error += (y_temp[i] - array_out[i]) * (y_temp[i] - array_out[i]);
      }
      error = sqrt(error);
      y_temp = array_out;
    }

    // use the final state to calculate the new state
    KinematicState ks_out{};
    ArrayToState(y_temp, ks_out);
    // update state, calculate the displacement along the way
    setArrayToState(ks_out);
  }

  /// BDF5, fall back to BDF4 when there is not enough data
  void BDF5(const ZEROPTV &force, const ZEROPTV &torque, const double dt) {
    if (validKinamatic <= 4) {
      BDF4(force, torque, dt);
      return;
    }
    assert(validKinamatic >= 5);
    // use a fixed point iteration
    std::vector<double> y_n = StateToArray(ks_.rbegin()[1]);  // y_n
    std::vector<double> y_n2_term = stdVectorMultiply(StateToArray(ks_.rbegin()[2]), -300.0 / 137.0);
    std::vector<double> y_n3_term = stdVectorMultiply(StateToArray(ks_.rbegin()[3]), 200.0 / 137.0);
    std::vector<double> y_n4_term = stdVectorMultiply(StateToArray(ks_.rbegin()[4]), -75.0 / 137.0);
    std::vector<double> y_n5_term = stdVectorMultiply(StateToArray(ks_.rbegin()[5]), 12.0 / 137.0);
    // guess
    std::vector<double> y_temp = stdVectorAdd(y_n, stdVectorMultiply(dydt(y_n, force, torque), dt));

    double error = 1.0;
    const double tol = 1e-10;
    while (error > tol) {
      std::vector<double> t1 = stdVectorMultiply(dydt(y_temp, force, torque), dt * 60.0 / 137.0);
      std::vector<double> t2 = stdVectorAdd(stdVectorMultiply(y_n, 300.0 / 137.0), t1);
      std::vector<double> t3 = stdVectorAdd(t2, y_n2_term);
      std::vector<double> t4 = stdVectorAdd(t3, y_n3_term);
      std::vector<double> t5 = stdVectorAdd(t4, y_n4_term);
      std::vector<double> array_out = stdVectorAdd(t5, y_n5_term);
      error = 0.0;
      for (int i = 0; i < 18; i++) {
        error += (y_temp[i] - array_out[i]) * (y_temp[i] - array_out[i]);
      }
      error = sqrt(error);
      y_temp = array_out;
    }

    // use the final state to calculate the new state
    KinematicState ks_out{};
    ArrayToState(y_temp, ks_out);
    // update state, calculate the displacement along the way
    setArrayToState(ks_out);
  }

  /// BDF6, fall back to BDF5 when there is not enough data
  void BDF6(const ZEROPTV &force, const ZEROPTV &torque, const double dt) {
#if (DIM==2)
    // force.z() should be zero
    // torque.x() nad torque.y() should be zero
    ZEROPTV force_mod = force;
    ZEROPTV torque_mod = torque;
    force_mod.z() = 0.0;
    torque_mod.x() = 0.0;
    torque_mod.y() = 0.0;
#endif
    if (validKinamatic <= 5) {
      BDF5(force, torque, dt);
      return;
    }
    assert(validKinamatic >= 6);
    // use a fixed point iteration
    std::vector<double> y_n = StateToArray(ks_.rbegin()[1]);  // y_n
    std::vector<double> y_n2_term = stdVectorMultiply(StateToArray(ks_.rbegin()[2]), -450.0 / 147.0);
    std::vector<double> y_n3_term = stdVectorMultiply(StateToArray(ks_.rbegin()[3]), 400.0 / 147.0);
    std::vector<double> y_n4_term = stdVectorMultiply(StateToArray(ks_.rbegin()[4]), -225.0 / 147.0);
    std::vector<double> y_n5_term = stdVectorMultiply(StateToArray(ks_.rbegin()[5]), 72.0 / 147.0);
    std::vector<double> y_n6_term = stdVectorMultiply(StateToArray(ks_.rbegin()[6]), -10.0 / 147.0);
    // guess
    std::vector<double> y_temp = stdVectorAdd(y_n, stdVectorMultiply(dydt(y_n, force, torque), dt));

    double error = 1.0;
    const double tol = 1e-10;
    while (error > tol) {
      std::vector<double> t1 = stdVectorMultiply(dydt(y_temp, force, torque), dt * 60.0 / 147.0);
      std::vector<double> t2 = stdVectorAdd(stdVectorMultiply(y_n, 360.0 / 147.0), t1);
      std::vector<double> t3 = stdVectorAdd(t2, y_n2_term);
      std::vector<double> t4 = stdVectorAdd(t3, y_n3_term);
      std::vector<double> t5 = stdVectorAdd(t4, y_n4_term);
      std::vector<double> t6 = stdVectorAdd(t5, y_n5_term);
      std::vector<double> array_out = stdVectorAdd(t6, y_n6_term);
      error = 0.0;
      for (int i = 0; i < 18; i++) {
        error += (y_temp[i] - array_out[i]) * (y_temp[i] - array_out[i]);
      }
      error = sqrt(error);
      y_temp = array_out;
    }

    // use the final state to calculate the new state
    KinematicState ks_out{};
    ArrayToState(y_temp, ks_out);
    // update state, calculate the displacement along the way
    setArrayToState(ks_out);
  }

  /**
   * store the previous fi_ and ks_
   */
  void storePreviousInfo() {
    for (int i = 0; i < (fi_.size() - 1); i++) {
      fi_[i] = fi_[i + 1];
      ks_[i] = ks_[i + 1];
    }
  }

  void PrintStateDebug() {
    PrintStatus("-------------- 6DOF debug -----------------");
    PrintStatus("force: ", fi_.back().Force);
    PrintStatus("torque: ", fi_.back().Torque);
    PrintStatus("x: ", ks_.back().x);
    PrintStatus("R");
    for (int i = 0; i < 3; i++) {
      PrintStatus(ks_.back().R[i](0), "\t", ks_.back().R[i](1), "\t", ks_.back().R[i](2));
    }
    PrintStatus("P: ", ks_.back().P);
    PrintStatus("L: ", ks_.back().L);
    PrintStatus("InverseMOIbody");
    for (int i = 0; i < 3; i++) {
      PrintStatus(ks_.back().InverseMOIbody[i](0),
                  "\t",
                  ks_.back().InverseMOIbody[i](1),
                  "\t",
                  ks_.back().InverseMOIbody[i](2));
    }
    PrintStatus("v: ", ks_.back().v);
    PrintStatus("omega: ", ks_.back().omega);

    PrintStatus("mass: ", kp_.mass);
    PrintStatus("--------- MOI: --------- ");
    for (int i = 0; i < 3; i++) {
      PrintStatus(kp_.MOI[i](0), "\t", kp_.MOI[i](1), "\t", kp_.MOI[i](2));
    }
    PrintStatus("--------- InverseMOI: --------- ");
    for (int i = 0; i < 3; i++) {
      PrintStatus(kp_.InverseMOI[i](0), "\t", kp_.InverseMOI[i](1), "\t", kp_.InverseMOI[i](2));
    }
    PrintStatus("----------Rotation_Matrix_update: --------- ");
    for (int i = 0; i < 3; i++) {
      PrintStatus(ks_.back().RotationMatrixUpdate[i](0),
                  "\t",
                  ks_.back().RotationMatrixUpdate[i](1),
                  "\t",
                  ks_.back().RotationMatrixUpdate[i](2));
    }
    PrintStatus("displacement update: ", ks_.back().LinearDisplacementUpdate);
    PrintStatus("-------------- 6DOF debug -----------------");
  };



  /// The following function is hidden from the user as their combination should be called.
  /// Theses function should not be called alone
 private:

  /**
   * transfrom the kinematic state into an array of 18 by 1, this will be used for the ode as y.
   * @param ks_in
   * @return
   */
  std::vector<double> StateToArray(const KinematicState &ks_in) {
    // set up a array of length 18 to hold the data
    std::vector<double> array(18);
    for (int i = 0; i < 3; i++) {
      array[i] = ks_in.x(i);
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        array[i * 3 + j + 3] = ks_in.R[i](j);
      }
    }
    for (int i = 0; i < 3; i++) {
      array[i + 12] = ks_in.P(i);
    }
    for (int i = 0; i < 3; i++) {
      array[i + 15] = ks_in.L(i);
    }
    return array;
  }

  /**
   * transfrom an array into a kinematic state struct of 18 by 1
   * @param array_in
   * @param ks_out
   */
  void ArrayToState(const std::vector<double> &array_in, KinematicState &ks_out) {
    assert(array_in.size() == 18);
    for (int i = 0; i < 3; i++) {
      ks_out.x(i) = array_in[i];
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        ks_out.R[i](j) = array_in[i * 3 + j + 3];
      }
    }
    for (int i = 0; i < 3; i++) {
      ks_out.P(i) = array_in[i + 12];
    }
    for (int i = 0; i < 3; i++) {
      ks_out.L(i) = array_in[i + 15];
    }

    ks_out.v = ks_out.P * (1.0 / kp_.mass);

    std::vector<ZEROPTV> R_transpose = matTranspose(ks_out.R);

    ks_out.InverseMOIbody = mat_mat_multi(ks_out.R, mat_mat_multi(kp_.InverseMOI, R_transpose));

    ks_out.omega = mat_vec_multi(ks_out.InverseMOIbody, ks_out.L);

  }

  /**
   * Update the state in ks_, use the difference between the previous state and the new state to calculate the displacement
   * @param state_in
   */
  void setArrayToState(KinematicState &state_in) {
    ZEROPTV LinearDisplacementUpdate = state_in.x - this->ks_.back().x;
    std::vector<ZEROPTV> RotationMatrixUpdate = mat_mat_multi(state_in.R, inverse_matrix(this->ks_.back().R));
    this->ks_.back() = state_in;
    this->ks_.back().LinearDisplacementUpdate = LinearDisplacementUpdate;
    this->ks_.back().RotationMatrixUpdate = RotationMatrixUpdate;
  }

  /**
   * Given an array and force, return the time derivative of the array.
   * @param y_in
   * @param force
   * @param torque
   * @return
   */
  std::vector<double> dydt(const std::vector<double> &y_in, const ZEROPTV &force, const ZEROPTV &torque) {
    /// set v, InverseMOIBody, and omega
    KinematicState ks_out{};
    ArrayToState(y_in, ks_out);
    // not ideal, try harder
    setForceInfo(force, torque);
    std::vector<double> y_dot = ddtStateToArray(ks_out);
    return y_dot;
  }

  /**
   * Used only in dydt when the force and torque is updated, this actually calculates the y_dot
   * @param ks_in
   * @return
   */
  std::vector<double> ddtStateToArray(const KinematicState &ks_in) {
    std::vector<double> y_dot(18);
    for (int i = 0; i < 3; i++) {
      y_dot[i] = ks_in.v(i);
    }
    std::vector<ZEROPTV> Rdot;
    Rdot.resize(3);
    Rdot = mat_mat_multi(Star(ks_in.omega), ks_in.R);

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        y_dot[i * 3 + j + 3] = Rdot[i](j);
      }
    }

    for (int i = 0; i < 3; i++) {
      y_dot[i + 12] = fi_.back().Force(i);
    }
    for (int i = 0; i < 3; i++) {
      y_dot[i + 15] = fi_.back().Torque(i);
    }
    return y_dot;
  }

  /**
   * multiply a vector with a constant
   * @tparam T
   * @param vec_in
   * @param multiply
   * @return
   */
  template<typename T>
  std::vector<T> stdVectorMultiply(const std::vector<T> &vec_in, const T multiply) {
    std::vector<T> vec_out(vec_in);
    std::transform(vec_out.begin(), vec_out.end(), vec_out.begin(),
                   bind2nd(std::multiplies<T>(), multiply));
    return vec_out;
  }

  /**
   * return the sum of two vectors
   * @tparam T
   * @param vec_in_1
   * @param vec_in_2
   * @return
   */
  template<typename T>
  std::vector<T> stdVectorAdd(const std::vector<T> &vec_in_1, const std::vector<T> &vec_in_2) {
    assert(vec_in_1.size() == vec_in_2.size());
    std::vector<T> vec_out;
    vec_out.reserve(vec_in_1.size());
    std::transform(vec_in_1.begin(), vec_in_1.end(), vec_in_2.begin(), std::back_inserter(vec_out), std::plus<T>());
    return vec_out;
  }

  void setForceInfo(const ZEROPTV &force, const ZEROPTV &torque) {
    fi_.back().Force = force;
    fi_.back().Torque = torque;
  }

};
