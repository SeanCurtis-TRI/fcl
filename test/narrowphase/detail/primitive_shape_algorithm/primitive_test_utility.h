/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2019. Toyota Research Institute
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of CNRS-LAAS and AIST nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/** @author Sean Curtis (sean@tri.global) (2019) */

#ifndef FCL_PRIMITIVE_TEST_UTILITY_H
#define FCL_PRIMITIVE_TEST_UTILITY_H

#include <functional>
#include <type_traits>
#include <vector>

#include "eigen_matrix_compare.h"
#include "fcl/narrowphase/contact_point.h"
#include "fcl/math/constants.h"

namespace fcl {
namespace detail {

// Templated functions to acquire the primitive-primitive function to evaluate.
// In each of the primitive-primitive tests, these should be specialized to
// provide a pointer to the appropriate function.
template <typename ShapeA, typename ShapeB>
using CollideFunc = std::function<bool(
    const ShapeA& shapeA, const Transform3<typename ShapeA::S>& X_WA,
    const ShapeB& shapeB, const Transform3<typename ShapeA::S>& X_WB,
    std::vector<ContactPoint<typename ShapeA::S>>* contacts)>;

template <typename ShapeA, typename ShapeB>
CollideFunc<ShapeA, ShapeB> GetCollide() {
  throw std::logic_error(
      "No collision function defined; specialize this function for your "
      "shapes");
}

template <typename ShapeA, typename ShapeB>
using DistanceFunc = std::function<bool(
    const ShapeA& shapeA, const Transform3<typename ShapeA::S>& X_WA,
    const ShapeB& shapeB, const Transform3<typename ShapeA::S>& X_WB,
    typename ShapeA::S* dist, Vector3<typename ShapeA::S>* p1,
    Vector3<typename ShapeA::S>* p2)>;

template <typename ShapeA, typename ShapeB>
DistanceFunc<ShapeA, ShapeB> GetDistance() {
  throw std::logic_error(
      "No distance function defined; specialize this function for your "
      "shapes");
}

template <typename ShapeA, typename ShapeB>
using SignedDistFunc = std::function<typename ShapeA::S(
    const ShapeA& shapeA, const Transform3<typename ShapeA::S>& X_WA,
    const ShapeB& shapeB, const Transform3<typename ShapeA::S>& X_WB,
    Vector3<typename ShapeA::S>* p1, Vector3<typename ShapeA::S>* p2)>;

template <typename ShapeA, typename ShapeB>
SignedDistFunc<ShapeA, ShapeB> GetSignedDistance() {
  throw std::logic_error(
      "No signed distance function defined; specialize this function for your "
      "shapes");
}

// Defines the test configuration for a single test. It includes the geometry
// and the pose of the geometry A in the geometry B's (B). It also includes the
// expected answers in that same frame. It does not include those quantities
// that vary from test invocation to invocation (e.g., the pose of the B
// in the world frame or the *orientation* of the ).
//
// Collision and distance are complementary queries -- two objects in collision
// have no defined distance because they are *not* separated and vice versa.
// These configurations allow for the test of the complementarity property.
template <typename ShapeA_, typename ShapeB_>
struct TestConfiguration {
  static_assert(std::is_same<typename ShapeA_::S, typename ShapeB_::S>::value,
                "The shapes must use the same scalar type");
  using Scalar = typename ShapeA_::S;
  using ShapeA = ShapeA_;
  using ShapeB = ShapeB_;

  TestConfiguration(const std::string& name_in, const ShapeA_& shapeA_in,
                    const ShapeB_& shapeB_in, const Transform3<Scalar>& X_BA_in,
                    bool colliding)
      : name(name_in),
        shapeA(shapeA_in),
        shapeB(shapeB_in),
        X_BA(X_BA_in),
        expected_colliding(colliding) {}

  // Descriptive name of the test configuration.
  std::string name;
  // Shape A in the test.
  ShapeA shapeA;
  // Shape B in the test.
  ShapeB shapeB;
  // The pose of geometry A in the frame of geometry B.
  Transform3<Scalar> X_BA;

  // TODO(SeanCurtis-TRI): Figure out if this has value; should I provide a
  // different value and then derive this?

  // Indicates if this test configuration is expected to be in collision.
  bool expected_colliding{false};

  // TODO(SeanCurtis-TRI): How do I generalize these? I should provide methods
  // to compute these and then let individual tests specialize them. Or some
  // such thing.

  // Collision values; only valid if expected_colliding is true.
  Scalar expected_depth{-1};
  // The normal, pointing from A to B, expressed in the B frame.
  Vector3<Scalar> n_AB_B_expected;
  // The position vector from B's origin to the contact point P, expressed in B.
  Vector3<Scalar> p_BP_expected;

  // Distance values; only valid if expected_colliding is false.
  Scalar expected_distance{-1};
  // The points on geometry A and B, respectively, closest to the other shape
  // measured and expressed in frame B. These should always be defined
  // (regardless of whether the geometries are colliding or separated).
  Vector3<Scalar> p_BNa_expected;
  Vector3<Scalar> p_BNb_expected;
};


// Utility for creating a copy of the input configurations and appending more
// labels to the configuration name -- aids in debugging.
template <typename Configuration>
std::vector<Configuration> AppendLabel(
    const std::vector<Configuration>& configurations,
    const std::string& label) {
  std::vector<Configuration> configs;
  for (const auto& config : configurations) {
    configs.push_back(config);
    configs.back().name += " - " + label;
  }
  return configs;
}

// The callback type to evaluate on a test configuration given the
// test configuration, the pose of geometry B in the world, and a tolerance.
// TODO(SeanCurtis-TRI): Nail down these parameters in this more generalized
// framework.
template <typename Configuration>
using EvalFunc = std::function<void(const Configuration&,
                                    const Transform3<typename Configuration::Scalar>&,
                                    const typename Configuration::Scalar&)>;

// An Evaluation function that tests for *collision* on test configurations.
// This evaluates an instance of a test configuration and confirms the results
// match the expected data. The test configuration is defined in the geometry
// B's frame. The answers are framed in terms of this known relative pose.
// The X_WB value allows us to asses the error that creeps in when the relative
// pose is in some arbitrary configuration in the world frame.
//
// Evaluates the collision query twice. Once as the boolean "is colliding" test
// and once with the collision characterized with depth, normal, and position.
template <typename Configuration>
void EvalCollisionForTestConfiguration(
    const Configuration& config,
    const Transform3<typename Configuration::Scalar>& X_WB,
    const typename Configuration::Scalar& eps) {
  using S = typename Configuration::Scalar;
  using ShapeA = typename Configuration::ShapeA;
  using ShapeB = typename Configuration::ShapeB;

  auto collide = GetCollide<ShapeA, ShapeB>();

  // Collide with no witness data returned.
  const Transform3<S> X_WA = X_WB * config.X_BA;
  bool colliding = collide(config.shapeA, X_WA, config.shapeB, X_WB, nullptr);
  EXPECT_EQ(config.expected_colliding, colliding) << config.name;

  // Collide with witness data.
  std::vector<ContactPoint<S>> contacts;
  colliding = collide(config.shapeA, X_WA, config.shapeB, X_WB, &contacts);
  EXPECT_EQ(colliding, config.expected_colliding) << config.name;
  if (config.expected_colliding) {
    EXPECT_EQ(1u, contacts.size()) << config.name;
    const ContactPoint<S>& contact = contacts[0];
    EXPECT_NEAR(config.expected_depth, contact.penetration_depth, eps)
              << config.name;
    EXPECT_TRUE(CompareMatrices(contact.normal,
                                X_WB.linear() * config.n_AB_B_expected, eps,
                                MatrixCompareType::absolute))
              << config.name;
    EXPECT_TRUE(CompareMatrices(contact.pos, X_WB * config.p_BP_expected, eps,
                                MatrixCompareType::absolute))
              << config.name;
  } else {
    EXPECT_EQ(contacts.size(), 0u) << config.name;
  }
}

// An Evaluation function that tests for *separating distance* on test
// configurations.
// This evaluates an instance of a test configuration and confirms the results
// match the expected data. The test configuration is defined in the geometry
// B's frame. The answers are framed in terms of this known relative pose.
// The X_WB value allows us to asses the error that creeps in when the relative
// pose is in some arbitrary configuration in the world frame.
//
// Evaluates the distance query twice. Once as the boolean "is separated" test
// and once with the separation characterized with witness points.
template <typename Configuration>
void EvalDistanceForTestConfiguration(
    const Configuration& config,
    const Transform3<typename Configuration::Scalar>& X_WB,
    const typename Configuration::Scalar& eps) {
  using S = typename Configuration::Scalar;
  using ShapeA = typename Configuration::ShapeA;
  using ShapeB = typename Configuration::ShapeB;

  auto distance_func = GetDistance<ShapeA, ShapeB>();

  const Transform3<S> X_WA = X_WB * config.X_BA;
  bool separated = distance_func(config.shapeA, X_WA, config.shapeB, X_WB,
                                 nullptr, nullptr, nullptr);
  EXPECT_NE(separated, config.expected_colliding) << config.name;

  // Initializing this to -2, to confirm that a colliding scenario sets
  // distance to -1.
  // TODO(SeanCurtis-TRI): Confirm if this behavior is consistent and logical.
  // After all, one would assume the returned bool would serve this purpose.
  S distance{-2};
  Vector3<S> p_WNa{0, 0, 0};
  Vector3<S> p_WNb{0, 0, 0};

  separated = distance_func(config.shapeA, X_WA, config.shapeB, X_WB, &distance,
                            &p_WNa, &p_WNb);
  EXPECT_NE(separated, config.expected_colliding) << config.name;
  if (!config.expected_colliding) {
    EXPECT_NEAR(distance, config.expected_distance, eps) << config.name;
    EXPECT_TRUE(CompareMatrices(p_WNa, X_WB * config.p_BNa_expected, eps,
                                MatrixCompareType::absolute))
              << config.name;
    EXPECT_TRUE(CompareMatrices(p_WNb, X_WB * config.p_BNb_expected, eps,
                                MatrixCompareType::absolute))
              << config.name;
  } else {
    EXPECT_EQ(distance, S(-1)) << config.name;
    EXPECT_TRUE(CompareMatrices(p_WNa, Vector3<S>::Zero(), 0,
                                MatrixCompareType::absolute));
    EXPECT_TRUE(CompareMatrices(p_WNb, Vector3<S>::Zero(), 0,
                                MatrixCompareType::absolute));
  }
}

// An Evaluation function that tests for *signed distance* on test
// configurations.
// This evaluates an instance of a test configuration and confirms the results
// match the expected data. The test configuration is defined in the geometry
// B's frame. The answers are framed in terms of this known relative pose.
// The X_WB value allows us to asses the error that creeps in when the relative
// pose is in some arbitrary configuration in the world frame.
//
// Evaluates the distance query twice. Once to get just the signed distance and
// once to get the witness points.
template <typename Configuration>
void EvalSignedDistanceForTestConfiguration(
    const Configuration& config,
    const Transform3<typename Configuration::Scalar>& X_WB,
    const typename Configuration::Scalar& eps) {
  using S = typename Configuration::Scalar;
  using ShapeA = typename Configuration::ShapeA;
  using ShapeB = typename Configuration::ShapeB;

  auto distance_func = GetSignedDistance<ShapeA, ShapeB>();

  const Transform3<S> X_WA = X_WB * config.X_BA;
  S distance =
      distance_func(config.shapeA, X_WA, config.shapeB, X_WB, nullptr, nullptr);
  const S expected_distance = config.expected_colliding
                                  ? -config.expected_depth
                                  : config.expected_distance;
  EXPECT_NEAR(distance, expected_distance, eps)
      << "Simple query: " << config.name;

  Vector3<S> p_WNa{0, 0, 0};
  Vector3<S> p_WNb{0, 0, 0};
  distance =
      distance_func(config.shapeA, X_WA, config.shapeB, X_WB, &p_WNa, &p_WNb);

  EXPECT_NEAR(distance, expected_distance, eps)
      << "Query with witness: " << config.name;
  // The points p_BNa and p_BNb should be valid colliding or not.
  EXPECT_TRUE(CompareMatrices(p_WNa, X_WB * config.p_BNa_expected, eps,
                              MatrixCompareType::absolute))
            << config.name;
  EXPECT_TRUE(CompareMatrices(p_WNb, X_WB * config.p_BNb_expected, eps,
                              MatrixCompareType::absolute))
            << config.name;
}

// This test defines the transforms for performing the single collision test.
template <typename Configuration>
void QueryWithVaryingWorldFrames(
    const std::vector<Configuration>& configurations,
    EvalFunc<Configuration> query_eval,
    const typename Configuration::Scalar& eps) {
  using S = typename Configuration::Scalar;

  // Evaluate all the configurations with the given pose X_FB.
  auto evaluate_all = [&eps, query_eval](
      const std::vector<Configuration>& configs,
      const Transform3<S>& X_FB) {
    for (const auto config : configs) {
      query_eval(config, X_FB, eps);
    }
  };

  // Frame F is coincident and aligned with the geometry B's frame.
  Transform3<S> X_FB = Transform3<S>::Identity();
#if 0
  evaluate_all(AppendLabel(configurations, "X_FB = I"), X_FB);

  // Simple arbitrary translation away from the origin.
  X_FB.translation() << 1.3, 2.7, 6.5;
  evaluate_all(AppendLabel(configurations, "X_FB is translation"), X_FB);
#endif
  std::string axis_name[] = {"x", "y", "z"};
  // 90 degree rotation around each axis.
  for (int axis = 1; axis < 2; ++axis) {
    std::string label = "X_FB is 90-degree rotation around " + axis_name[axis];
    AngleAxis<S> angle_axis{constants<S>::pi() / 2, Vector3<S>::Unit(axis)};
    X_FB.linear() << angle_axis.matrix();
    evaluate_all(AppendLabel(configurations, label), X_FB);
  }
#if 0
  // Arbitrary orientation.
  {
    AngleAxis<S> angle_axis{constants<S>::pi() / 3,
                            Vector3<S>{1, 2, 3}.normalized()};
    X_FB.linear() << angle_axis.matrix();
    evaluate_all(AppendLabel(configurations, "X_FB is arbitrary rotation"),
                 X_FB);
  }

  // Near axis aligned.
  {
    AngleAxis<S> angle_axis{constants<S>::eps_12(), Vector3<S>::UnitX()};
    X_FB.linear() << angle_axis.matrix();
    evaluate_all(AppendLabel(configurations, "X_FB is near identity"),
                 X_FB);
  }
#endif
}

}  // namespace detail
}  // namespace fcl


#endif //FCL_PRIMITIVE_TEST_UTILITY_H
