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

// Tests the custom sphere-cylinder tests: distance, signed distance, and
// collision.

// TODO(SeanCurtis-TRI): Add tests for collision and distance; but, preferably,
// reformulate the code so that *one* function covers all three.

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "eigen_matrix_compare.h"
#include "fcl/geometry/shape/sphere.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_sphere-inl.h"
#include "fcl/narrowphase/collision.h"
#include "fcl/narrowphase/collision_object.h"
#include "fcl/narrowphase/collision_request.h"
#include "fcl/narrowphase/collision_result.h"
#include "fcl/narrowphase/distance.h"
#include "fcl/narrowphase/distance_request.h"
#include "fcl/narrowphase/distance_result.h"
#include "primitive_test_utility.h"

namespace fcl {
namespace detail {

// Specializations of the `GetCollide()` method so
// EvalCollisionForTestConfiguration() can pull up the right function.
template <>
CollideFunc<Sphere<double>, Sphere<double>> GetCollide() {
  return &sphereSphereIntersect<double>;
}

template <>
CollideFunc<Sphere<float>, Sphere<float>> GetCollide() {
  return &sphereSphereIntersect<float>;
}

template <>
DistanceFunc<Sphere<double>, Sphere<double>> GetDistance() {
  return &sphereSphereDistance<double>;
}

template <>
DistanceFunc<Sphere<float>, Sphere<float>> GetDistance() {
  return &sphereSphereDistance<float>;
}

template <>
SignedDistFunc<Sphere<double>, Sphere<double>> GetSignedDistance() {
  return &sphereSphereSignedDistance<double>;
}

template <>
SignedDistFunc<Sphere<float>, Sphere<float>> GetSignedDistance() {
  return &sphereSphereSignedDistance<float>;
}

namespace {

// Sphere-sphere analysis is very precise and the tests indicate that we only
// lose a few bits of precision due to arbitrary rotations of the relative
// position of spheres.
template <typename S>
struct Eps {
  using Real = typename constants<S>::Real;
  static Real value() { return 8 * constants<S>::eps(); }
};

template <typename S>
using Config = TestConfiguration<Sphere<S>, Sphere<S>>;

// Returns a collection of configurations of two spheres. The two spheres are
// designed to be *slightly* different sizes. The argument `scale_A` can be
// used to change sphere A's radius significantly -- *must* be greater than or
// equal to one.
template <typename S>
std::vector<Config<S>> GetTestConfigurations(const S& scale_A = S(1)) {
  std::vector<Config<S>> configurations;

  EXPECT_GE(scale_A, S(1));
  const S radiusA(0.75 * scale_A);
  const S radiusB(0.85);
  Vector3<S> dir = Vector3<S>{1, 2, 3}.normalized();
  Sphere<S> sphereA(radiusA);
  Sphere<S> sphereB(radiusB);
  using std::to_string;
  std::string prefix =
      scale_A == S(1)
          ? std::string("")
          : (std::string("A scaled to ") + to_string(scale_A) + " - ");

  {
    // Case: Completely separated in arbitrary direction.
    const S distance = 0.125;
    const Vector3<S> p_BA = (distance + radiusA + radiusB) * dir;
    configurations.emplace_back(prefix + "Separated in arbitrary direction", sphereA,
                                sphereB, Transform3<S>{Translation3<S>{p_BA}},
                                false);
    auto& config = configurations.back();
    // Not colliding --> no collision values.
    config.expected_distance = distance;
    config.p_BNa_expected = p_BA - (radiusA * dir);
    config.p_BNb_expected = dir * radiusB;
  }

  // The penetration point is not the *mid-point* between contact surfaces, but
  // the proportional to the relative radii of the spheres.
  // TODO(SeanCurtis-TRI): Determine if this is consistent with other collision
  // algorithms -- this is *very* worrying.
  auto compute_p_BP = [radiusA, radiusB, &dir](const S penetration_depth) {
    return dir * (radiusB - penetration_depth * radiusB / (radiusA + radiusB));
  };

  {
    // Case: Touching (interpreted as colliding).
    // Penetrate a bit less than machine precision and call that "contact".
    const S distance = -constants<S>::eps() * 0.5;
    const Vector3<S> p_BA = (distance + radiusA + radiusB) * dir;
    configurations.emplace_back("Touching in arbitrary direction", sphereA,
                                sphereB, Transform3<S>{Translation3<S>{p_BA}},
                                true);
    auto& config = configurations.back();
    // Colliding!
    config.expected_depth = -distance;
    config.n_AB_B_expected = -dir;
    config.p_BP_expected = compute_p_BP(-distance);
    // We also need p_BNa and p_BNb for signed distance.
    config.p_BNa_expected = p_BA - (radiusA * dir);
    config.p_BNb_expected = dir * radiusB;
  }

  {
    // Case: Simple penetration.
    using std::min;
    const S distance = -min(radiusA, radiusB) * 0.5;
    const Vector3<S> p_BA = (distance + radiusA + radiusB) * dir;
    configurations.emplace_back("Penetrating half the smaller radius", sphereA,
                                sphereB, Transform3<S>{Translation3<S>{p_BA}},
                                true);
    auto& config = configurations.back();
    // Colliding!
    config.expected_depth = -distance;
    config.n_AB_B_expected = -dir;
    config.p_BP_expected = compute_p_BP(-distance);
    // We also need p_BNa and p_BNb for signed distance.
    config.p_BNa_expected = p_BA - (radiusA * dir);
    config.p_BNb_expected = dir * radiusB;
  }

  {
    // Case: A arbitrarily rotated.
    const S distance = -0.1;
    const Vector3<S> p_BA = (distance + radiusA + radiusB) * dir;
    Transform3<S> X_BA(
        Translation3<S>{p_BA} *
        AngleAxis<S>(constants<S>::pi() / 4, Vector3<S>{1, 1, 2}.normalized()));
    configurations.emplace_back(
        "Penetration with non-zero relative orientation", sphereA, sphereB,
        X_BA, true);
    auto& config = configurations.back();
    // Colliding!
    config.expected_depth = -distance;
    config.n_AB_B_expected = -dir;
    config.p_BP_expected = compute_p_BP(-distance);
    // We also need p_BNa and p_BNb for signed distance.
    config.p_BNa_expected = p_BA - (radiusA * dir);
    config.p_BNb_expected = dir * radiusB;
  }

  {
    // Case: Coincident origins.
    using std::min;
    const S distance = -(radiusA + radiusB);
    const Vector3<S> p_BA{0, 0, 0};
    configurations.emplace_back("Coincident origins", sphereA,
                                sphereB, Transform3<S>{Translation3<S>{p_BA}},
                                true);
    auto& config = configurations.back();
    // Colliding!
    config.expected_depth = -distance;
    // TODO(SeanCurtis-TRI): This is *not* desirable behavior; an arbitrary
    // normal would be better.
    config.n_AB_B_expected << 0, 0, 0;
    config.p_BP_expected << 0, 0, 0;
    // We also need p_BNa and p_BNb for signed distance.
    config.p_BNa_expected << -radiusA, 0, 0;
    config.p_BNb_expected << radiusB, 0, 0;
  }

  {
    // Case: Spheres with radically (three orders of magnitude) different scales
    // colliding.
    Sphere<S> sphereA(radiusA * 1000);
    Sphere<S> sphereB(radiusB);
  }

  return configurations;
}

GTEST_TEST(SphereSpherePrimitiveTest, CollisionAcrossVaryingWorldFramesDouble) {
  for (const double scale_A : {1., 1000.}) {
    QueryWithVaryingWorldFrames<Config<double>>(
        GetTestConfigurations<double>(scale_A),
        &EvalCollisionForTestConfiguration<Config<double>>,
        Eps<double>::value());
  }
}

GTEST_TEST(SphereSpherePrimitiveTest, CollisionAcrossVaryingWorldFramesFloat) {
  for (const double scale_A : {1., 1000.}) {
    QueryWithVaryingWorldFrames<Config<float>>(
        GetTestConfigurations<float>(scale_A),
        &EvalCollisionForTestConfiguration<Config<float>>, Eps<float>::value());
  }
}

GTEST_TEST(SphereSpherePrimitiveTest, DistanceAcrossVaryingWorldFramesDouble) {
  for (const double scale_A : {1., 1000.}) {
    QueryWithVaryingWorldFrames<Config<double>>(
        GetTestConfigurations<double>(scale_A),
        &EvalDistanceForTestConfiguration<Config<double>>,
        Eps<double>::value());
  }
}

GTEST_TEST(SphereSpherePrimitiveTest, DistanceAcrossVaryingWorldFramesFloat) {
  for (const double scale_A : {1., 1000.}) {
    QueryWithVaryingWorldFrames<Config<float>>(
        GetTestConfigurations<float>(scale_A),
        &EvalDistanceForTestConfiguration<Config<float>>, Eps<float>::value());
  }
}

GTEST_TEST(SphereSpherePrimitiveTest,
           SignedDistanceAcrossVaryingWorldFramesDouble) {
  for (const double scale_A : {1., 1000.}) {
    QueryWithVaryingWorldFrames<Config<double>>(
        GetTestConfigurations<double>(scale_A),
        &EvalSignedDistanceForTestConfiguration<Config<double>>,
        Eps<double>::value());
  }
}

GTEST_TEST(SphereSpherePrimitiveTest,
           SignedDistanceAcrossVaryingWorldFramesFloat) {
  for (const double scale_A : {1., 1000.}) {
    QueryWithVaryingWorldFrames<Config<float>>(
        GetTestConfigurations<float>(scale_A),
        &EvalSignedDistanceForTestConfiguration<Config<float>>,
        Eps<float>::value());
  }
}

// This test confirms that the primitive test is exercised by the
// fcl::distance() method. We create to coincident spheres and rely on the fact
// that we'll only get the points aligned with the x-axis because of the
// primitive test.
GTEST_TEST(SphereSpherePrimitiveTest, CollidePrimitiveActiveLibccd) {
  const double expected_depth = 0.25;
  using CollisionGeometryPtr_t = std::shared_ptr<fcl::CollisionGeometryd>;
  CollisionGeometryPtr_t geo_A(new fcl::Sphered(0.5));
  CollisionObjectd obj_A(geo_A,
                         Transform3d{Translation3d{1 - expected_depth, 0, 0}});
  CollisionGeometryPtr_t geo_B(new fcl::Sphered(0.5));
  CollisionObjectd obj_B(geo_B, Transform3d::Identity());

  CollisionRequestd request;
  request.gjk_solver_type = GJKSolverType::GST_LIBCCD;
  // Huge tolerance to allow GJK to provide a really sloppy answer.
  request.gjk_tolerance = 1e-1;
  request.enable_contact = true;  // Return witnesses to penetration.
  CollisionResultd result;
  collide(&obj_A, &obj_B, request, result);
  // No epsilon; these tests should come out *perfectly*.
  ASSERT_TRUE(result.isCollision());
  EXPECT_EQ(result.numContacts(), 1u);
  const auto& contact = result.getContact(0);
  EXPECT_EQ(contact.penetration_depth, expected_depth);
  EXPECT_TRUE(CompareMatrices(contact.normal, Vector3d{-1, 0, 0}, 0,
                              MatrixCompareType::absolute));
  EXPECT_TRUE(CompareMatrices(contact.pos,
                              Vector3d{(1 - expected_depth) * 0.5, 0, 0}, 0,
                              MatrixCompareType::absolute));
}

// This test confirms that the primitive test is exercised by the
// fcl::distance() method. We create to coincident spheres and rely on the fact
// that we'll only get the points aligned with the x-axis because of the
// primitive test.
GTEST_TEST(SphereSpherePrimitiveTest, DistancePrimitiveActiveLibccd) {
  const double expected_dist = 0.25;
  using CollisionGeometryPtr_t = std::shared_ptr<fcl::CollisionGeometryd>;
  CollisionGeometryPtr_t geo_A(new fcl::Sphered(0.5));
  CollisionObjectd obj_A(geo_A,
                         Transform3d{Translation3d{1 + expected_dist, 0, 0}});
  CollisionGeometryPtr_t geo_B(new fcl::Sphered(0.5));
  CollisionObjectd obj_B(geo_B, Transform3d::Identity());

  DistanceRequestd request;
  request.gjk_solver_type = GJKSolverType::GST_LIBCCD;
  request.enable_nearest_points = true;
  request.enable_signed_distance = false;
  // Huge tolerance to allow GJK to provide a really sloppy answer.
  request.distance_tolerance = 1e-1;
  DistanceResultd result;
  double dist = distance(&obj_A, &obj_B, request, result);
  // No epsilon; these tests should come out *perfectly*.
  EXPECT_EQ(dist, expected_dist);
  EXPECT_TRUE(CompareMatrices(result.nearest_points[0],
                              Vector3d{0.5 + expected_dist, 0, 0}, 0,
                              MatrixCompareType::absolute));
  EXPECT_TRUE(CompareMatrices(result.nearest_points[1], Vector3d{0.5, 0, 0},
                              0, MatrixCompareType::absolute));
}

// This test confirms that the primitive test is exercised by the
// fcl::distance() method. We create to coincident spheres and rely on the fact
// that we'll only get the points aligned with the x-axis because of the
// primitive test.
GTEST_TEST(SphereSpherePrimitiveTest, SignedDistancePrimitiveActiveLibccd) {
  using CollisionGeometryPtr_t = std::shared_ptr<fcl::CollisionGeometryd>;
  CollisionGeometryPtr_t geo_A(new fcl::Sphered(0.5));
  CollisionObjectd obj_A(geo_A, Transform3d::Identity());
  CollisionGeometryPtr_t geo_B(new fcl::Sphered(0.5));
  CollisionObjectd obj_B(geo_B, Transform3d::Identity());

  DistanceRequestd request;
  request.gjk_solver_type = GJKSolverType::GST_LIBCCD;
  request.enable_nearest_points = true;
  request.enable_signed_distance = true;
  // Huge tolerance to allow GJK to provide a really sloppy answer.
  request.distance_tolerance = 1e-1;
  DistanceResultd result;
  double dist = distance(&obj_A, &obj_B, request, result);
  // No epsilon; these tests should come out *perfectly*.
  EXPECT_EQ(dist, -1.0);
  EXPECT_TRUE(CompareMatrices(result.nearest_points[0], Vector3d{-0.5, 0, 0},
                              0, MatrixCompareType::absolute));
  EXPECT_TRUE(CompareMatrices(result.nearest_points[1], Vector3d{0.5, 0, 0},
                              0, MatrixCompareType::absolute));
}

// TODO(SeanCurtis-TRI): Add pipeline for GJKSolverType::GST_INDEP.

}  // namespace
}  // namespace detail
}  // namespace fcl

//==============================================================================
int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
