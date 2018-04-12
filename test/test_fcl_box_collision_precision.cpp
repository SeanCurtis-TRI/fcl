/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2018. Toyota Research Institute
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

/** @author Sean Curtis */
#include <memory>
#include <sstream>
#include <utility>

#include <gtest/gtest.h>

#include "eigen_matrix_compare.h"
#include "fcl/narrowphase/collision.h"
#include "fcl/narrowphase/collision_object.h"

using namespace fcl;
template <typename S>
using CollisionGeometryPtr_t = std::shared_ptr<fcl::CollisionGeometry<S>>;

// This performs a very specific test. It collides an oriented box with a
// surface that is tangent to the z = 0 plane. The box is a cube with unit size.
// It is oriented and positioned so that a single corner, placed on the z axis,
// is the point on the cube that most deeply penetrates the tangent plane.
//
// For a given axis, we should have *this* picture
//        ┆  ╱╲
//        ┆ ╱  ╲
//        ┆╱    ╲
//       ╱┆╲    ╱
//      ╱ ┆ ╲  ╱
//     ╱  ┆  ╲╱
//     ╲  ┆  ╱
//      ╲ ┆ ╱
//  _____╲┆╱_____    ╱____ With small penetration depth of d
//        ┆          ╲
//        ┆
//        ┆
//        ┆
//
// We can use this against various objects to determine uniformity of behavior.
// As long as the object osculates the z = 0 plane at (0, 0, 0), then we should
// get a fixed, known contact (assuming the tangent object is A and the box is
// B):
//   Normal is (0, 0, 1) from plane into the box
//   penetration depth is the specified size
//   contact position is (0, 0, -depth / 2.
//
// More particularly, these tests are parameterized in the following way:
//
//  1. We want to confirm that the precision with which we get the targeted
//     answer scales with the scalar.
//  2. We want to make sure that we get the same precision across all tangent
//     geometries.
//
// This test has been introduced because it is known that we *don't* get the
// same answer across primitives.  See below for details.
//
// NOTE: This *isn't* a ::testing::Test class because we need a bit more control
// over the implementations based on scalar type.
template <typename S>
class BoxPrecisionTest {
 public:
  // Constructs and initializes the box in its penultimate test configuration.
  // The actual test will move it along the z-axis to determine the desired
  // penetration depth.
  BoxPrecisionTest(GJKSolverType solver_type) : solver_type_(solver_type) {
    const S pi = constants<S>::pi();
    // Set up box; this performs the orientation, determines the position of the
    // corner to apply the appropriate offset and then confirms the effect.
    Transform3<S> box_pose{AngleAxis<S>(-pi / 4, Vector3<S>::UnitY()) *
                           AngleAxis<S>(pi / 4, Vector3<S>::UnitX())};
    Vector3<S> corner{-1. / 2, -1. / 2, -1. / 2};
    Vector3<S> rotated_corner = box_pose.linear() * corner;
    // This will require one last bump in the (0, 0, -depth) direction once
    // target penetration depth is specified.
    box_offset_ = -rotated_corner;
    box_pose.translation() = box_offset_;

    // Confirm the transform moves the corner (-1, -1, -1) corner to the origin.
    Vector3<S> test_corner = box_pose * corner;
    EXPECT_NEAR(test_corner(0), 0, kDefaultTolerance);
    EXPECT_NEAR(test_corner(1), 0, kDefaultTolerance);
    EXPECT_NEAR(test_corner(2), 0, kDefaultTolerance);

    CollisionGeometryPtr_t<S> box_geometry(new Box<S>(1, 1, 1));
    box_.reset(new CollisionObject<S>(box_geometry, box_pose));
    InitializeTangentSphere();
    InitializeTangentPlane();
    InitializeTangentBox();
  }

  CollisionObject<S>* tangent_sphere() { return tangent_sphere_.get(); }
  CollisionObject<S>* tangent_plane() { return tangent_plane_.get(); }
  CollisionObject<S>* tangent_box() { return tangent_box_.get(); }

  // Given a tangent object and a desired penetration, confirms that a collision
  // is (or is not) reported. The bool `expect_collision` reports whether a
  // collision is expected or not.
  void AssertCollision(CollisionObject<S>* tangent_object, S depth,
                       bool expect_collision, const char* source) {
    CollisionResult<S> collisionResult;
    PerformCollision(tangent_object, depth, &collisionResult);
    if (expect_collision) {
      ASSERT_TRUE(collisionResult.isCollision())
          << source << " at depth " << depth << " failed to report a collision";
    } else {
      ASSERT_FALSE(collisionResult.isCollision())
          << source << " at depth " << depth
          << " actually reported a collision";
    }
  }

  // This test asserts that the tests fail with the given threshold -- this
  // should be the *biggest* threshold value for which this is true.
  void AssertInsufficientPrecision(CollisionObject<S>* tangent_object, S depth,
                                   S test_tolerance) {
    CollisionResult<S> collisionResult;
    PerformCollision(tangent_object, depth, &collisionResult);

    // If there are no collisions, or the wrong number, then something more
    // fundamental went wrong.
    if (!collisionResult.isCollision()) {
      GTEST_FAIL() << "Fundamental error! No collisions reported";
    }

    std::vector<Contact<S>> contacts;
    collisionResult.getContacts(contacts);
    if (contacts.size() != 1u) {
      GTEST_FAIL() << "\nFundamental error! "
                   << "Reported wrong number of collisions; expected 1 got "
                   << contacts.size();
    }

    // For this test to pass, at least one quantity must lie outside the
    // threshold. Each following tests test for *non-compliant* values.
    bool passed = false;

    const Contact<S>& contact = contacts[0];
    passed |= std::abs(contact.penetration_depth - depth) > test_tolerance;

    const Vector3<S> normal{0, 0, 1};
    passed |= !CompareMatrices(contact.normal, normal, test_tolerance);

    // The contact position should be 1/2 the penetration depth below the
    // sphere's extent along the +z axis. This maximal extent is (0, 0, 0).
    const Vector3<S> contact_point{0, 0, -depth / 2};
    passed |= !CompareMatrices(contact.pos, contact_point, test_tolerance);

    if (!passed) {
      GTEST_FAIL() << "All quantities passed within the given threshold: "
                   << test_tolerance
                   << "\n\tNormal - expected: " << normal.transpose()
                   << ", tested: " << contact.normal.transpose()
                   << "\n\tPos - expected: " << contact_point.transpose()
                   << ", tested: " << contact.pos.transpose()
                   << "\n\tDepth - expected: " << depth
                   << ", tested: " << contact.penetration_depth;
    }
  }

  // This asserts that the collision values all lie within the given threshold.
  void AssertCollisionPrecision(CollisionObject<S>* tangent_object, S depth,
                                S test_tolerance) {
    CollisionResult<S> collisionResult;
    PerformCollision(tangent_object, depth, &collisionResult);

    ASSERT_TRUE(collisionResult.isCollision());
    std::vector<Contact<S>> contacts;
    collisionResult.getContacts(contacts);
    GTEST_ASSERT_EQ(contacts.size(), 1u);

    const Contact<S>& contact = contacts[0];
    EXPECT_NEAR(contact.penetration_depth, depth, test_tolerance);
    const Vector3<S> normal{0, 0, 1};
    EXPECT_TRUE(CompareMatrices(contact.normal, normal, test_tolerance))
        << "\texpected normal: " << normal.transpose()
        << "\n\ttested normal: " << contact.normal.transpose()
        << "\n\ttolerance: " << test_tolerance;
    // The contact position should be 1/2 the penetration depth below the
    // sphere's extent along the +z axis. This maximal extent is (0, 0, 0).
    const Vector3<S> contact_point{0, 0, -depth / 2};
    EXPECT_TRUE(CompareMatrices(contact.pos, contact_point, test_tolerance))
        << "\texpected point: " << contact_point.transpose()
        << "\n\ttested point: " << contact.pos.transpose()
        << "\n\ttolerance: " << test_tolerance;
  }

 private:
  // Initialize a unit sphere tangent to the contact plane.
  void InitializeTangentSphere() {
    const S radius = 1;
    Transform3<S> sphere_pose{Translation3<S>{0, 0, -radius}};
    CollisionGeometryPtr_t<S> sphere_geometry(new Sphere<S>(radius));
    tangent_sphere_.reset(new CollisionObject<S>(sphere_geometry, sphere_pose));
  }

  // Initializes a plane coincident with the abstract tangent plane.
  void InitializeTangentPlane() {
    CollisionGeometryPtr_t<S> plane_geometry(
        new Halfspace<S>(Vector3<S>::UnitZ(), 0));
    tangent_plane_.reset(
        new CollisionObject<S>(plane_geometry, Transform3<S>::Identity()));
  }

  // Initialize a large, axis-aligned box whose +z face lies on the z = 0 plane.
  void InitializeTangentBox() {
    const S size{30};
    Transform3<S> box_pose = Transform3<S>::Identity();
    box_pose.translation() << 0, 0, -size / 2;
    CollisionGeometryPtr_t<S> box_geometry(new Box<S>(size, size, size));
    tangent_box_.reset(new CollisionObject<S>(box_geometry, box_pose));
  }

  // Does the work of configuring the box to penetrate the plane with the given
  // depth and colliding the tangent object against the box.
  void PerformCollision(CollisionObject<S>* tangent_object, S depth,
                        CollisionResult<S>* results) {
    // Offset the box by the depth to create penetration.
    Vector3<S> offset{0, 0, -depth};
    box_->setTranslation(box_offset_ + offset);

    // Compute collision - single contact and enable contact.
    CollisionRequest<S> collisionRequest(1 /* num contacts */,
                                         true /* contacts enabled */);
    collisionRequest.gjk_solver_type = solver_type_;
    collide(tangent_object, box_.get(), collisionRequest, *results);
  }

  // The offset that puts the target vertex at (0, 0, 0);
  GJKSolverType solver_type_;
  std::unique_ptr<CollisionObject<S>> box_{nullptr};
  // The tangent geometries.
  std::unique_ptr<CollisionObject<S>> tangent_sphere_{nullptr};
  std::unique_ptr<CollisionObject<S>> tangent_plane_{nullptr};
  std::unique_ptr<CollisionObject<S>> tangent_box_{nullptr};
  // The offset to the rotated box so that it's colliding corner is sitting
  // on the origin. Moving below the z=0 plane will lead to penetrations.
  Vector3<S> box_offset_;
  // Default tolerance - used to determine that the penetrating corner is
  // *truly* at the origin.
  static const S kDefaultTolerance;
};

// Specify per-scalar default test tolerances.
template <typename S>
const S BoxPrecisionTest<S>::kDefaultTolerance = S(1e-6);

template <>
const double BoxPrecisionTest<double>::kDefaultTolerance = 1e-14;

template <>
const float BoxPrecisionTest<float>::kDefaultTolerance = 1e-7;

// This tests collision under shallow collision. All tangent geometries should
// report the same collision.
template <typename S>
void test_box_shallow_collision(GJKSolverType solver_type) {
  const std::vector<S> depths{1e-10, 1e-8, 1e-6};
  // Perform the tests.
  BoxPrecisionTest<S> test(solver_type);
  for (const auto& depth : depths) {
    test.AssertCollision(test.tangent_plane(), depth, true, "plane");
    test.AssertCollision(test.tangent_box(), depth, true, "box");
    test.AssertCollision(test.tangent_sphere(), depth, true, "sphere");
  }
}

// This test illustrates the disparity in precision in the contact solutions
// across solutions. For a given depth, it brackets the precision of the
// reported contact data.
template <typename S>
void test_box_precision_threshold(GJKSolverType solver_type) {
  BoxPrecisionTest<S> test(solver_type);
  const double depth = 1e-6;

  // Pairs of tests per geometry
  //  1. The largest precision threshold for which the reported answer does
  //     *not* match the expected answer.
  //  2. The smallest precision threshold for which the reported answer does
  //     match.

  // The sphere exhibits *very* low precision. NOTE: this level of precision
  // bounds the *worst* of the two solvers and doesn't suggest that the solvers
  // are necessarily equivalent.
  // NOTE: Both solvers will fail for a precision of 1e-4 or smaller.
  test.AssertInsufficientPrecision(test.tangent_sphere(), depth, 1e-4);
  test.AssertCollisionPrecision(test.tangent_sphere(), depth, 3e-2);

  // The box and plane have full precision in this scenario.
  test.AssertInsufficientPrecision(test.tangent_plane(), depth, 1e-16);
  test.AssertCollisionPrecision(test.tangent_plane(), depth, 1e-15);
  test.AssertInsufficientPrecision(test.tangent_box(), depth, 1e-16);
  test.AssertCollisionPrecision(test.tangent_box(), depth, 1e-15);
}

GTEST_TEST(FCL_BOX_COLLISION_PRECISION, shallow_collision_ccd) {
  test_box_shallow_collision<double>(GJKSolverType::GST_LIBCCD);
}

GTEST_TEST(FCL_BOX_COLLISION_PRECISION, bracketed_precision_ccd) {
  test_box_precision_threshold<double>(GJKSolverType::GST_LIBCCD);
}

GTEST_TEST(FCL_BOX_COLLISION_PRECISION, shallow_collision_indep) {
  test_box_shallow_collision<double>(GJKSolverType::GST_INDEP);
}

GTEST_TEST(FCL_BOX_COLLISION_PRECISION, bracketed_precision_indep) {
  test_box_precision_threshold<double>(GJKSolverType::GST_INDEP);
}

//==============================================================================
int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
