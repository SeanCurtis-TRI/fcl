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

#include <gtest/gtest.h>

#include "fcl/math/constants.h"
#include "fcl/narrowphase/collision.h"
#include "fcl/narrowphase/collision_object.h"
#include "fcl/narrowphase/distance.h"

/** @author Sean Curtis */

// This is a scenario that arose in exercising FCL in another program.
// https://github.com/robotlocomotion/drake/pull/10578
// In short: with v2.1 of libccd, the collision between a touching box and
// cylinder reports a functionally zero distance, a contact point on the
// surface, but a *zero* normal. This is captured in issue
// https://github.com/flexible-collision-library/fcl/issues/375
GTEST_TEST(FCL_BOX_CYLINDER, IncompleteCollisionResult) {
  using CollisionGeometryPtr_t = std::shared_ptr<fcl::CollisionGeometryd>;
  CollisionGeometryPtr_t box_geom{
      new fcl::Boxd(fcl::Vector3d{0.49, 0.63, 0.015})};
  CollisionGeometryPtr_t cylinder_geom{new fcl::Cylinderd{0.08, 0.002}};

  fcl::Transform3d X_WB, X_WC;
  // clang-format off
  X_WB.matrix() << 0, -1, 0, -0.25,
                   1,  0, 0,     0,
                   0,  0, 1,     0,
                   0,  0, 0,     1;
  X_WC.matrix() << 1, 0, 0,      0,
                   0, 1, 0,      0,
                   0, 0, 1, 0.0085,
                   0, 0, 0,      1;
  // clang-format on
  fcl::CollisionObjectd box(box_geom, X_WB);
  fcl::CollisionObjectd cylinder(cylinder_geom, X_WC);

  fcl::CollisionRequestd collisionRequest(1, true);
  collisionRequest.gjk_solver_type = fcl::GJKSolverType::GST_LIBCCD;
  fcl::CollisionResultd collisionResult;

  fcl::collide(&box, &cylinder, collisionRequest, collisionResult);

  const double kEps = fcl::constants<double>::eps();
  ASSERT_TRUE(collisionResult.isCollision());
  std::vector<fcl::Contactd> contacts;
  collisionResult.getContacts(contacts);
  GTEST_ASSERT_EQ(contacts.size(), 1u);
  EXPECT_NEAR(contacts[0].penetration_depth, 0, kEps);
  EXPECT_NEAR(contacts[0].pos(2), 0.0075, kEps);
  // This is the error! This should *not* be a zero normal.
  EXPECT_NEAR(contacts[0].normal.norm(), 0, kEps);
}

//==============================================================================
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
