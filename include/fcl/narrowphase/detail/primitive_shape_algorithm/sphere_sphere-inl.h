/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011-2014, Willow Garage, Inc.
 *  Copyright (c) 2014-2016, Open Source Robotics Foundation
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
 *   * Neither the name of Open Source Robotics Foundation nor the names of its
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

/** @author Jia Pan */

#ifndef FCL_NARROWPHASE_DETAIL_SPHERESPHERE_INL_H
#define FCL_NARROWPHASE_DETAIL_SPHERESPHERE_INL_H

#include "fcl/narrowphase/detail/primitive_shape_algorithm/sphere_sphere.h"

namespace fcl
{

namespace detail
{

//==============================================================================
extern template
FCL_EXPORT
bool sphereSphereIntersect(const Sphere<double>& s1, const Transform3<double>& tf1,
                           const Sphere<double>& s2, const Transform3<double>& tf2,
                           std::vector<ContactPoint<double>>* contacts);

//==============================================================================
extern template
FCL_EXPORT
bool sphereSphereDistance(const Sphere<double>& s1, const Transform3<double>& tf1,
                          const Sphere<double>& s2, const Transform3<double>& tf2,
                          double* dist, Vector3<double>* p1, Vector3<double>* p2);

//==============================================================================
extern template
FCL_EXPORT
double sphereSphereSignedDistance(
    const Sphere<double>& s1, const Transform3<double>& tf1,
    const Sphere<double>& s2, const Transform3<double>& tf2,
    Vector3<double>* p1, Vector3<double>* p2);

//==============================================================================
template <typename S>
FCL_EXPORT
bool sphereSphereIntersect(const Sphere<S>& s1, const Transform3<S>& tf1,
                           const Sphere<S>& s2, const Transform3<S>& tf2,
                           std::vector<ContactPoint<S>>* contacts)
{
  Vector3<S> diff = tf2.translation() - tf1.translation();
  S len = diff.norm();
  if(len > s1.radius + s2.radius)
    return false;

  if(contacts)
  {
    // If the centers of two sphere are at the same position, the normal is (0, 0, 0).
    // Otherwise, normal is pointing from center of object 1 to center of object 2
    const Vector3<S> normal = len > 0 ? (diff / len).eval() : diff;
    const Vector3<S> point = tf1.translation() + diff * s1.radius / (s1.radius + s2.radius);
    const S penetration_depth = s1.radius + s2.radius - len;
    contacts->emplace_back(normal, point, penetration_depth);
  }

  return true;
}

//==============================================================================
template <typename S>
FCL_EXPORT
bool sphereSphereDistance(const Sphere<S>& s1, const Transform3<S>& tf1,
                          const Sphere<S>& s2, const Transform3<S>& tf2,
                          S* dist, Vector3<S>* p1, Vector3<S>* p2)
{
  Vector3<S> o1 = tf1.translation();
  Vector3<S> o2 = tf2.translation();
  Vector3<S> diff = o1 - o2;
  S len = diff.norm();
  if(len > s1.radius + s2.radius)
  {
    if(dist) *dist = len - (s1.radius + s2.radius);
    if(p1) *p1 = (o1 - diff * (s1.radius / len));
    if(p2) *p2 = (o2 + diff * (s2.radius / len));
    return true;
  }

  if(dist) *dist = -1;
  return false;
}

template <typename S>
FCL_EXPORT
S sphereSphereSignedDistance(const Sphere<S>& s1, const Transform3<S>& X_WS1,
                             const Sphere<S>& s2, const Transform3<S>& X_WS2,
                             Vector3<S>* p_WN1, Vector3<S>* p_WN2)
{
  Vector3<S> p_WS1 = X_WS1.translation();
  Vector3<S> p_WS2 = X_WS2.translation();
  Vector3<S> p_S2S1_W = p_WS1 - p_WS2;
  S center_distance = p_S2S1_W.norm();

  // Negative for penetration, positive for separated, zero for touching.
  S distance = center_distance - (s1.radius + s2.radius);
  if (p_WN1 || p_WN2) {
    const S eps = constants<S>::eps();
    // Unit normal pointing in the direction from S2's center to S1's center.
    const Vector3<S> n_S2S1 = center_distance < eps
                                  ? Vector3<S>{1, 0, 0}
                                  : (p_S2S1_W / center_distance).eval();
    if (p_WN1) *p_WN1 = p_WS1 - n_S2S1 * s1.radius;
    if (p_WN2) *p_WN2 = p_WS2 + n_S2S1 * s2.radius;
  }
  return distance;
}

} // namespace detail
} // namespace fcl

#endif
