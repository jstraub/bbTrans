/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <vector>
#include <Eigen/Dense>

namespace OptRot {
Eigen::Vector4d normed(const Eigen::Vector4d& x);

class Tetrahedron4D {
 public:
  Tetrahedron4D(const Eigen::Matrix<double, 4, 4>& vertices);
  Tetrahedron4D(const Eigen::Vector4d& a, const Eigen::Vector4d& b,
      const Eigen::Vector4d& c, const Eigen::Vector4d& d);
  ~Tetrahedron4D() = default;

  Eigen::Vector4d GetCenter() const;
  Eigen::Quaterniond GetCenterQuaternion() const;
  Eigen::Vector4d GetVertex(uint32_t i) const;
  Eigen::Quaterniond GetVertexQuaternion(uint32_t i) const;
  std::vector<Tetrahedron4D> Subdivide() const;
 private:
  /// One 4D vertex per column. 4 vertices in total to describe the 4D
  /// Tetrahedron.
  Eigen::Matrix<double, 4, 4> vertices_;
};
}
