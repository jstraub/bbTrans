/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include "optRot/tetrahedron.h"

Eigen::Vector4d normed(const Eigen::Vector4d& x) {
  return x / x.norm();
}

Tetrahedron4D::Tetrahedron4D(const Eigen::Matrix<double, 4, 4>& vertices) :
  vertices_(vertices) {}

Tetrahedron4D::Tetrahedron4D(const Eigen::Vector4d& a, const
    Eigen::Vector4d& b, const Eigen::Vector4d& c, const
    Eigen::Vector4d& d) {
  vertices_ << a, b, c, d;
}

Eigen::Vector4d Tetrahedron4D::GetCenter() {
  return normed(vertices_.rowwise().sum());
}

std::vector<Tetrahedron4D> Tetrahedron4D::Subdivide() {
  std::vector<Tetrahedron4D> tetrahedra;  
  tetrahedra.reserve(8);
  // Compute new vertices and "pop" them out to the sphere.
  Eigen::Matrix<double, 4, 6> vertices;
  vertices << normed(vertices_.col(0) + vertices_.col(1)),
    normed(vertices_.col(1) + vertices_.col(2)),
    normed(vertices_.col(2) + vertices_.col(0)),
    normed(vertices_.col(0) + vertices_.col(3)),
    normed(vertices_.col(1) + vertices_.col(3)),
    normed(vertices_.col(2) + vertices_.col(3));
  // Corner tetrahedron at 0th corner of parent.
  tetrahedra.push_back(Tetrahedron4D(vertices_.col(0), vertices.col(0),
        vertices.col(2), vertices.col(3)));
  // Corner tetrahedron at 1th corner of parent.
  tetrahedra.push_back(Tetrahedron4D(vertices_.col(1), vertices.col(0),
        vertices.col(1), vertices.col(4)));
  // Corner tetrahedron at 2th corner of parent.
  tetrahedra.push_back(Tetrahedron4D(vertices_.col(2), vertices.col(1),
        vertices.col(2), vertices.col(5)));
  // Corner tetrahedron at 3th corner of parent.
  tetrahedra.push_back(Tetrahedron4D(vertices_.col(3), vertices.col(3),
        vertices.col(4), vertices.col(5)));
  Eigen::Vector3d dots;
  dots[0] = vertices.col(0).transpose() * vertices.col(5);
  dots[1] = vertices.col(2).transpose() * vertices.col(4);
  dots[2] = vertices.col(3).transpose() * vertices.col(1);
  uint32_t skewEdgeId = 0;
  dots.maxCoeff(&skewEdgeId);
  if (skewEdgeId == 0) {
    tetrahedra.push_back(Tetrahedron4D(vertices.col(0), vertices.col(5),
        vertices.col(3), vertices.col(2)));
    tetrahedra.push_back(Tetrahedron4D(vertices.col(0), vertices.col(5),
        vertices.col(3), vertices.col(4)));
    tetrahedra.push_back(Tetrahedron4D(vertices.col(0), vertices.col(5),
        vertices.col(1), vertices.col(4)));
    tetrahedra.push_back(Tetrahedron4D(vertices.col(0), vertices.col(5),
        vertices.col(1), vertices.col(2)));
  } else if (skewEdgeId == 1) {
    tetrahedra.push_back(Tetrahedron4D(vertices.col(2), vertices.col(4),
        vertices.col(3), vertices.col(0)));
    tetrahedra.push_back(Tetrahedron4D(vertices.col(2), vertices.col(4),
        vertices.col(0), vertices.col(1)));
    tetrahedra.push_back(Tetrahedron4D(vertices.col(2), vertices.col(4),
        vertices.col(1), vertices.col(5)));
    tetrahedra.push_back(Tetrahedron4D(vertices.col(2), vertices.col(4),
        vertices.col(3), vertices.col(5)));
  } else if (skewEdgeId == 2) {
    tetrahedra.push_back(Tetrahedron4D(vertices.col(3), vertices.col(1),
        vertices.col(0), vertices.col(2)));
    tetrahedra.push_back(Tetrahedron4D(vertices.col(3), vertices.col(1),
        vertices.col(0), vertices.col(4)));
    tetrahedra.push_back(Tetrahedron4D(vertices.col(3), vertices.col(1),
        vertices.col(5), vertices.col(4)));
    tetrahedra.push_back(Tetrahedron4D(vertices.col(3), vertices.col(1),
        vertices.col(5), vertices.col(2)));
  }
}
