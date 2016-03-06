/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include "manifold/S.h"
#include "bbTrans/s3_tessellation.h"
using namespace bb;

int main(int argc, char**argv) {
  std::vector<Tetrahedron4D> tetrahedra = TessellateS3();
  std::cout << tetrahedra.size() << std::endl;

  uint32_t N = 0;
  uint32_t Nall = 100000;
  for (uint32_t i=0; i<Nall; ++i) {
    S4d q = S4d::Random();
    for (const auto& tetra : tetrahedra) 
      if (tetra.Intersects(q.vector())) {
        ++N;
        break;
      }
  }
  std::cout << "fraction all over sphere intersected with tessellation: "
    << static_cast<double>(N)/static_cast<double>(Nall)
    << std::endl;
  N = 0.;
  for (uint32_t i=0; i<Nall; ++i) {
    S4d q = S4d::Random();
    q.vector()(0) = q.vector()(0) < 0. ? -q.vector()(0) : q.vector()(0);
    for (const auto& tetra : tetrahedra) 
      if (tetra.Intersects(q.vector())) {
        ++N;
        break;
      }
  }
  std::cout << "fraction on top half-sphere intersected with tessellation: "
    << static_cast<double>(N)/static_cast<double>(Nall)
    << std::endl;
  N = 0.;
  for (uint32_t i=0; i<Nall; ++i) {
    S4d q = S4d::Random();
    q.vector()(0) = q.vector()(0) < 0. ? q.vector()(0) : -q.vector()(0);
    for (const auto& tetra : tetrahedra) 
      if (tetra.Intersects(q.vector())) {
        ++N;
        break;
      }
  }
  std::cout << "fraction on bottom half-sphere intersected with tessellation: "
    << static_cast<double>(N)/static_cast<double>(Nall)
    << std::endl;
}
