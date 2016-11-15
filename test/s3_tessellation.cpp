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

  double V = 0.;
  for (const auto& tetra : tetrahedra) {
    V += tetra.GetVolume();
//    std::cout << tetra.GetCenter().transpose() << std::endl;
  }
  std::cout << "Tessellation l0 volume: " << V << std::endl;

//  V = 0.;
//  for (const auto& tetra : tetrahedra) {
//    auto l1 = tetra.Subdivide();
//    for (const auto& t : l1) {
//      V += t.GetVolume();
//    }
//  }
//  std::cout << "Tessellation l1 volume: " << V << std::endl;
//
//  V = 0.;
//  for (const auto& tetra : tetrahedra) {
//    auto l1 = tetra.Subdivide();
//    for (const auto& t1 : l1) {
//      auto l2 = t1.Subdivide();
//      for (const auto& t2 : l2) {
//        V += t2.GetVolume();
//      }
//    }
//  }
//  std::cout << "Tessellation l2 volume: " << V << std::endl;
//
//  V = 0.;
//  for (const auto& tetra : tetrahedra) {
//    auto l1 = tetra.Subdivide();
//    for (const auto& t1 : l1) {
//      auto l2 = t1.Subdivide();
//      for (const auto& t2 : l2) {
//        auto l3 = t2.Subdivide();
//        for (const auto& t3 : l3) {
//          V += t3.GetVolume();
//        }
//      }
//    }
//  }
//  std::cout << "Tessellation l3 volume: " << V << std::endl;
  std::cout << "true volume of half sphere: " << M_PI*M_PI << std::endl;

  uint32_t N = 0;
  uint32_t Nall = 100000;
  for (uint32_t i=0; i<Nall; ++i) {
    S4d q = S4d::Random();
//    if (i%100 == 0) 
//      std::cout << q.vector().transpose() << std::endl;
    for (const auto& tetra : tetrahedra) 
      if (tetra.Intersects(q.vector())) {
        ++N;
//        break;
      }
  }
  std::cout << "fraction all over sphere intersected with tessellation: "
    << static_cast<double>(N)/static_cast<double>(Nall)
    << std::endl;

  N = 0.;
  for (uint32_t i=0; i<Nall; ++i) {
    S4d q = S4d::Random();
    q.vector()(3) = q.vector()(3) < 0. ? -q.vector()(3) : q.vector()(3);
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
    q.vector()(3) = q.vector()(3) < 0. ? q.vector()(3) : -q.vector()(3);
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
