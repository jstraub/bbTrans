
#include <Eigen/Dense>

#include "tetrahedron.h"

std::vector<Tetrahedron4D> TesselateS3() {
  std::vector<Tetrahedron4D> tetrahedra;
  tetrahedra.reserve(600);
  
  Eigen::Matrix<double, 4, 120> vertices;
  uint32_t i = 0;
  for (double a=0; a<2; ++a) 
    for (double b=0; b<2; ++b) 
      for (double c=0; c<2; ++c) 
        for (double d=0; d<2; ++d) {
          vertices(0,i) = a*2. - 1.;
          vertices(1,i) = b*2. - 1.;
          vertices(2,i) = c*2. - 1.;
          vertices(3,i++) = d*2. - 1.;
        }
  for (uint32_t j=0; j<4; ++j) {
    vertices(j,i++) = 2.;
    vertices(j,i++) = -2.;
  }
  // Golden Ratio 
  double phi = 0.5 * (1. + sqrt(5.));
  // All even permutations
  // http://mathworld.wolfram.com/EvenPermutation.html
  Eigen::Matrix<uint32_t, 4, 12> evenPerms;
  evenPerms << 
    0,1,2,3,
    0,2,3,1,
    0,3,1,2,
    1,0,3,2,
    1,2,0,3,
    1,3,2,0,
    2,0,1,3,
    2,1,3,0,
    2,3,0,1,
    3,0,2,1,
    3,1,0,2,
    3,2,1,0;
  for (uint32_t j=0; j<12; ++j) 
    for (double a=0; a<2; ++a) 
      for (double b=0; b<2; ++b) 
        for (double c=0; c<2; ++c) {
          vertices(evenPerms(0,j),i) = (a*2.-1.)*phi;
          vertices(evenPerms(1,j),i) = (b*2.-1.);
          vertices(evenPerms(2,j),i) = (c*2.-1.)/phi;
          vertices(evenPerms(3,j),i++) = 0.;
        }
  vertices *= 0.5;
  assert(i == 120);
  // Filter out half of the sphere.
  Eigen::Vector4d north;
  north << 1., 0., 0., 0.;
  uint32_t j = 0;
  for (i = 0; i < 120; ++i) 
    if (acos(north.transpose() * vertices.col(i)) <= 120.*M_PI/180.){
      vertices.col(j++) = vertices.col(i); 
    }

  return tetrahedra;
}
