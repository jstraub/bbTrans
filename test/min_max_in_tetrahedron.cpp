/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "optRot/node.h"
#include "optRot/lower_bound_S3.h"
#include "optRot/upper_bound_indep_S3.h"
#include "optRot/upper_bound_convex_S3.h"

using namespace OptRot;

bool VectorsClose(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
  return ((a-b).array().abs() < 1.e-6).all();
}

int main(int argc, char ** argv) {
  std::cout << " Testing ComputeExtremumOnGeodesic() " << std::endl;

  Eigen::Vector3d q1, q2, p;
  q1 << 1., 0., 0.;
  q2 << 0., 1., 0.;
  p = q1;
  Eigen::Vector3d a = ComputeExtremumOnGeodesic(q1, q2, p);
  if (!VectorsClose(a,q1)) 
    std::cout << " ERROR: " << p.transpose();
  p = q2;
  a = ComputeExtremumOnGeodesic(q1, q2, p);
  if (!VectorsClose(a,q2)) 
    std::cout << " ERROR: " << p.transpose();
  p << -1.,0.,0.;
  a = ComputeExtremumOnGeodesic(q1, q2, p);
  if (!VectorsClose(a,q2)) 
    std::cout << " ERROR: " << p.transpose();
  p << 0.,-1.,0.;
  a = ComputeExtremumOnGeodesic(q1, q2, p);
  if (!VectorsClose(a,q1)) 
    std::cout << " ERROR: " << p.transpose();
}
