
/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "bbTrans/node_S3.h"
#include "bbTrans/numeric_helpers.h"
#include "bbTrans/vmf.h"
#include "bbTrans/vmf_mm.h"
#include "bbTrans/bound.h"

namespace bb {

class UpperBoundIndepS3 : public Bound<NodeS3> {
 public:
  UpperBoundIndepS3(const vMFMM<3>& vmf_mm_A, const vMFMM<3>& vmf_mm_B);
  virtual double Evaluate(const NodeS3& node);
  virtual double EvaluateAndSet(NodeS3& node);
  virtual double EvaluateRotationSet(const
      std::vector<Eigen::Quaterniond>& qs) const;
 protected:
  const vMFMM<3>& vmf_mm_A_;
  const vMFMM<3>& vmf_mm_B_;

};

template<uint32_t k> 
void ComputeMaxJofSubset(const Eigen::Matrix<double,3,Eigen::Dynamic>& M, 
    const Eigen::Vector3d& mu, bool verbose,
    double& Jmax,
    Eigen::Vector3d& pMax ) {
  Combinations combNKs(M.cols(),k);
  for (auto comb : combNKs.Get()) {
    Eigen::Matrix<double,3,k> M_; 
    for (uint32_t i=0; i<k; ++i) {
      M_.col(k) = M.col(comb[k]);
    }
    Eigen::Matrix<double,3,3> MTM = M.transpose()*M;
    Eigen::Matrix<double,3,1> MTmu = M.transpose()*mu;
    Eigen::Vector3d alpha = MTM.ldlt().solve(MTmu);
    double sigma = std::numeric_limits<double>::lowest();
    if ((MTmu.array() >= 0.).all())
      sigma = 1.;
    else if ((MTmu.array() < 0.).all())
      sigma = -1.;
    double J =  sigma*sqrt(MTmu.dot(alpha)) ;
    if (J > Jmax) {
      Jmax = J;
      pMax = M * alpha;
    }
  }
}

Eigen::Vector3d ComputeExtremumOnGeodesic(const Eigen::Vector3d& q1,
    const Eigen::Vector3d& q2, const Eigen::Vector3d& p, bool verbose);

Eigen::Vector3d ClosestPointInRotationSetOld(const vMF<3>& vmf_A, const
    vMF<3>& vmf_B, const std::vector<Eigen::Quaterniond>& qs, bool
    furthest=false, bool verbose=false);

Eigen::Vector3d ClosestPointInRotationSet(const vMF<3>& vmf_A, const
    vMF<3>& vmf_B, const std::vector<Eigen::Quaterniond>& qs, bool
    furthest=false, bool verbose=false);

Eigen::Vector3d FurthestPointInRotationSet(const vMF<3>& vmf_A, const
    vMF<3>& vmf_B, const std::vector<Eigen::Quaterniond>& qs, 
    bool verbose);

/// This function just calls ClosestPointInRotationSet with the
/// rotations implied by Tetrahedron.
Eigen::Vector3d ClosestPointInTetrahedron(const vMF<3>& vmf_A, const
    vMF<3>& vmf_B, const Tetrahedron4D& tetrahedron, bool
    furthest=false, bool verbose = false);

Eigen::Vector3d FurthestPointInTetrahedron(const vMF<3>& vmf_A, const
    vMF<3>& vmf_B, const Tetrahedron4D& tetrahedron, bool verbose = false);

}
