
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
  Eigen::Matrix<double,3,k> M_ = Eigen::Matrix<double,3,k>::Zero(); 
  for (const auto& comb : combNKs.Get()) {
    for (uint32_t i=0; i<comb.size(); ++i) {
      M_.col(i) = M.col(comb[i]);
    }

    Eigen::Matrix<double,k,k> MTM = M_.transpose()*M_;
    Eigen::Matrix<double,k,1> MTmu = M_.transpose()*mu;
    // slow and stable
//    Eigen::Matrix<double,k,1> x = MTM.fullPivHouseholderQr().solve(MTmu);
//    Eigen::Matrix<double,k,1> x = MTM.colPivHouseholderQr().solve(MTmu);
//    Eigen::Matrix<double,k,1> x = MTM.householderQr().solve(MTmu);
    // fastest not as stable
    Eigen::Matrix<double,k,1> x = MTM.ldlt().solve(MTmu);
    if (verbose)
      std::cout << M_ << std::endl
        << MTM << std::endl
        << MTmu.transpose() << std::endl;
    double sigma = std::numeric_limits<double>::lowest();
    if ((x.array() >= 0.).all())
      sigma = 1.;
    else if ((x.array() < 0.).all())
      sigma = -1.;
    double J = sigma*sqrt(MTmu.dot(x)) ;
    if (verbose)
      std::cout << "J=" << J << " x: " << x.transpose() << std::endl;
    if (J > Jmax) {
      Jmax = J;
//      pMax = ((M_ * x) / J).normalized();
      pMax = (sigma*M_*x).normalized();
    }
  }
}

Eigen::Vector3d ComputeExtremumOnGeodesic(const Eigen::Vector3d& q1,
    const Eigen::Vector3d& q2, const Eigen::Vector3d& p, bool verbose);


Eigen::Vector3d ClosestPointInRotationSetNew(const vMF<3>& vmf_A, const
    vMF<3>& vmf_B, const std::vector<Eigen::Quaterniond>& qs, bool
    furthest=false, bool verbose=false);

Eigen::Vector3d ClosestPointInRotationSetOld(const vMF<3>& vmf_A, const
    vMF<3>& vmf_B, const std::vector<Eigen::Quaterniond>& qs, bool
    furthest=false, bool verbose=false);

inline Eigen::Vector3d ClosestPointInRotationSet(const vMF<3>& vmf_A, const
    vMF<3>& vmf_B, const std::vector<Eigen::Quaterniond>& qs, bool
    furthest=false, bool verbose=false) {

//  verbose = true;
  Eigen::Vector3d pNew;
//  Eigen::Vector3d pOld;
//  if (verbose) 
//    std::cout << " ---- old method " << std::endl;
//  pOld = ClosestPointInRotationSetOld(vmf_A, vmf_B, qs,
//      furthest, verbose);
//  if (verbose) 
//    std::cout << " ---- new method " << std::endl;
  pNew = ClosestPointInRotationSetNew(vmf_A, vmf_B, qs,
      furthest, verbose);
//  if (verbose && pOld.dot(pNew) < 0.99)  {
//    std::cout << "old: " << pOld.transpose() 
//      << " new " << pNew.transpose() 
//      << " dot: " << pOld.dot(pNew)
//      << (furthest ? " furthest" : " closest")
//      << std::endl;
//  }
  return pNew;
//  return pOld;
}

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
