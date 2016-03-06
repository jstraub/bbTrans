/* Copyright (c) 2016, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "manifold/S.h"
#include "bbTrans/combinations.h"
#include "bbTrans/tetrahedron.h"

namespace bb {

std::vector<Tetrahedron4D> TessellateS3();
void TessellationTest(std::vector<Tetrahedron4D>& tetrahedra, uint32_t Nsamples);

}
