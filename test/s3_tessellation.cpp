/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include "optRot/s3_tessellation.h"
using namespace OptRot;

int main(int argc, char**argv) {
  std::vector<Tetrahedron4D> tetrahedra = TessellateS3();
  std::cout << tetrahedra.size() << std::endl;
}
