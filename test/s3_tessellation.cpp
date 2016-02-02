/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include "bbTrans/s3_tessellation.h"
using namespace bb;

int main(int argc, char**argv) {
  std::vector<Tetrahedron4D> tetrahedra = TessellateS3();
  std::cout << tetrahedra.size() << std::endl;
}
