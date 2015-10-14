
#include <iostream>
#include "optRot/s3_tessellation.h"

int main(int argc, char**argv) {
  std::vector<Tetrahedron4D> tetrahedra = TessellateS3();
  std::cout << tetrahedra.size() << std::endl;
}
