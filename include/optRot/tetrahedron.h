
#include <vector>
#include <Eigen/Dense>

Vector4d normed(const Vector4d& x);

class Tetrahedron4D {
 public:
  Tetrahedron4D(const Eigen::Matrix<double, 4, 4>& vertices);
  Tetrahedron4D(const Eigen::Vector4d& a, const Eigen::Vector4d& b,
      const Eigen::Vector4d& c, const Eigen::Vector4d& d);
  ~Tetrahedron4D() = default;
  Vector4d GetCenter();
  std::vector<Tetrahedron4D> Subdivide();
 private:
  /// One 4D vertex per column. 4 vertices in total to describe the 4D
  /// Tetrahedron.
  Eigen::Matrix<double, 4, 4> vertices_;
};
