#include <gtest/gtest.h>

#include "chan3d.hpp"

// Just in case
std::istream &operator>>(std::istream &in, std::vector<Point> &points) {
  for (std::size_t i = 0; i < points.size(); ++i) {
    in >> points[i].x >> points[i].y >> points[i].z;
  }
  return in;
}

// Just in case
std::ostream &operator<<(std::ostream &out, DistancesHull3D &dist) {
  std::vector<double> distances = dist.FindDistancesToHull();

  out << std::fixed << std::setprecision(4);

  for (auto distance : distances) {
    out << distance << '\n';
  }

  return out;
}

TEST(Simple, Tetrahedron) {
  std::vector<Point> tetrahedron = {{0, 0, 0}
                                  , {100, 0, 0}
                                  , {0, 100, 0}
                                  , {0, 0, 100}
                                  , {20, 20, 20}
                                  , {30, 20, 10}};
  std::vector<Point> requests = {{1, 1, 1}
                               , {30, 30, 35}
                               , {7, 8, 9}
                               , {90, 2, 2}};
  DistancesHull3D distances(tetrahedron, requests);
  
  std::vector<double> result = distances.FindDistancesToHull();
  std::vector<double> correct = {1.0000000000
                               , 2.8867513459
                               , 7.0000000000
                               , 2.0000000000};
  
  for (std::size_t i = 0; i < correct.size(); ++i) {
    ASSERT_NEAR(result[i], correct[i], 1e-4);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  
  return RUN_ALL_TESTS();
}
