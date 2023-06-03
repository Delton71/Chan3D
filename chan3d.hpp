#ifndef CHAN3D_HPP
#define CHAN3D_HPP

// A Minimalist’s Implementation of the 3-d Divide-and-Conquer Convex Hull
// Algorithm Source: http://tmc.web.engr.illinois.edu/ch3d/ch3d.pdf

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

struct Point {
  double x{0.0};
  double y{0.0};
  double z{0.0};

  // prev and next points in tuple {prev, this, next}
  Point *prev{nullptr};
  Point *next{nullptr};

  void ModifyRelationship() {
    if (prev->next != this) {
      // INSERTION
      prev->next = this;
      next->prev = this;
    } else {
      // DELETION
      prev->next = next;
      next->prev = prev;
    }
  }
};

// to find the convex hull
class Hull3D {
public:
  explicit Hull3D(const std::vector<Point> &points)
      : points_(points), inv_points_(0) {}
  Hull3D(const Hull3D &other) : points_(other.points_), inv_points_(0) {}

  constexpr std::size_t Size() const { return points_.size(); }
  std::vector<Point> FindHull() {
    hull_.clear();
    CreateHull();
    return hull_;
  }

  // To rotate all points around the x, y, z axes by degrees
  constexpr static const double kRotateDegree = 1.5;
  static void RotatePoints(std::vector<Point> &points, double x_deg,
                           double y_deg, double z_deg) {
    double cos_x = std::cos(x_deg);
    double sin_x = std::sin(x_deg);
    double cos_y = std::cos(y_deg);
    double sin_y = std::sin(y_deg);
    double cos_z = std::cos(z_deg);
    double sin_z = std::sin(z_deg);

    double rot_m[3][3] = {{cos_x * cos_y, cos_x * sin_y * sin_z,
                           cos_x * sin_y * cos_z + sin_x * sin_z},
                          {sin_x * cos_y, sin_x * sin_y * sin_z + cos_x * cos_z,
                           sin_x * sin_y * cos_z - cos_x * sin_z},
                          {-sin_y, cos_y * sin_z, cos_y * cos_z}};

    double t_x = 0;
    double t_y = 0;
    double t_z = 0;
    for (std::size_t i = 0; i < points.size(); ++i) {
      t_x = points[i].x;
      t_y = points[i].y;
      t_z = points[i].z;
      points[i].x = rot_m[0][0] * t_x + rot_m[0][1] * t_y + rot_m[0][2] * t_z;
      points[i].y = rot_m[1][0] * t_x + rot_m[1][1] * t_y + rot_m[1][2] * t_z;
      points[i].z = rot_m[2][0] * t_x + rot_m[2][1] * t_y + rot_m[2][2] * t_z;
    }
  }

private:
  std::vector<Point> points_;
  std::vector<Point> inv_points_;
  std::vector<Point> hull_;

  const double k_inf_ = 1e9;
  Point k_inf_point_ = {k_inf_, k_inf_, k_inf_, nullptr, nullptr};
  Point *k_pointer_inf_point_ = &k_inf_point_;

  enum times_enum { LEFT_REC, RIGHT_REC, U_UN_V, UP_U_V, U_VP_V, U_V_VN };

  // To find top of the hull
  void InverseZCoord(std::vector<Point> &points) {
    for (std::size_t i = 0; i < points.size(); ++i) {
      points[i].z = -points[i].z;
    }
  }

  // if result < 0 -> clockwise; else -> counterclockwise
  double TurnPoints(Point *left, Point *mid, Point *right) {
    if (left == k_pointer_inf_point_ || mid == k_pointer_inf_point_ ||
        right == k_pointer_inf_point_) {
      return 1.0;
    }

    return ((mid->x - left->x) * (right->y - left->y) -
            (right->x - left->x) * (mid->y - left->y));
  }

  // Time when clockwise -> counterclockwise && <-
  double TimeWhenTurnPointsChanges(Point *left, Point *mid, Point *right) {
    if (left == k_pointer_inf_point_ || mid == k_pointer_inf_point_ ||
        right == k_pointer_inf_point_) {
      return k_inf_;
    }

    return (((mid->x - left->x) * (right->z - left->z) -
             (right->x - left->x) * (mid->z - left->z)) /
            TurnPoints(left, mid, right));
  }

  Point *MergeSort(Point *points, int size) {
    if (size == 1) {
      points[0].next = k_pointer_inf_point_;
      return points;
    }

    Point *points_first_half = MergeSort(points, (size >> 1));
    Point *points_second_half =
        MergeSort(points + (size >> 1), size - (size >> 1));

    Point head = k_inf_point_;
    Point *curr = &head;
    do {
      if (points_first_half->x < points_second_half->x) {
        curr->next = points_first_half;
        curr = points_first_half;
        points_first_half = points_first_half->next;
      } else {
        curr->next = points_second_half;
        curr = points_second_half;
        points_second_half = points_second_half->next;
      }
    } while (curr != k_pointer_inf_point_);

    return head.next;
  }

  void CreateHull() {
    std::vector<Point> copy_points = points_;

    RotatePoints(points_, kRotateDegree, kRotateDegree, kRotateDegree);
    inv_points_ = points_;
    InverseZCoord(inv_points_);

    int size = static_cast<int>(points_.size());
    std::vector<Point *> base(size << 1);
    std::vector<Point *> inv_base(size << 1);
    std::vector<Point *> tmp(size << 1);

    Hull(MergeSort(&(*points_.begin()), size), size, &(*base.begin()),
         &(*tmp.begin()));
    Hull(MergeSort(&(*inv_points_.begin()), size), size, &(*inv_base.begin()),
         &(*tmp.begin()));

    // Normalize points (from RotatePoints)
    for (std::size_t i = 0; i < points_.size(); ++i) {
      points_[i].x = copy_points[i].x;
      points_[i].y = copy_points[i].y;
      points_[i].z = copy_points[i].z;
      inv_points_[i].x = copy_points[i].x;
      inv_points_[i].y = copy_points[i].y;
      inv_points_[i].z = copy_points[i].z;
    }

    for (int i = 0; base[i] != k_pointer_inf_point_;
         base[i++]->ModifyRelationship()) {
      hull_.push_back(*base[i]);
    }
    for (int i = 0; inv_base[i] != k_pointer_inf_point_;
         inv_base[i++]->ModifyRelationship()) {
      hull_.push_back(*inv_base[i]);
    }
  }

  void Hull(Point *start_point, int size, Point **base, Point **temp) {
    if (size == 1) {
      start_point->next = k_pointer_inf_point_;
      start_point->prev = k_pointer_inf_point_;
      base[0] = k_pointer_inf_point_;
      return;
    }

    Point *u = start_point;
    int half_size = (size >> 1);
    for (int i = 0; i < (half_size - 1); ++i) {
      u = u->next;
    }

    Point *v = u->next;
    Point *mid = v;

    // Recursive on the Left and Right sides of Hull
    Hull(start_point, half_size, temp, base);
    Hull(mid, size - half_size, temp + (half_size << 1),
         base + (half_size << 1));

    // Find start bridge
    for (;;) {
      if (TurnPoints(u, v, v->next) < 0) { // clockwise
        v = v->next;
      } else if (TurnPoints(u->prev, u, v) < 0) { // clockwise
        u = u->prev;
      } else {
        break;
      }
    }

    double times[6];
    double old_time = -k_inf_;
    double new_time = -k_inf_;
    int case_flag = 6;
    int k = 0;
    // Merge Left & Right sides by tracking found bridge UV
    for (int i = 0, j = (half_size << 1);; old_time = new_time) {
      times[LEFT_REC] =
          TimeWhenTurnPointsChanges(temp[i]->prev, temp[i], temp[i]->next);
      times[RIGHT_REC] =
          TimeWhenTurnPointsChanges(temp[j]->prev, temp[j], temp[j]->next);
      times[U_UN_V] = TimeWhenTurnPointsChanges(u, u->next, v);
      times[UP_U_V] = TimeWhenTurnPointsChanges(u->prev, u, v);
      times[U_VP_V] = TimeWhenTurnPointsChanges(u, v->prev, v);
      times[U_V_VN] = TimeWhenTurnPointsChanges(u, v, v->next);

      new_time = k_inf_;
      for (int time_ind = 0; time_ind < 6; ++time_ind) {
        if (times[time_ind] > old_time && times[time_ind] < new_time) {
          case_flag = time_ind;
          new_time = times[time_ind];
        }
      }
      if (new_time == k_inf_) {
        // No changes
        break;
      }

      switch (case_flag) {
      case LEFT_REC:
        if (temp[i]->x < u->x) {
          base[k++] = temp[i];
        }
        temp[i++]->ModifyRelationship();
        break;
      case RIGHT_REC:
        if (temp[j]->x > v->x) {
          base[k++] = temp[j];
        }
        temp[j++]->ModifyRelationship();
        break;
      case U_UN_V:
        u = u->next;
        base[k++] = u;
        break;
      case UP_U_V:
        base[k++] = u;
        u = u->prev;
        break;
      case U_VP_V:
        v = v->prev;
        base[k++] = v;
        break;
      case U_V_VN:
        base[k++] = v;
        v = v->next;
        break;
      default:
        std::cout << "HULL SWITCH ERROR" << std::endl;
        break;
      }
    }
    base[k] = k_pointer_inf_point_;

    // Going back in time to update pointers
    u->next = v;
    v->prev = u;
    for (k--; k >= 0; --k) {
      if (base[k]->x <= u->x || base[k]->x >= v->x) {
        base[k]->ModifyRelationship();
        if (base[k] == u) {
          u = u->prev;
        } else if (base[k] == v) {
          v = v->next;
        }
      } else {
        u->next = base[k];
        base[k]->prev = u;
        v->prev = base[k];
        base[k]->next = v;
        if (base[k]->x < mid->x) {
          u = base[k];
        } else {
          v = base[k];
        }
      }
    }
  }
};

// to find distances from the specified points to the convex hull
class DistancesHull3D {
  friend class Hull3D;

public:
  DistancesHull3D(const std::vector<Point> &points,
                  const std::vector<Point> &requests)
      : hull_3d_(Hull3D(points)), requests_(requests),
        distances_(requests.size(), k_inf_) {}

  std::vector<double> FindDistancesToHull() {
    std::vector<Point> hull = hull_3d_.FindHull();

    for (std::size_t i = 0; i < hull.size(); ++i) {
      FindDistToFacet(distances_, requests_, hull[i].prev, &hull[i],
                      hull[i].next);
    }

    return distances_;
  }

private:
  const double k_inf_ = 1e9;
  Hull3D hull_3d_;
  std::vector<Point> requests_;
  std::vector<double> distances_;

  // Find distance from points to the hull facet
  void FindDistToFacet(std::vector<double> &distances,
                       std::vector<Point> &requests, Point *facet_point_0,
                       Point *facet_point_1, Point *facet_point_2) {
    double a = facet_point_0->y * (facet_point_1->z - facet_point_2->z) +
               facet_point_1->y * (facet_point_2->z - facet_point_0->z) +
               facet_point_2->y * (facet_point_0->z - facet_point_1->z);
    double b = facet_point_0->z * (facet_point_1->x - facet_point_2->x) +
               facet_point_1->z * (facet_point_2->x - facet_point_0->x) +
               facet_point_2->z * (facet_point_0->x - facet_point_1->x);
    double c = facet_point_0->x * (facet_point_1->y - facet_point_2->y) +
               facet_point_1->x * (facet_point_2->y - facet_point_0->y) +
               facet_point_2->x * (facet_point_0->y - facet_point_1->y);
    double d = -(facet_point_0->x * (facet_point_1->y * facet_point_2->z -
                                     facet_point_2->y * facet_point_1->z) +
                 facet_point_1->x * (facet_point_2->y * facet_point_0->z -
                                     facet_point_0->y * facet_point_2->z) +
                 facet_point_2->x * (facet_point_0->y * facet_point_1->z -
                                     facet_point_1->y * facet_point_0->z));
    for (std::size_t i = 0; i < requests.size(); ++i) {
      double distance = std::abs(
          (a * requests[i].x + b * requests[i].y + c * requests[i].z + d) /
          std::sqrt(a * a + b * b + c * c));
      distances[i] = std::min(distances[i], distance);
    }
  }
};

#endif /* CHAN3D_HPP */