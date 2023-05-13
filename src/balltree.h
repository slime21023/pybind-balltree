#ifndef BALL_TREE_HPP
#define BALL_TREE_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <utility>
#include <vector>

using namespace std;
typedef vector<double> Point;

class BallTree {
private:
  double radius;
  int dim;
  int leaf_size;
  Point centroid;
  BallTree *left;
  BallTree *right;
  void split();
  void get_search_leaves(vector<pair<double, BallTree *>> &inserts,
                         const Point &q);

  vector<pair<double, BallTree *>> get_search_leaves(const Point &q,
                                                     const int k);

public:
  vector<Point> points;

  BallTree(int n_dim, int n_leaf, vector<Point> ps);
  ~BallTree();
  vector<Point> kneighbors(const Point &q, const int k);
};

Point calculate_centroid(const vector<Point> &ps, int n_dim) {
  Point p(n_dim, 0);

  for (auto point : ps) {
    for (int i = 0; i < n_dim; i++) {
      p[i] += point[i];
    }
  }

  for (int i = 0; i < n_dim; i++) {
    p[i] /= ps.size();
  }
  return p;
}

double distance(const Point &a, const Point &b, int dim) {
  double sum = 0.0;
  for (int i = 0; i < dim; i++) {
    sum += pow(a[i] - b[i], 2);
  }
  return sqrt(sum);
}

double calculate_radius(Point &centroid, vector<Point> &ps, int dim) {
  double radius = 0.0;
  for (auto &p : ps) {
    double d = distance(p, centroid, dim);
    if (d > radius)
      radius = d;
  }
  return radius;
}

BallTree::BallTree(int n_dim, int n_leaf, vector<Point> ps) {
  dim = n_dim;
  leaf_size = n_leaf;
  left = nullptr;
  right = nullptr;
  points = ps;

  centroid = calculate_centroid(points, dim);
  radius = calculate_radius(centroid, points, dim);

  if (points.size() > leaf_size) {
    split();
  }
}

BallTree::~BallTree() {
  delete left;
  delete right;
}

void BallTree::split() {
  // Select the farest node from the centroid
  Point first_node = points[0];
  double max_d = distance(centroid, first_node, dim);
  for (auto &p : points) {
    double d = distance(centroid, p, dim);
    if (d > max_d) {
      max_d = d;
      first_node = p;
    }
  }

  // Select the farest node from the first node
  Point second_node = points[0];
  max_d = 0.0;
  for (auto &p : points) {
    double d = distance(first_node, p, dim);
    if (d > max_d) {
      max_d = d;
      second_node = p;
    }
  }

  vector<Point> left_points, right_points;
  for (auto &p : points) {
    double fd = distance(p, first_node, dim);
    double sd = distance(p, second_node, dim);

    if (sd > fd) {
      left_points.push_back(p);
    } else {
      right_points.push_back(p);
    }
  }

  left = new BallTree(dim, leaf_size, left_points);
  right = new BallTree(dim, leaf_size, right_points);
  points.clear();
}

void BallTree::get_search_leaves(vector<pair<double, BallTree *>> &inserts,
                                 const Point &q) {
  bool is_leaf = (left == nullptr) && (right == nullptr);
  if (is_leaf) {
    double d = distance(q, centroid, dim);
    inserts.push_back(pair<double, BallTree *>(d, this));
  } else {
    left->get_search_leaves(inserts, q);
    right->get_search_leaves(inserts, q);
  }
}

vector<pair<double, BallTree *>> BallTree::get_search_leaves(const Point &q,
                                                             const int k) {
  vector<pair<double, BallTree *>> leaves;
  bool is_leaf = (left == nullptr) && (right == nullptr);

  if (is_leaf) {
    double d = distance(q, centroid, dim);
    leaves.push_back(pair<double, BallTree *>(d, this));
  } else {
    left->get_search_leaves(leaves, q);
    right->get_search_leaves(leaves, q);
  }

  sort(leaves.begin(), leaves.end(),
       [](pair<double, BallTree *> &a, pair<double, BallTree *> &b) {
         return a.first < b.first;
       });

  return leaves;
}

void query_list_insert(list<pair<double, Point>> &ql, const Point &q, Point p,
                       int dim) {
  double d = distance(q, p, dim);
  if (ql.size() == 0) {
    ql.push_back(pair<double, Point>(d, p));
    return;
  }

  bool has_insert = false;

  for (auto it = ql.begin(); it != ql.end(); it++) {
    double compare = it->first;
    if (d < compare) {
      ql.insert(it, pair<double, Point>(d, p));
      has_insert = true;
      break;
    }
  }

  if (!has_insert) {
    ql.insert(ql.end(), pair<double, Point>(d, p));
  }
}

vector<Point> BallTree::kneighbors(const Point &q, const int k) {
  list<pair<double, Point>> finding;
  vector<pair<double, BallTree *>> leaves = get_search_leaves(q, k);

  auto end = leaves.end();
  if (leaves.size() > k) {
    end = leaves.begin();
    end += k;
  }

  for (auto it = leaves.begin(); it != end; it++) {
    BallTree *tree = it->second;
    for (auto &point : tree->points) {
      query_list_insert(finding, q, point, dim);
    }
  }

  vector<Point> result;
  for (auto it = finding.begin(); it != finding.end(); it++) {
    result.push_back(it->second);
    if (result.size() == k) {
      break;
    }
  }
  return result;
}

#endif