#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "./balltree.h"

namespace py = pybind11;

struct Tree
{
    Tree(const py::array_t<double> points, int leaf_num)
    {
        py::buffer_info buf = points.request();

        int rows = buf.shape[0];
        int cols = buf.shape[1];

        vector<Point> ps;
        for (int idx = 0; idx < rows; idx++)
        {
            double *row = static_cast<double *>(buf.ptr) + idx * cols;
            Point p(row, row + cols);

            ps.push_back(p);
        }

        t = new BallTree(cols, leaf_num, ps);
    };

    vector<Point> kneighbors(const Point &q, const int k)
    {
        return t->kneighbors(q, k);
    };

    BallTree *t;
};    

PYBIND11_MODULE(pybind_balltree, m)
{
    py::class_<Tree>(m, "BallTree")
        .def(py::init<const py::array_t<double> &, int>())
        .def("kneighbors", &Tree::kneighbors);
}