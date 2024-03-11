// c++ -O3 -Wall -shared -std=c++17 -fPIC $(python3 -m pybind11 --includes) example.cpp -o example$(python3-config --extension-suffix)

#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int i, int j)
{
    return i + j;
}

PYBIND11_MODULE(example, m)
{
    m.doc() = "pybind11 example"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers");
}