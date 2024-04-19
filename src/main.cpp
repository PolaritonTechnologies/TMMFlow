#include "FilterStack.h"
#include <iostream>
#include <cmath>

int main()
{
    // Create a FilterStack
    FilterStack filter = FilterStack("../src/temp/temp_cpp_order.json");

    std::vector<double> d_list_vector({45.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 100.0, 45.0});
    filter.change_material_thickness(d_list_vector);

    filter.get_material_order();
    filter.get_thicknesses();

    std::vector<double> d_list_vector2({111.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 111.0, 45.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 65.0, 45.0, 100.0, 45.0});
    filter.change_material_thickness(d_list_vector2);

    std::vector<int> order_list({21, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 0, 22});
    filter.change_material_order(order_list);

    filter.get_material_order();
    filter.get_thicknesses();

    // change_material_order(fs, material_order, size_material_order);

    // Destroy the FilterStack when done

    return 0;
}