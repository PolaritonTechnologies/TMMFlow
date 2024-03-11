#include <iostream>

extern "C"
{
    int add(int a, int b)
    {
        int result = a + b;
        std::cout << "Adding " << a << " and " << b << ". Result: " << result << std::endl;
        return result;
    }
}