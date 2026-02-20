#include <capd/capdlib.h>
#include <iostream>

int main() {
    capd::DVector v1(2), v2(2);
    v1[0] = 1.0; v1[1] = 2.0;
    v2[0] = 1.0; v2[1] = 2.0;

    // Check if subsetInterior exists
    // bool b = v1.subsetInterior(v2); // This will likely fail to compile if method doesn't exist
    // std::cout << "subsetInterior: " << b << std::endl;
    return 0;
}
