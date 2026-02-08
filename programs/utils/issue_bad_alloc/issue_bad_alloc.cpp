#include <iostream>
#include "capd/capdlib.h"
#include "capd/ddes/storage/SharedDoubleton.h"

int main() {
    try {
        std::cout << "Starting SharedDoubleton std::bad_alloc reproduction..." << std::endl;

        typedef capd::ddes::SharedDoubleton<capd::IMatrix> Doubleton;
        typedef Doubleton::VectorType VectorType;

        VectorType x1(2); x1[0]=1.0; x1[1]=2.0;
        Doubleton db1(x1);

        VectorType x2(2); x2[0]=0.5; x2[1]=0.5;
        Doubleton db2(x2);

        std::cout << "Adding db2 to db1..." << std::endl;
        db1.add(db2);

        std::cout << "Addition successful. Result x: " << db1.get_x() << std::endl;

        Doubleton::ScalarType c = 2.0;
        std::cout << "Multiplying db1 by " << c << "..." << std::endl;
        db1.mul(c);

        std::cout << "Multiplication successful. Result x: " << db1.get_x() << std::endl;

        std::cout << "Test PASSED (No bad_alloc)" << std::endl;
        return 0;
    } catch (std::bad_alloc& e) {
        std::cerr << "Caught std::bad_alloc: " << e.what() << std::endl;
        return 1;
    } catch (std::exception& e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Caught unknown exception" << std::endl;
        return 1;
    }
}
