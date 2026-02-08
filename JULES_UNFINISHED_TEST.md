# Notes on Tests

## BasicDoubleton.h (DONE)
- Tests implemented in `tests/capd-ddes-storage-BasicDoubleton.cpp`. Coverage > 94%.
- **Known Bug 1:** Constructor `BasicDoubleton(size_type d, size_type N0 = -1)` fails due to unsigned underflow of `N0 = -1`.
- **Known Bug 2:** Constructor `BasicDoubleton(VectorType* x, ...)` fails to compile due to passing pointers instead of references to `setupFromData` and ambiguous `VectorType(0)`.
- These bugs are marked in the header file and corresponding tests are disabled or avoided.

## SharedDoubleton.h (UNFINISHED)
- Tests implemented in `tests/capd-ddes-storage-SharedDoubleton.cpp`. Coverage ~27.6% (low).
- **Known Bug 1:** `SharedDoubleton::add(SharedDoubleton const &)` fails with `std::bad_alloc` for small test vectors (size 2), likely due to memory corruption or huge allocation in `capd::vectalg::split` or internal logic.
- **Test Disabled:** `AddSetTest` is commented out/disabled. `MulThenAddTest` is also disabled. `MulTest` catches `std::bad_alloc` but warns.
- **Verification:** The bug was verified to persist even after a clean rebuild with coverage enabled.
- Future work: Debug memory issue in `add`. Improve test coverage for complex ownership logic (`rawSetup`, `assureOwner`, `sanityCheck`, `deallocate`).
