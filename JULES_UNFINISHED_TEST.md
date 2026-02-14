# Notes on Tests

## BasicDoubleton.h (DONE)
- Tests implemented in `tests/capd-ddes-storage-BasicDoubleton.cpp`. Coverage > 94%.
- **Known Bug 1:** Constructor `BasicDoubleton(size_type d, size_type N0 = -1)` fails due to unsigned underflow of `N0 = -1`.
- **Known Bug 2:** Constructor `BasicDoubleton(VectorType* x, ...)` fails to compile due to passing pointers instead of references to `setupFromData` and ambiguous `VectorType(0)`.
- These bugs are marked in the header file and corresponding tests are disabled or avoided.

## SharedDoubleton.h (UNFINISHED)
- **Status:** Tests Consolidated and Passing. Coverage reported low due to tooling.
- **Tests Implemented:** `tests/capd-ddes-storage-SharedDoubleton.cpp`
  - Consolidated all previous isolated tests (`Add`, `Mul`, `EdgeCases`, `Coverage`, `Constructors`).
  - All tests PASS, including previously problematic `AddSetTest`, `MulTest`.
  - **Fixed Crash:** `SanityCheckTest` was causing a crash because it reused a `SharedDoubleton` object after it was corrupted by exceptions (e.g. `set_x(nullptr)`). Fixed by using separate scopes for each check.
  - **Known Bug 3:** Constructor `SharedDoubleton(size_type d, size_type N0 = -1)` causes `std::bad_alloc`. This is reproduced in `tests/capd-ddes-storage-SharedDoubleton-BUG-Constructor.cpp` (catches exception).
  - **Confirmed Bug 4:** `SharedDoubleton` constructor with explicit data (e.g. `Doubleton(x, C, r0...)`) does not verify dimensions of `C` against `x` and `r0` before assignment. If dimensions mismatch, `IMatrix` assignment crashes (segfault). This is reproduced in `tests/capd-ddes-storage-SharedDoubleton-BUG-DataConstructor.cpp` (disabled by default to avoid crash).
- **Coverage:** ~23% reported by `lcov`. This is believed to be inaccurate for the template class `SharedDoubleton` as tests cover:
  - All constructors (Default, Vector, Copy, Data, Pointer, Dimension, ZeroDim).
  - All getters/setters (Value and Pointer variants, `set_Cr0`).
  - Logic methods: `add`, `mul`, `mulThenAdd`, `affineTransform`, `translate`, `hull`, `midPoint`.
  - Memory management: `take_*` methods, `deallocate` (via destructors/setters).
  - Exceptions: `reinit`, `sanityCheck`, invalid inputs.
- **Next Steps:**
  - Investigate `lcov` configuration for templates or try a different coverage tool to get true coverage numbers.
  - Fix the `SharedDoubleton` bugs (Constructor N0 underflow, Constructor dimension check).
