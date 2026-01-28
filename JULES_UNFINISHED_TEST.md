# Notes on Tests

## BasicDoubleton.h
- Found a bug in constructor `BasicDoubleton(size_type d, size_type N0 = -1)`.
- `size_type` is usually unsigned (from `capd::IMatrix::size_type`), so `N0 = -1` becomes `MAX_SIZE`.
- The condition `if (N0 < 0)` is always false for unsigned types.
- This causes huge allocation in `setupDimension(d, N0)`, leading to `std::bad_alloc` (or process kill).
- Workaround used in `tests/capd-ddes-storage-DoubletonInterface.cpp`: explicitly pass `N0` (e.g., 0).
- Suggestion: Fix `BasicDoubleton.h` constructor logic (e.g., use signed type for argument or different sentinel).
