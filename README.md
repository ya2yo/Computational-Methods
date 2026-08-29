# Computational Methods

A modern C++ implementation of algorithms from a computational methods course. The project uses standard CMake and supports Linux, macOS, and Windows with MSVC or MinGW.

## Topics

- Linear systems: Gaussian elimination, LU decomposition, Cholesky decomposition, and the Thomas algorithm
- Iterative linear solvers: Jacobi, Gauss-Seidel, and SOR
- Nonlinear equations: fixed-point iteration, Newton's method, and Steffensen acceleration
- Interpolation and fitting: Lagrange, Newton, Hermite, cubic spline, and least squares
- Numerical integration: Romberg integration
- Ordinary differential equations: fourth-order Runge-Kutta

The project follows a separated layout: public headers are in `include/computational_methods/`, implementations are in `src/`, the demonstration program is in `examples/`, and automated tests are in `tests/`. File names use lowercase snake_case consistently.

## Build

Requirements: CMake 3.20 or newer, a C++17-compatible compiler, and Ninja (recommended).

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

Run the demonstration program:

```bash
./build/debug/computational_methods_demo
```

To build without tests:

```bash
cmake -S . -B build -DCOMPUTATIONAL_METHODS_BUILD_TESTS=OFF
cmake --build build
```

## Development conventions

- Use UTF-8 text and C++17.
- Keep numerical algorithms separate from command-line demonstrations.
- Run the build and tests before committing changes.
- To enable GCC/Clang coverage instrumentation, configure with `-DCOMPUTATIONAL_METHODS_ENABLE_COVERAGE=ON`, then run `ctest` and use `gcovr` or `lcov` to generate a report.
- Do not commit build directories, IDE settings, or generated binaries.
