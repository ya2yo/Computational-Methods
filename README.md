# Computational Methods

计算方法课程算法的现代 C++ 实现与实验代码。项目使用标准 CMake 管理构建，可在 Linux、macOS 和 Windows（MSVC/MinGW）上使用。

## 内容

- 线性方程组：高斯消元、LU、Cholesky、追赶法
- 线性方程组迭代：Jacobi、Gauss-Seidel、SOR
- 非线性方程：简单迭代、Newton、Steffensen
- 插值与拟合：Lagrange、Newton、Hermite、三次样条、最小二乘
- 数值积分：Romberg 积分
- 常微分方程：四阶 Runge-Kutta

项目目录遵循头文件与源文件分离原则：公共头文件位于 `include/computational_methods/`，算法实现位于 `src/`，演示程序位于 `examples/`，自动化测试位于 `tests/`。现有文件名暂保留，以避免破坏既有学习笔记。

## 构建

需要 CMake 3.20+、支持 C++17 的编译器；推荐使用 Ninja。

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

运行演示：

```bash
./build/debug/computational_methods_demo
```

不需要测试时可关闭：

```bash
cmake -S . -B build -DCOMPUTATIONAL_METHODS_BUILD_TESTS=OFF
cmake --build build
```

## 开发约定

- 源文件使用 UTF-8 编码，C++ 标准为 C++17。
- 算法实现与命令行演示分离；库目标为 `computational_methods`。
- 提交前执行构建和 `ctest --test-dir build/debug --output-on-failure`。
- 不提交 `build/`、IDE 配置和编译产物。
