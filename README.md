# B-Splines
C++ implementation of fitting B-splines to an Eigen array of data in 3D. 2D data may also be fit by homogenising the data points.

Note: much more documentation and tidying up required, however usage should be along the lines of:

```cpp
  BSpline spline;
  
  spline.fit(<3xn Eigen matrix>, <degree>, <number control points>);
  
  Eigen::Matrix<double,3,Dynamic> spline_samples = spline.sample(<number of samples>, <degree>);
```
Much of the class if left public if you are working with B-spline coefficients, knot vectors etc.

I hope you find some use from this.


## Dependencies
- Eigen (>=v3.4.0)
