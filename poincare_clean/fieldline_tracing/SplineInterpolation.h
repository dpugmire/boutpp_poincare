#pragma once

#include <vector>
#include <stdexcept>
#include <iostream>

double interpolate2D(
    const std::vector<double>& xarray,
    const std::vector<double>& zarray,
    const std::vector<std::vector<double>>& dxdyp,
    double x, double z);

class SplineInterpolation_OLD
{
public:
    // Constructor
    SplineInterpolation_OLD(const std::vector<double>& x, const std::vector<double>& y)
    {
        if (x.size() < 2 || y.size() < 2)
        {
            throw std::invalid_argument("Input vectors must contain at least two points.");
        }

        this->x = x;
        this->y = y;
        computeCoefficients();
    }

    // Evaluate the spline at a given point
    double evaluate(double t) const
    {
        if (t < x.front() || t > x.back()) {
            throw std::out_of_range("Interpolation point is outside the range of input data.");
        }

        // Find the interval [x[i], x[i+1]] containing t
        auto it = std::lower_bound(x.begin(), x.end(), t);
        size_t i = std::max(static_cast<size_t>(0), static_cast<size_t>(it - x.begin() - 1));

        double h = t - x[i];
        return a[i] + b[i] * h + c[i] * h * h + d[i] * h * h * h;
    }

private:
    std::vector<double> x, y;       // Input data points
    std::vector<double> a, b, c, d; // Spline coefficients

    // Compute spline coefficients
    void computeCoefficients() {
        size_t n = x.size() - 1;
        std::vector<double> h(n), alpha(n), l(n + 1), mu(n), z(n + 1);

        a = y;
        for (size_t i = 0; i < n; ++i) {
            h[i] = x[i + 1] - x[i];
        }

        for (size_t i = 1; i < n; ++i) {
            alpha[i] = (3.0f / h[i] * (a[i + 1] - a[i])) - (3.0f / h[i - 1] * (a[i] - a[i - 1]));
        }

        l[0] = 1.0f;
        mu[0] = 0.0f;
        z[0] = 0.0f;

        for (size_t i = 1; i < n; ++i) {
            l[i] = 2.0f * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[n] = 1.0f;
        z[n] = 0.0f;
        c.resize(n + 1, 0.0f);
        b.resize(n, 0.0f);
        d.resize(n, 0.0f);

        for (size_t j = n - 1; j < n; --j) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0f * c[j]) / 3.0f;
            d[j] = (c[j + 1] - c[j]) / (3.0f * h[j]);
        }
    }
};

class SplineInterpolation
{
public:
    // Constructor for 1D interpolation
    SplineInterpolation(const std::vector<double>& x, const std::vector<double>& y);

    // Constructor for 2D interpolation
    SplineInterpolation(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& z);

    // Evaluate 1D interpolation
    double evaluate(double x) const;

    // Evaluate 2D interpolation
    double evaluate(double x, double y) const;

private:
    // 1D data
    std::vector<double> x1D, y1D;
    std::vector<double> a1D, b1D, c1D, d1D;

    // 2D data
    std::vector<double> x2D, y2D;
    std::vector<std::vector<double>> z2D;

    // Helper for 1D spline coefficients
    void compute1DSplineCoefficients();

    // Helper for 1D evaluation
    size_t findInterval(const std::vector<double>& x, double val) const;

    // 2D spline coefficients
    std::vector<std::vector<double>> a2D, b2D, c2D, d2D;

    // Helper for 2D spline coefficients
    void compute2DSplineCoefficients();
};
