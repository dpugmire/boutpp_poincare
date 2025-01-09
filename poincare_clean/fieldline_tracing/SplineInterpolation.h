#pragma once

#include <vector>
#include <stdexcept>
#include <iostream>

float interpolate2D(
    const std::vector<float>& xarray,
    const std::vector<float>& zarray,
    const std::vector<std::vector<float>>& dxdyp,
    float x, float z);

class SplineInterpolation
{
public:
    // Constructor
    SplineInterpolation(const std::vector<float>& x, const std::vector<float>& y)
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
    float evaluate(float t) const
    {
        if (t < x.front() || t > x.back()) {
            throw std::out_of_range("Interpolation point is outside the range of input data.");
        }

        // Find the interval [x[i], x[i+1]] containing t
        auto it = std::lower_bound(x.begin(), x.end(), t);
        size_t i = std::max(static_cast<size_t>(0), static_cast<size_t>(it - x.begin() - 1));

        float h = t - x[i];
        return a[i] + b[i] * h + c[i] * h * h + d[i] * h * h * h;
    }

private:
    std::vector<float> x, y;       // Input data points
    std::vector<float> a, b, c, d; // Spline coefficients

    // Compute spline coefficients
    void computeCoefficients() {
        size_t n = x.size() - 1;
        std::vector<float> h(n), alpha(n), l(n + 1), mu(n), z(n + 1);

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