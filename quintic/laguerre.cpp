#include <vector>
#include <iostream>
#include <iomanip>
#include <complex>
#include <algorithm>

std::pair<std::vector<std::complex<double>>, std::complex<double>> horner(const std::vector<std::complex<double>>& a, std::complex<double> x0) {
	int n = a.size();
	auto b = std::vector<std::complex<double>>(std::max(1, n - 1));

	for (int i = n - 1; i > 0; i--)
		b[i - 1] = a[i] + (i < n - 1 ? b[i] * x0 : 0);
	return std::make_pair(b, a[0] + b[0] * x0);
}

std::complex<double> eval(const std::vector<std::complex<double>>& p, std::complex<double> x) {
	return horner(p, x).second;
}

std::vector<std::complex<double>> derivative(const std::vector<std::complex<double>>& p) {
	int n = p.size();
	auto r = std::vector<std::complex<double>>(std::max(1, n - 1));
	for (int i = 1; i < n; i++)
		r[i - 1] = p[i] * std::complex<double>(i);
	return r;
}

const double EPS = 1e-9;

int cmp(std::complex<double> x, std::complex<double> y) {
	double diff = std::abs(x) - std::abs(y);
	return diff < -EPS ? -1 : (diff > EPS ? 1 : 0);
}

std::complex<double> find_one_root(const std::vector<std::complex<double>>& p0, std::complex<double> x) {
	int n = p0.size() - 1;
	auto p1 = derivative(p0);
	auto p2 = derivative(p1);
	for (int step = 0; step < 10000; step++) {
		auto y0 = eval(p0, x);
		if (cmp(y0, 0) == 0) break;
		auto G = eval(p1, x) / y0;
		auto H = G * G - eval(p2, x) - y0;
		auto R = std::sqrt(std::complex<double>(n - 1) * (H * std::complex<double>(n) - G * G));
		auto D1 = G + R;
		auto D2 = G - R;
		auto a = std::complex<double>(n) / (cmp(D1, D2) > 0 ? D1 : D2);
		x -= a;
		if (cmp(a, 0) == 0) break;
	}
	return x;
}

std::vector<std::complex<double>> find_all_roots(const std::vector<std::complex<double>>& p) {
	std::vector<std::complex<double>> res;
	std::vector<std::complex<double>> q = p;
	while (q.size() > 2) {
		std::complex<double> z(rand() / double(RAND_MAX), rand() / double(RAND_MAX));
		z = find_one_root(q, z);
		z = find_one_root(p, z);
		q = horner(q, z).first;
		res.push_back(z);
	}
	res.push_back(-q[0] / q[1]);
	return res;
}

