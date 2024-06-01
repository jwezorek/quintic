#include <iostream>
#include <complex>
#include <tuple>
#include <limits>
#include <numeric>
#include <algorithm> 
#include <vector>
#include <chrono>
#include <random>
#include "quintic.hpp"


std::vector<std::complex<double>> find_all_roots(const std::vector<std::complex<double>>& p);

struct Results {
	bool new_is_faster;
	bool new_is_more_accurate;
	int new_fails;
	int old_fails;
	double new_error;
	double old_error;
	double time_taken_by_new;
	double time_taken_by_old;
};

std::random_device rd;

std::vector<std::complex<double>> getRandomQuintic(int lower, int upper)
{
	std::mt19937 e2(rd());
	std::uniform_real_distribution<> dist(lower, upper);
	std::vector<std::complex<double>> quintic(5);
	std::generate(quintic.begin(), quintic.end(), [&]() {return dist(e2); });
	quintic.push_back(1.0);
	return quintic;
}

Results DoTest(int lower, int upper) {
	std::chrono::high_resolution_clock::time_point t1, t2;
	auto coefs = getRandomQuintic(lower, upper);

	quin::details::GeneralQuintic<long double> gp(coefs[4], coefs[3], coefs[2], coefs[1], coefs[0]);
	t1 = std::chrono::high_resolution_clock::now();
	auto [a, b, c, d, e] = gp.getParameters();
	auto [root1, root2, root3, root4, root5] = quin::solveQuintic(a,b,c,d,e);
	t2 = std::chrono::high_resolution_clock::now();
	double new_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
	double new_error = 0.0;
	int new_fails = 0;
	for (auto root : { root1, root2, root3, root4, root5 }) {
		auto err = std::abs(gp.evaluate(root));
		new_error += err;
		new_fails += (err > 0.005) ? 1 : 0;
	}

	t1 = std::chrono::high_resolution_clock::now();
	auto roots = find_all_roots(coefs);
	t2 = std::chrono::high_resolution_clock::now();
	double old_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
	double old_error = 0.0;
	int old_fails = 0;
	for (auto root : roots) {
		auto err = std::abs(gp.evaluate(root));
		old_error += std::abs(gp.evaluate(root));
		old_fails += (err > 0.005) ? 1 : 0;
	}

	return { new_time < old_time, new_error < old_error, new_fails, old_fails, new_error, old_error, new_time, old_time };
}

int main()
{
	int n = 1000;
	int range = 10;

	int new_is_faster = 0;
	int new_is_more_accurate = 0;
	int new_fails = 0;
	int old_fails = 0;
	double new_error = 0.0;
	double old_error = 0.0;
	double time_taken_by_new = 0.0;
	double time_taken_by_old = 0.0;

	for (int i = 0; i < n; i++) {
		auto results = DoTest(-range, range);
		new_is_faster += (results.new_is_faster) ? 1 : 0;
		new_is_more_accurate += (results.new_is_more_accurate) ? 1 : 0;
		new_fails += results.new_fails;
		old_fails += results.old_fails;
		new_error += results.new_error;
		old_error += results.old_error;
		time_taken_by_new += results.time_taken_by_new;
		time_taken_by_old += results.time_taken_by_old;
		if (i % 100 == 0)
			std::cout << i << std::endl;
	}

	std::cout << "speed: " << 100.0 * static_cast<double>(new_is_faster) / static_cast<double>(n) << std::endl;
	std::cout << "accuracy: " << 100.0 * static_cast<double>(new_is_more_accurate) / static_cast<double>(n) << std::endl;
	std::cout << "new fails: " << new_fails << std::endl;
	std::cout << "old fails " << old_fails << std::endl;
	std::cout << "new error: " << new_error << std::endl;
	std::cout << "old error: " << old_error << std::endl;
	std::cout << "total time new: " << time_taken_by_new << std::endl;
	std::cout << "total time old: " << time_taken_by_old << std::endl;
}
