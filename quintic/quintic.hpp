#pragma once

#include <iostream>
#include <complex>
#include <tuple>
#include <limits>
#include <numeric>
#include <algorithm> 
#include <vector>

using namespace std::complex_literals;

namespace quin {

	template<typename F>
	using Z = std::complex<F>;

	template<typename F>
	const F EPS = F(1e-7);

	template<typename F>
	struct options {
		F epsilon;
		int max_iterations;
		int max_refinement_iterations;
		F root_refinement_thresh;

		options(F eps = EPS<F>, int iter = 10000, int r_iter = 1000, F r_thresh = 0.00005) :
			epsilon(eps),
			max_iterations(iter),
			max_refinement_iterations(r_iter),
			root_refinement_thresh(r_thresh)
		{}
	};

	namespace details {

		template<typename F>
		const F EPS = 1e-7;

		template<typename F>
		const F nan = std::numeric_limits<F>::quiet_NaN();

		template<typename F>
		Z<F> operator*(int c, Z<F> z) {
			return Z<F>(c) * z;
		}

		template<typename F>
		Z<F> operator-(Z<F> z, int a) {
			return z - Z<F>(a);
		}

		template<typename F>
		Z<F> operator-(int a, Z<F> z) {
			return Z<F>(a) - z;
		}

		template<typename F>
		Z<F> operator+(Z<F> z, int a) {
			return z + Z<F>(a);
		}

		template<typename F>
		Z<F> operator/(Z<F> z, int a) {
			return z / Z<F>(a);
		}

		template<typename F>
		Z<F> operator^(Z<F> z, int exponent) {
			Z<F> val = F(1);
			for (int i = 0; i < exponent; i++)
				val *= z;
			return val;
		}

		template <typename T> int sgn(T val) {
			return (T(0) < val) - (val < T(0));
		}

		template <typename T,typename F>
		bool is_zero(T z, F epsilon = EPS) {
			return std::abs(z) < epsilon;
		}

		template <typename F>
		bool is_real(Z<F> z, F epsilon = EPS)
		{
			return is_zero(z.imag(), epsilon);
		}

		template <typename F>
		bool is_nan(Z<F> z)
		{
			return isnan(z.real()) || isnan(z.imag());
		}

		template<typename T>
		T const pi = std::acos(-T(1));

		template<typename F>
		std::tuple<Z<F>, Z<F>> performDoyleMcMullenIteration(Z<F> C, const options<F>& opts) {

			auto g = [](Z<F> z, Z<F> w) -> Z<F> {
				return 91125 * (z ^ 6) + ((-133650) * (w ^ 2) + 61560 * w - 193536) * (z ^ 5)
					+ ((-66825) * (w ^ 4) + 142560 * (w ^ 3) + 133056 * (w ^ 2) + (-61440) * w + 102400) * (z ^ 4)
					+ (5940 * (w ^ 6) + 4752 * (w ^ 5) + 63360 * (w ^ 4) + (-140800) * (w ^ 3)) * (z ^ 3)
					+ ((-1485) * (w ^ 8) + 3168 * (w ^ 7) + (-10560) * (w ^ 6)) * (z ^ 2)
					+ (440 * (w ^ 9) - 66 * (w ^ 10)) * z
					+ (w ^ 12);
			};

			auto h = [](Z<F> z, Z<F> w) -> Z<F> {
				return ((1215 * w - 648) * (z ^ 4)
					+ ((-540) * (w ^ 3) - 216 * (w ^ 2) + (-1152) * w + 640) * (z ^ 3)
					+ (378 * (w ^ 5) - 504 * (w ^ 4) + 960 * (w ^ 3)) * (z ^ 2)
					+ (36 * (w ^ 7) - 168 * (w ^ 6)) * z
					- (w ^ 9));
			};

			auto  g_prime = [](Z<F> z, Z<F> w) -> Z<F> {
				return 12 * (w ^ 11) - 1620 * (165 * w - 38) * (z ^ 5)
					- 12 * (22275 * (w ^ 3) - 35640 * (w ^ 2) - 22176 * w + 5120) * (z ^ 4)
					+ 1320 * (27 * (w ^ 5) + 18 * (w ^ 4) + 192 * (w ^ 3) - 320 * (w ^ 2)) * (z ^ 3)
					- 792 * (15 * (w ^ 7) - 28 * (w ^ 6) + 80 * (w ^ 5)) * (z ^ 2)
					- 660 * ((w ^ 9) - 6 * (w ^ 8)) * z;
			};

			auto z = 1 - 1728 * C;
			auto T = [g, g_prime, z](Z<F> w) {
				return w - 12 * (g(z, w) / g_prime(z, w));
			};
			auto u = [h, g, z](Z<F> wi) {
				return (100 * z * (z - 1) * h(z, wi)) / g(z, wi);
			};

			int count = 0;
			Z<F> w(1), old_old_w, old_w, t;
			do {
				old_old_w = old_w;
				old_w = w;
				t = T(w);
				if (is_nan(t)) {
					// This occurs if g_prime(w) == 0 causing a division by zero.
					// We will just bail out with the last good w and T(w) pair we have.
					w = old_old_w;
					break;
				}
				w = T(t);
				if (is_nan(w)) {
					w = old_w;
					break;
				}
				count++;
			} while (std::abs(w - old_w) > opts.epsilon && count < opts.max_iterations);

			auto u1 = u(w), u2 = u(T(w));
			return {
				(Z<F>(9, std::sqrt(F(15))) * u1 + Z<F>(9, -std::sqrt(F(15))) * u2) / 90,
				(Z<F>(9, std::sqrt(F(15))) * u2 + Z<F>(9, -std::sqrt(F(15))) * u1) / 90
			};
		}

		template<typename F>
		std::tuple<std::vector<Z<F>>, std::vector<Z<F>>> performSyntheticDivision(
			const std::vector<Z<F>>& dividend,
			const std::vector<Z<F>>& divisor, F epsilon = EPS)
		{
			std::vector<Z<F>> output(dividend);
			auto normalizer = divisor[0];

			for (size_t i = 0; i < dividend.size() - (divisor.size() - 1); i++) {
				output[i] = output[i] / normalizer;
				auto coef = output[i];
				if (std::abs(coef) > epsilon) {
					for (size_t j = 1; j < divisor.size(); j++)
						output[i + j] += -divisor[j] * coef;
				}
			}

			auto separator = output.size() - (divisor.size() - 1);
			return {
				std::vector(output.begin(), output.begin() + separator),
				std::vector(output.begin() + separator, output.end())
			};
		}

		template<typename F>
		std::vector<Z<F>> multiplySimpleBinomials(Z<F> z1, Z<F> z2)
		{
			return { F(1), z1 + z2, z1 * z2 };
		}

		template<typename F>
		std::tuple<Z<F>, Z<F>, Z<F>> solveCubic(Z<F> a, Z<F> b, Z<F> c, Z<F> d, F epsilon = EPS)
		{
			b /= a;
			c /= a;
			d /= a;

			auto Q = (b * b - 3 * c) / 9;
			auto R = (2 * b * b * b - 9 * b * c + 27 * d) / 54;

			if (is_real(Q, epsilon) && is_real(R, epsilon)) {
				auto r = R.real();
				auto q = Q.real();
				if (r * r < q * q * q) {
					// this equation has 3 real roots.
					auto theta = std::acos(r / std::sqrt(q * q * q));
					return {
						-F(2) * std::sqrt(q) * std::cos(theta / F(3)) - b.real() / F(3),
						-F(2) * std::sqrt(q) * std::cos((theta + F(2) * pi<F>) / F(3)) - b.real() / F(3),
						-F(2) * std::sqrt(q) * std::cos((theta - F(2) * pi<F>) / F(3)) - b.real() / F(3)
					};
				}
			}

			auto sqrt_term = std::sqrt(R * R - Q * Q * Q);
			auto sign = sgn((std::conj(R) * sqrt_term).real());
			auto A = -std::pow(R + sign * sqrt_term, F(1) / F(3));
			auto B = (!is_zero(A, epsilon)) ? Q / A : 0.0;

			return {
				(A + B) - (b / 3),
				-Z<F>(0.5) * (A + B) - (b / 3) + Z<F>(1i) * (std::sqrt(F(3)) / F(2)) * (A - B),
				-Z<F>(0.5) * (A + B) - (b / 3) - Z<F>(1i) * (std::sqrt(F(3)) / F(2)) * (A - B)
			};
		}

		template<typename F>
		std::tuple<Z<F>, Z<F>> solveQuadratic(Z<F> a, Z<F> b, Z<F> c)
		{
			auto determinant = b * b - 4 * a * c;
			return {
				(-b + std::sqrt(determinant)) / (2 * a),
				(-b - std::sqrt(determinant)) / (2 * a)
			};
		}

		template<typename F>
		class GeneralQuintic
		{
		private:
			Z<F> a, b, c, d, e;
		public:
			GeneralQuintic(Z<F> a, Z<F> b, Z<F> c, Z<F> d, Z<F> e) :
				a(a), b(b), c(c), d(d), e(e)
			{}

			std::tuple<Z<F>, Z<F>, Z<F>, Z<F>, Z<F>> getParameters() const
			{
				return { a,b,c,d,e };
			}

			Z<F> evaluate(Z<F> x) const {
				return (x ^ 5) + a * (x ^ 4) + b * (x ^ 3) + c * (x ^ 2) + d * x + e;
			}

			Z<F> evalDerivative(Z<F> x) const {
				return 5 * (x ^ 4) + 4 * a * (x ^ 3) + 3 * b * (x ^ 2) + 2 * c * x + d;
			}

			Z<F> evalSecondDerivative(Z<F> x) const {
				return 20 * x * x * x + 12 * a * x * x + 6 * b * x + 2 * c;
			}
		};

		template<typename F>
		class PrincipalForm
		{
		private:
			GeneralQuintic<F> gq;
			Z<F> m, n, c3, c4, c5;

			std::tuple<Z<F>, Z<F>> getGqRootCandidates(Z<F> y) const
			{
				return {
					-F(0.5) * m - F(0.5) * std::sqrt((m ^ 2) - 4 * n + 4 * y),
					-F(0.5) * m + F(0.5) * std::sqrt((m ^ 2) - 4 * n + 4 * y)
				};
			}

		public:
			PrincipalForm(const GeneralQuintic<F>& gq) :
				gq(gq)
			{
				auto [a, b, c, d, e] = gq.getParameters();

				// where a, b, c, d, and e are coefficients of a general quintic:
				// x^5 + ax^4 + bx^3 + cx^2 + dx + e = 0

				m = F(0.5) * (4 * (a ^ 3) - 13 * a * b + 15 * c + std::sqrt(-15 * (a ^ 2) * (b ^ 2) + 60 * (b ^ 3) + 10 * (4 * (a ^ 3) - 19 * a * b) * c
					+ 225 * (c ^ 2) + 40 * (2 * (a ^ 2) - 5 * b) * d)) / (2 * (a ^ 2) - 5 * b);

				n = F(0.1) * (5 * (a ^ 2) * b - 20 * (b ^ 2) + 15 * a * c + std::sqrt(-15 * (a ^ 2) * (b ^ 2) + 60 * (b ^ 3)
					+ 10 * (4 * (a ^ 3) - 19 * a * b) * c + 225 * (c ^ 2) + 40 * (2 * (a ^ 2) - 5 * b) * d) * a) / (2 * (a ^ 2) - 5 * b);

				c3 = c * (m ^ 3) - (a * c - 4 * d) * (m ^ 2) - 6 * ((a ^ 2) - a * m - 2 * b) * (n ^ 2) - 10 * (n ^ 3) - (c ^ 2) + 2 * b * d - 2 * a * e
					+ (b * c - 3 * a * d + 5 * e) * m - 3 * (b * (m ^ 2) + (b ^ 2) - 2 * a * c - (a * b - 3 * c) * m + 2 * d) * n;

				c4 = d * (m ^ 4) - (a * d - 5 * e) * (m ^ 3) + 4 * ((a ^ 2) - a * m - 2 * b) * (n ^ 3) + 5 * (n ^ 4) + (b * d - 4 * a * e) * (m ^ 2)
					+ 3 * (b * (m ^ 2) + (b ^ 2) - 2 * a * c - (a * b - 3 * c) * m + 2 * d) * (n ^ 2) + (d ^ 2) - 2 * c * e
					- (c * d - 3 * b * e) * m - 2 * (c * (m ^ 3) - (a * c - 4 * d) * (m ^ 2) - (c ^ 2) + 2 * b * d - 2 * a * e + (b * c - 3 * a * d + 5 * e) * m) * n;

				c5 = -a * e * (m ^ 4) + e * (m ^ 5) + b * e * (m ^ 3) - ((a ^ 2) - a * m - 2 * b) * (n ^ 4) - (n ^ 5)
					- c * e * (m ^ 2) - (b * (m ^ 2) + (b ^ 2) - 2 * a * c - (a * b - 3 * c) * m + 2 * d) * (n ^ 3)
					+ d * e * m + (c * (m ^ 3) - (a * c - 4 * d) * (m ^ 2) - (c ^ 2) + 2 * b * d - 2 * a * e + (b * c - 3 * a * d + 5 * e) * m) * (n ^ 2)
					- (e ^ 2) - (d * (m ^ 4) - (a * d - 5 * e) * (m ^ 3) + (b * d - 4 * a * e) * (m ^ 2) + (d ^ 2) - 2 * c * e - (c * d - 3 * b * e) * m) * n;
			}

			std::tuple<Z<F>, Z<F>, Z<F>> getParameters() const
			{
				return { c3 / 5, c4 / 5, c5 };
			}

			Z<F> convertRootToGeneralQuintic(Z<F> y)
			{
				auto [possible_root_1, possible_root_2] = getGqRootCandidates(y);
				return (std::abs(gq.evaluate(possible_root_1)) < std::abs(gq.evaluate(possible_root_2))) ?
					possible_root_1 :
					possible_root_2;
			}

			Z<F> evaluate(Z<F> y) {
				return (y ^ 5) + c3 * y * y + c4 * y + c5;
			}
		};

		template<typename F>
		class BrioschiForm
		{
		private:
			Z<F> lambda;
			Z<F> p;
			Z<F> q;
			Z<F> j;

		public:
			BrioschiForm(const PrincipalForm<F>& pq)
			{
				auto [a, b, c] = pq.getParameters();
				// a, b, and c are the parameters of a principal quintic:
				// y^5 + 5ay^2 + 5by + c = 0...

				auto [lambda_candidate_1, lambda_candidate_2] = solveQuadratic(
					((a ^ 4) + a * b * c - (b ^ 3)),
					-(11 * (a ^ 3) * b - a * (c ^ 2) + 2 * (b ^ 2) * c),
					(64 * (a ^ 2) * (b ^ 2) - 27 * (a ^ 3) * c - b * (c ^ 2))
				);

				lambda = (std::abs(lambda_candidate_1) < std::abs(lambda_candidate_2)) ? lambda_candidate_1 : lambda_candidate_2;

				j = ((a * (lambda ^ 2) - 3 * b * lambda - 3 * c) ^ 3) / ((a ^ 2) * (a * c * lambda - (b ^ 2) * lambda - b * c));
				p = (-8 * a * (lambda ^ 3) - 72 * b * (lambda ^ 2) - 72 * c * lambda + j * (a ^ 2)) / (a * (lambda ^ 2) + b * lambda + c);
				q = F(1.0) / (1728 - j);
			}

			Z<F> getParameter() const {
				return q;
			}

			Z<F> convertRootToPrincipal(Z<F> w) {
				return (p * w + lambda) / ((w * w) / q - 3);
			}

			Z<F> evaluate(Z<F> w) const {
				return (w ^ 5) - 10 * (q ^ 3) + 45 * q * q * w - q * q;
			}
		};

		template<typename F>
		std::tuple<Z<F>, Z<F>, Z<F>, Z<F>> solveQuartic(Z<F> G, Z<F> a, Z<F> b, Z<F> c, Z<F> d) {
			a /= G;
			b /= G;
			c /= G;
			d /= G;

			auto one_third = (F(1) / F(3));
			auto [v1, v2] = solveQuadratic(
				Z<F>(1),
				-2 * (b ^ 3) + 9 * a * b * c - 27 * (c ^ 2) - 27 * (a ^ 2) * d + 72 * b * d,
				((b ^ 2) - 3 * a * c + 12 * d) ^ 3
			);
			auto cube_root_of_v1 = std::pow(v1, one_third);
			auto u = ((3 * a * a - 8 * b) / 12) + one_third * (cube_root_of_v1 + (b * b - 3 * a * c + 12 * d) / cube_root_of_v1);
			auto sqrt_of_u = std::sqrt(u);
			auto radical_1 = std::sqrt(
				3 * a * a - 8 * b - 4 * u + (-a*a*a + 4*a*b - 8*c) / sqrt_of_u
			);
			auto radical_2 = std::sqrt(
				3 * a * a - 8 * b - 4 * u - (-a * a * a + 4 * a * b - 8 * c) / sqrt_of_u
			);
			return {
				-F(0.25) * a + F(0.5) * sqrt_of_u + F(0.25) * radical_1,
				-F(0.25) * a + F(0.5) * sqrt_of_u - F(0.25) * radical_1,
				-F(0.25) * a - F(0.5) * sqrt_of_u + F(0.25) * radical_2,
				-F(0.25) * a - F(0.5) * sqrt_of_u - F(0.25) * radical_2
			};
		}

		template<typename F>
		Z<F> refineRoot(const GeneralQuintic<F>& q, Z<F> root, double thresh, int iter)
		{
			// Halley's method...

			Z<F> z = root;
			int count = 0;
			do {
				auto f = q.evaluate(z);
				auto f_prime = q.evalDerivative(z);
				auto f_prime_prime = q.evalSecondDerivative(z);

				z = z - (2 * f * f_prime) / (2 * f_prime * f_prime - f * f_prime_prime);
				count++;
			} while (std::abs(q.evaluate(z)) > thresh && count < iter);

			return z;
		}

		template<typename F>
		std::tuple<Z<F>, Z<F>> refineRoots(const GeneralQuintic<F>& quintic, Z<F> root1, Z<F> root2, F thresh, int iter) {
			if (std::abs(quintic.evaluate(root1)) > thresh)
				root1 = refineRoot(quintic, root1, thresh, iter);

			if (std::abs(quintic.evaluate(root2)) > thresh)
				root2 = refineRoot(quintic, root2, thresh, iter);

			auto error_root1 = std::abs(quintic.evaluate(root1));
			auto error_root2 = std::abs(quintic.evaluate(root2));

			if (error_root1 < thresh && error_root2 < thresh)
				return { root1, root2 };

			if (error_root1 < thresh)
				return { root1, nan<F> };

			if (error_root2 < thresh)
				return { root2, nan<F> };

			return { root1, root2 };
		}

		template<typename F>
		std::tuple<Z<F>, Z<F>> findOneOrTwoRoots(Z<F> a, Z<F> b, Z<F> c, Z<F> d, Z<F> e, const options<F>& options)
		{
			GeneralQuintic<F> quintic(a, b, c, d, e);
			PrincipalForm<F> principal(quintic);
			BrioschiForm<F> broischi(principal);

			auto [broischi_root1, broischi_root2] = performDoyleMcMullenIteration(broischi.getParameter(), options);

			auto root1 = principal.convertRootToGeneralQuintic(broischi.convertRootToPrincipal(broischi_root1));
			auto root2 = principal.convertRootToGeneralQuintic(broischi.convertRootToPrincipal(broischi_root2));

			auto error_root1 = std::abs(quintic.evaluate(root1));
			auto error_root2 = std::abs(quintic.evaluate(root2));

			if (error_root1 < options.root_refinement_thresh && error_root2 < options.root_refinement_thresh)
				return { root1, root2 };

			if (error_root1 < options.root_refinement_thresh)
				return { root1, nan<F> };

			if (error_root2 < options.root_refinement_thresh)
				return { root2, nan<F> };

			return (options.max_refinement_iterations > 0) ?
				refineRoots(quintic, root1, root2, options.root_refinement_thresh, options.max_refinement_iterations) :
				std::tuple<Z<F>, Z<F>>{root1, root2};
		}

		template<typename F>
		std::tuple<Z<F>, Z<F>, Z<F>, Z<F>, Z<F>> solveQuintic(Z<F> a, Z<F> b, Z<F> c, Z<F> d, Z<F> e, const options<F>& opts)
		{
			auto [root1, root2] = details::findOneOrTwoRoots(a, b, c, d, e, opts);

			if (!is_nan(root2)) {
				auto quadratic = details::multiplySimpleBinomials(-root1, -root2);
				auto quintic = std::vector<Z<F>>{ 1.0, a, b, c, d, e };

				auto [quotient, remainder] = details::performSyntheticDivision(quintic, quadratic, opts.epsilon);
				//	if (remainder.size() != 2 || !is_zero(remainder[0]) || !is_zero(remainder[1]))
				//		return { nan,nan,nan,nan,nan };

				auto [root3, root4, root5] = solveCubic(quotient[0], quotient[1], quotient[2], quotient[3], opts.epsilon);
				return { root1, root2, root3, root4, root5 };
			} else {
				auto quintic = std::vector<Z<F>>{ 1.0, a, b, c, d, e };
				auto linear = std::vector<Z<F>>{ 1.0, -root1 };
				auto [quotient, remainder] = details::performSyntheticDivision(quintic, linear, opts.epsilon);
				auto [root2, root3, root4, root5] = solveQuartic(quotient[0], quotient[1], quotient[2], quotient[3], quotient[4]);
				return { root1, root2, root3, root4, root5 };
			}
		}
	}

	template<typename F>
	std::tuple<Z<F>, Z<F>> solveQuadratic(Z<F> a, Z<F> b, Z<F> c)
	{
		return details::solveQuadratic(a, b, c);
	}

	template<typename F>
	std::tuple<Z<F>, Z<F>, Z<F>> solveCubic(Z<F> a, Z<F> b, Z<F> c, Z<F> d, F epsilon = EPS)
	{
		return details::solveCubic(a, b, c, d, epsilon);
	}

	template<typename F>
	std::tuple<Z<F>, Z<F>, Z<F>, Z<F>> solveQuartic(Z<F> a, Z<F> b, Z<F> c, Z<F> d, Z<F> e)
	{
		return details::solveQuartic(a, b, c, d, e);
	}

	template<typename F>
	std::tuple<Z<F>, Z<F>, Z<F>, Z<F>, Z<F>> solveQuintic(Z<F> a, Z<F> b, Z<F> c, Z<F> d, Z<F> e, const options<F>& opts = options<F>())
	{
		return details::solveQuintic(a, b, c, d, e, opts);
	}
}

