#ifndef _STATISTICS_TEST_H
#define _STATISTICS_TEST_H

#include <algorithm>
#include <numeric>
#include <cmath>

#define STATS_ENABLE_STDVEC_WRAPPERS
#include "../../../stats/include/stats.hpp"

namespace refresh
{
	// *************************************************************************************
	// 
	// *************************************************************************************
	class statistical_test
	{
	public:
		struct t_test_t
		{
			double statistic;
			double p_value;
			double df;
		};

		struct mann_whitney_t
		{
			double statistic_U1;
			double statistic_U2;
			double p_value;
		};

		struct snr_t
		{
			double ratio;
			double p_value;
		};

	private:
		std::vector<std::pair<double, int>> mem_mann_whitney;

		static double pow2(size_t x)
		{
			return pow2((double)x);
		}

		static double pow2(double x)
		{
			return x * x;
		}

		static double pow3(size_t x)
		{
			return pow3((double)x);
		}

		static double pow3(double x)
		{
			return x * x * x;
		}

		template<typename Iter>
		static double mean(Iter first, size_t n)
		{
			if (n == 0)
				return 0;

			return std::accumulate(first, first + n, (double)0.0) / (double) n;
		}

		template<typename Iter>
		static double std_dev(double x, Iter first, size_t n)
		{
			if (n < 2)
				return 0;

			double res = 0;

			for (size_t i = 0; i < n; ++i)
				res += pow2(*first++ - x);

			return sqrt(res / (double)(n - 1));
		}

	public:
		statistical_test() = default;

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static t_test_t t_test(Iter1 X_first, Iter1 X_last, Iter2 Y_first, bool equal_var = true)
		{
			return t_test_n(X_first, Y_first, std::distance(X_first, X_last), equal_var);
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static t_test_t t_test_n(Iter1 X_first, Iter2 Y_first, size_t n, bool equal_var = true)
		{
			t_test_t ret;

			double m1 = mean(X_first, n);
			double m2 = mean(Y_first, n);
			double sd1 = std_dev(m1, X_first, n);
			double sd2 = std_dev(m2, Y_first, n);

			if (equal_var)
			{
				double sp = sqrt((pow2(sd1) + pow2(sd2)) / 2.0);

				ret.statistic = (m1 - m2) / (sp * sqrt(2.0 / n));
				ret.df = 2.0 * n - 2.0;
			}
			else
			{
				double s_delta = sqrt(pow2(sd1) / (double)n + pow2(sd2) / (double)n);
				ret.statistic = (m1 - m2) / s_delta;

				double a1 = pow2(sd1) / n;
				double a2 = pow2(sd2) / n;

				ret.df = pow2(a1 + a2) / (pow2(a1) / (n - 1) + pow2(a2) / (n - 1));
			}

			ret.p_value = std::clamp<double>(2.0 * (1.0 - stats::pt(ret.statistic, ret.df, false)), 0, 1);

			return ret;
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		mann_whitney_t mann_whitney_U_test(Iter1 X_first, Iter1 X_last, Iter2 Y_first)
		{
			return mann_whitney_U_test_n(X_first, Y_first, std::distance(X_first, X_last));
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		mann_whitney_t mann_whitney_U_test_n(Iter1 X_first, Iter2 Y_first, size_t n)
		{
//			const size_t approx_method_thr = 8;		// as in SciPy
			const size_t approx_method_thr = 0;

			mann_whitney_t ret;

			mem_mann_whitney.clear();
			mem_mann_whitney.resize(2*n);

			for (size_t i = 0; i < n; ++i)
				mem_mann_whitney[i] = std::make_pair(*X_first++, 1);
			for (size_t i = 0; i < n; ++i)
				mem_mann_whitney[n+i] = std::make_pair(*Y_first++, 2);

			std::sort(mem_mann_whitney.begin(), mem_mann_whitney.end());

			double R1 = 0;
			double R2 = 0;

			double tie_corr = 0;

			for (size_t i = 0; i < 2 * n;)
			{
				size_t j;
				for (j = i + 1; j < 2 * n; ++j)
					if (mem_mann_whitney[j].first != mem_mann_whitney[i].first)
						break;

				double adder = 1.0 + (double)(i + j - 1) / 2.0;

				for(size_t k = i; k < j; ++k)
					if (mem_mann_whitney[k].second == 2)
						R2 += adder;
					else
						R1 += adder;

				if (j > i + 1)
				{
					double t = (double) (j - i);
					tie_corr += pow3(t) - (t);
				}

				i = j;
			}

			double V = pow2(n) + n * (n + 1) / 2.0;
			ret.statistic_U1 = V - R2;
			ret.statistic_U2 = V - R1;

			if (n > approx_method_thr)
			{
				double m_u = pow2(n) / 2.0;
				double s_u = sqrt(pow2(n) / 12.0 * ((2*n + 1) - tie_corr / (2 * n * (2 * n - 1))));

				double cont_term = 0.5;

				double z = (std::max(ret.statistic_U1, ret.statistic_U2) - m_u - cont_term) / s_u;
				ret.p_value = std::clamp<double>(2 * (1.0 - stats::pnorm(z, 0, 1)), 0, 1);
			}
			else
			{
				// !!! TODO
			}

			return ret;
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static snr_t SNR_test(Iter1 X_first, Iter1 X_last, Iter2 Y_first)
		{
			return SNR_test_n(X_first, Y_first, std::distance(X_first, X_last));
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static snr_t SNR_test_n(Iter1 X_first, Iter2 Y_first, size_t n)
		{
			snr_t ret;

			double m1 = mean(X_first, n);
			double m2 = mean(Y_first, n);
			double sd1 = std_dev(m1, X_first, n);
			double sd2 = std_dev(m2, Y_first, n);

			ret.ratio = (m1 - m2) / (sd1 + sd2);

			ret.p_value = -1;			// !!! TODO

			return ret;
		}
	};

	// *************************************************************************************
	// 
	// *************************************************************************************
	class p_val_correction
	{
		std::vector<std::pair<double, size_t>> pv_orig;
		std::vector<double> H_n;

		template<typename Iter>
		void load_and_sort(Iter first, size_t n)
		{
			pv_orig.resize(n);

			for (size_t i = 0; i < n; ++i)
			{
				pv_orig[i].first = (double)*first++;
				pv_orig[i].second = i;
			}

			std::sort(pv_orig.begin(), pv_orig.end());
		}

		void build_Hn(size_t n)
		{
			// H_0 = 0
			// H_1 = 1
			// H_2 = 1 + 1/2

			if (n < H_n.size())
				return;

			size_t n_ready = std::max<size_t>(1, H_n.size());
			H_n.resize(n + 1);

			for (size_t i = n_ready; i <= n; ++i)
				H_n[i] = H_n[i - 1] + 1.0 / i;
		}

		template<typename Iter1, typename Iter2>
		void benjamini_base_n(Iter1 S_first, Iter2 D_first, size_t n, double cm)
		{
			load_and_sort(S_first, n);

			double min_pv = 1.0;
			double mult = n * cm;

			for (size_t i = n; i > 0; --i)
			{
				double pv = pv_orig[i - 1].first * mult / i;
				min_pv = std::min(min_pv, pv);
				*(D_first + pv_orig[i - 1].second) = (decltype(*D_first))min_pv;
			}
		}

	public:
		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static void bonferroni(Iter1 S_first, Iter1 S_last, Iter2 D_first)
		{
			bonferroni_n(S_first, D_first, std::distance(S_first, S_last));
		}
			
		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static void bonferroni_n(Iter1 S_first, Iter2 D_first, size_t n)
		{
			for (size_t i = 0; i < n; ++i)
				*D_first++ = (decltype(*D_first))std::clamp<double>(*S_first++ * (double)n, 0, 1);
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void holm_bonferroni(Iter1 S_first, Iter1 S_last, Iter2 D_first)
		{
			holm_bonferroni_n(S_first, D_first, std::distance(S_first, S_last));
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void holm_bonferroni_n(Iter1 S_first, Iter2 D_first, size_t n)
		{
			load_and_sort(S_first, n);

			double max_pv = 0.0;
			for (size_t i = 0; i < n; ++i)
			{
				double pv = pv_orig[i].first * (double)(n - i);
				max_pv = std::min(std::max(max_pv, pv), 1.0);
				*(D_first + pv_orig[i].second) = (decltype(*D_first))max_pv;
			}
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void benjamini_hochberg(Iter1 S_first, Iter1 S_last, Iter2 D_first)
		{
			benjamini_hochberg_n(S_first, D_first, std::distance(S_first, S_last));
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void benjamini_hochberg_n(Iter1 S_first, Iter2 D_first, size_t n)
		{
			benjamini_base_n(S_first, D_first, n, 1.0);
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void benjamini_yekutieli(Iter1 S_first, Iter1 S_last, Iter2 D_first)
		{
			benjamini_yekutieli_n(S_first, D_first, std::distance(S_first, S_last));
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void benjamini_yekutieli_n(Iter1 S_first, Iter2 D_first, size_t n)
		{
			build_Hn(n);

			benjamini_base_n(S_first, D_first, n, H_n[n]);
		}
	};
}

#endif
