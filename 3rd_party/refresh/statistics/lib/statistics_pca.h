#ifndef _STATISTICS_PCA_H
#define _STATISTICS_PCA_H

#define EIGEN_NO_DEPRECATED_WARNING
#include <Eigen/Dense>
#include <vector>
#include <iterator>
#include <cinttypes>
#include <numeric>
#include <algorithm>
#include <cassert>

namespace refresh
{
	template<typename VALUE_T>
	class pca
	{
	public:
		enum class computation_mode_t { covariance, svd };

	protected:
		enum class input_mode_t { unknown, feature_oriented, object_oriented };

		using entry_t = std::vector<VALUE_T>;

		std::vector<entry_t> input_data;
		std::vector<entry_t> ret_data;

		std::vector<VALUE_T> explained_variance;
		std::vector<VALUE_T> explained_variance_ratio;

		input_mode_t input_mode = input_mode_t::unknown;
		size_t input_vector_size{};

		template<typename VEC>
		void vec_clear(VEC& v)
		{
			v.clear();
			v.shrink_to_fit();
		}

		template<typename Iter>
		bool try_add(input_mode_t requested_mode, Iter first, Iter last)
		{
			if (input_mode == input_mode_t::unknown)
			{
				input_mode = requested_mode;
				input_vector_size = std::distance(first, last);
			}
			else if (input_mode != requested_mode)
				return false;

			if (static_cast<size_t>(std::distance(first, last)) != input_vector_size)
				return false;

			input_data.emplace_back(first, last);
			return true;
		}

		void clear()
		{
			vec_clear(input_data);
			vec_clear(ret_data);
			vec_clear(explained_variance);
			vec_clear(explained_variance_ratio);
			input_mode = input_mode_t::unknown;
			input_vector_size = 0;
		}

		// Pomocnicza funkcja do budowania macierzy Eigen z input_data
		Eigen::MatrixXd build_matrix() const
		{
			size_t no_objects = 0;
			size_t no_features = 0;
			Eigen::MatrixXd X;

			if (input_mode == input_mode_t::object_oriented)
			{
				no_objects = input_data.size();
				no_features = input_data.front().size();
				X.resize(no_objects, no_features);
				for (size_t i = 0; i < no_objects; ++i)
					for (size_t j = 0; j < no_features; ++j)
						X(i, j) = static_cast<double>(input_data[i][j]);
			}
			else if (input_mode == input_mode_t::feature_oriented)
			{
				no_objects = input_data.front().size();
				no_features = input_data.size();
				X.resize(no_objects, no_features);
				for (size_t i = 0; i < no_objects; ++i)
					for (size_t j = 0; j < no_features; ++j)
						X(i, j) = static_cast<double>(input_data[j][i]);
			}
			return X;
		}

		bool run_covariance(size_t no_dimensions = 2)
		{
			if (input_mode == input_mode_t::unknown || input_data.empty())
				return false;

			Eigen::MatrixXd X = build_matrix();
			size_t no_objects = X.rows();

			auto centered = X.rowwise() - X.colwise().mean();

			auto covariance = (centered.adjoint() * centered) / double(no_objects - 1);

			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(covariance);
			if (eig.info() != Eigen::Success)
				return false;

			Eigen::VectorXd eigenvalues = eig.eigenvalues().reverse();
			Eigen::MatrixXd eigenvectors = eig.eigenvectors().colwise().reverse();

			explained_variance.resize(eigenvalues.size());
			for (int i = 0; i < eigenvalues.size(); ++i)
				explained_variance[i] = static_cast<VALUE_T>(eigenvalues(i));

			double total = std::accumulate(explained_variance.begin(), explained_variance.end(), 0.0);
			explained_variance_ratio.resize(explained_variance.size());
			for (size_t i = 0; i < explained_variance.size(); ++i)
				explained_variance_ratio[i] = static_cast<VALUE_T>(double(explained_variance[i]) / total);

			auto principal_components = eigenvectors.leftCols(std::min(no_dimensions, (size_t)eigenvectors.cols()));
			auto transformed = centered * principal_components;

			ret_data.assign(no_objects, entry_t(transformed.cols()));
			for (size_t i = 0; i < no_objects; ++i)
				for (size_t j = 0; j < (size_t)transformed.cols(); ++j)
					ret_data[i][j] = static_cast<VALUE_T>(transformed(i, j));

			return true;
		}

		bool run_svd(size_t no_dimensions = 2)
		{
			if (input_mode == input_mode_t::unknown || input_data.empty())
				return false;

			Eigen::MatrixXd X = build_matrix();
			size_t no_objects = X.rows();

			auto centered = X.rowwise() - X.colwise().mean();

			Eigen::JacobiSVD<Eigen::MatrixXd> svd(centered, Eigen::ComputeThinU | Eigen::ComputeThinV);

			Eigen::MatrixXd principal_components = svd.matrixV().leftCols(std::min(no_dimensions, (size_t)svd.matrixV().cols()));
			auto transformed = centered * principal_components;

			Eigen::VectorXd s = svd.singularValues();
			explained_variance.resize(s.size());
			for (int i = 0; i < s.size(); ++i)
				explained_variance[i] = static_cast<VALUE_T>((s(i) * s(i)) / (double(no_objects) - 1.0));

			double total = std::accumulate(explained_variance.begin(), explained_variance.end(), 0.0);
			explained_variance_ratio.resize(explained_variance.size());
			for (size_t i = 0; i < explained_variance.size(); ++i)
				explained_variance_ratio[i] = static_cast<VALUE_T>(double(explained_variance[i]) / total);

			ret_data.assign(no_objects, entry_t(transformed.cols()));
			for (size_t i = 0; i < no_objects; ++i)
				for (size_t j = 0; j < (size_t)transformed.cols(); ++j)
					ret_data[i][j] = static_cast<VALUE_T>(transformed(i, j));

			return true;
		}

	public:
		pca() = default;

		void reset() { clear(); }

		template<typename Iter>
		bool add_object(Iter first, Iter last) { return try_add(input_mode_t::object_oriented, first, last); }

		template<typename Iter>
		bool add_feature(Iter first, Iter last) { return try_add(input_mode_t::feature_oriented, first, last); }

		bool run(size_t no_dimensions = 2, computation_mode_t computation_mode = computation_mode_t::svd)
		{
			if (computation_mode == computation_mode_t::covariance)
				return run_covariance(no_dimensions);
			else if (computation_mode == computation_mode_t::svd)
				return run_svd(no_dimensions);
			return false;
		}

		const std::vector<std::vector<VALUE_T>>& result() const noexcept { return ret_data; }
		const std::vector<VALUE_T>& get_explained_variance() const noexcept { return explained_variance; }
		const std::vector<VALUE_T>& get_explained_variance_ratio() const noexcept { return explained_variance_ratio; }
	};

	template<typename VALUE_T>
	class pca_parallel_add : public pca<VALUE_T>
	{
		size_t no_features_val;
	public:
		pca_parallel_add(size_t no_objects, size_t no_features) : no_features_val(no_features)
		{
			this->input_vector_size = no_objects;
			this->input_mode = pca<VALUE_T>::input_mode_t::feature_oriented;
			this->input_data.assign(no_features, typename pca<VALUE_T>::entry_t(no_objects));
		}

		bool add(size_t object_idx, size_t feature_idx, const VALUE_T value)
		{
			if (feature_idx >= this->input_data.size() || object_idx >= this->input_vector_size) return false;
			this->input_data[feature_idx][object_idx] = value;
			return true;
		}
	};
}

#endif