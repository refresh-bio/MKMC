#ifndef _STATISTICS_PCA_H
#define _STATISTICS_PCA_H

#define EIGEN_NO_DEPRECATED_WARNING
#include "umappp/Umap.hpp"

#include <vector>
#include <iterator>
#include <cinttypes>

namespace refresh
{
	template<typename VALUE_T>
	class pca
	{
	public:
		enum class computation_mode_t { covariance, svd };

	private:
		enum class input_mode_t {unknown, feature_oriented, object_oriented};

		using entry_t = std::vector<VALUE_T>;

		std::vector<entry_t> input_data;
		std::vector<entry_t> ret_data;

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

			if (std::distance(first, last) != input_vector_size)
				return false;

			input_data.emplace_back(first, last);

			return true;
		}

		void clear()
		{
			vec_clear(input_data);
			vec_clear(ret_data);

			input_mode = input_mode_t::unknown;
			input_vector_size = 0;
		}

		bool run_covariance(size_t no_dimensions = 2)
		{
			if (input_mode == input_mode_t::unknown)
				return false;

			size_t k = 0;
			size_t no_objects = 0;
			size_t no_features = 0;

			Eigen::MatrixXd X;

			switch (input_mode)
			{
			case input_mode_t::object_oriented:
				no_objects = input_data.size();
				no_features = input_data.front().size();

				X.resize(no_objects, no_features);

				for (size_t i = 0; i < no_objects; ++i)
					for (size_t j = 0; j < no_features; ++j)
						X(i, j) = input_data[i][j];

				break;
			case input_mode_t::feature_oriented:
				no_objects = input_data.front().size();
				no_features = input_data.size();

				X.resize(no_objects, no_features);

				for (size_t i = 0; i < no_objects; ++i)
					for (size_t j = 0; j < no_features; ++j)
						X(i, j) = input_data[j][i];

				break;
			default:
				assert(0);						// Never should be here
			}

			// PCA calculations
			auto centered = X.rowwise() - X.colwise().mean();

			auto covariance = (centered.adjoint() * centered) / double(X.rows() - 1);

			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(covariance);

			if (eig.info() != Eigen::Success)
				return false;

			Eigen::VectorXd eigenvalues = eig.eigenvalues().reverse();
			Eigen::MatrixXd eigenvectors = eig.eigenvectors().rowwise().reverse();

			auto principal_components = eigenvectors.leftCols(no_dimensions);

			auto transformed = centered * principal_components;

			// Convert results to vector of vectors
			ret_data.clear();
			ret_data.resize(no_objects, entry_t(no_dimensions));

			k = 0;
			for (size_t i = 0; i < no_objects; ++i)
			{
				ret_data[i].resize(no_dimensions);
				for (size_t j = 0; j < no_dimensions; ++j)
					ret_data[i][j] = transformed(i, j);
			}

			return true;
		}

		bool run_svd(size_t no_dimensions = 2)
		{
			if (input_mode == input_mode_t::unknown)
				return false;

			size_t k = 0;
			size_t no_objects = 0;
			size_t no_features = 0;

			Eigen::MatrixXd X;

			switch (input_mode)
			{
			case input_mode_t::object_oriented:
				no_objects = input_data.size();
				no_features = input_data.front().size();

				X.resize(no_objects, no_features);

				for (size_t i = 0; i < no_objects; ++i)
					for (size_t j = 0; j < no_features; ++j)
						X(i, j) = input_data[i][j];

				break;
			case input_mode_t::feature_oriented:
				no_objects = input_data.front().size();
				no_features = input_data.size();

				X.resize(no_objects, no_features);

				for (size_t i = 0; i < no_objects; ++i)
					for (size_t j = 0; j < no_features; ++j)
						X(i, j) = input_data[j][i];

				break;
			default:
				assert(0);						// Never should be here
			}

			// PCA calculations
			auto centered = X.rowwise() - X.colwise().mean();

			Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);

			Eigen::MatrixXd principal_components = svd.matrixV().leftCols(no_dimensions);
			
			auto transformed = centered * principal_components;			

			// Convert results to vector of vectors
			ret_data.clear();
			ret_data.resize(no_objects, entry_t(no_dimensions));

			k = 0;
			for (size_t i = 0; i < no_objects; ++i)
			{
				ret_data[i].resize(no_dimensions);
				for (size_t j = 0; j < no_dimensions; ++j)
					ret_data[i][j] = transformed(i, j);
			}

			return true;
		}

	public:
		pca() = default;

		void reset()
		{
			clear();
		}

		template<typename Iter>
		bool add_object(Iter first, Iter last)
		{
			return try_add(input_mode_t::object_oriented, first, last);
		}

		template<typename Iter>
		bool add_feature(Iter first, Iter last)
		{
			return try_add(input_mode_t::feature_oriented, first, last);
		}

		bool run(size_t no_dimensions = 2, computation_mode_t computation_mode = computation_mode_t::svd)
		{
			if (computation_mode == computation_mode_t::covariance)
				return run_covariance(no_dimensions);
			else if (computation_mode == computation_mode_t::svd)
				return run_svd(no_dimensions);
			return false;
		}


		std::vector<std::vector<VALUE_T>>& result()
		{
			return ret_data;
		}
	};
}

#endif