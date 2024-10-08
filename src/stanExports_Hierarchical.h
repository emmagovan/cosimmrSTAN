// Generated by rstantools.  Do not edit by hand.

/*
    cosimmrSTAN is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cosimmrSTAN is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cosimmrSTAN.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.32.2
#include <stan/model/model_header.hpp>
namespace model_Hierarchical_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 35> locations_array__ =
  {" (found before start of program)",
  " (in 'Hierarchical', line 9, column 3 to column 31)",
  " (in 'Hierarchical', line 10, column 2 to column 17)",
  " (in 'Hierarchical', line 13, column 3 to column 24)",
  " (in 'Hierarchical', line 30, column 2 to column 20)",
  " (in 'Hierarchical', line 17, column 6 to column 67)",
  " (in 'Hierarchical', line 16, column 19 to line 18, column 5)",
  " (in 'Hierarchical', line 16, column 4 to line 18, column 5)",
  " (in 'Hierarchical', line 15, column 17 to line 19, column 3)",
  " (in 'Hierarchical', line 15, column 2 to line 19, column 3)",
  " (in 'Hierarchical', line 31, column 2 to column 26)",
  " (in 'Hierarchical', line 24, column 0 to column 57)",
  " (in 'Hierarchical', line 25, column 0 to column 30)",
  " (in 'Hierarchical', line 23, column 13 to line 26, column 1)",
  " (in 'Hierarchical', line 23, column 0 to line 26, column 1)",
  " (in 'Hierarchical', line 22, column 13 to line 27, column 1)",
  " (in 'Hierarchical', line 22, column 0 to line 27, column 1)",
  " (in 'Hierarchical', line 2, column 4 to column 19)",
  " (in 'Hierarchical', line 3, column 4 to column 19)",
  " (in 'Hierarchical', line 4, column 17 to column 18)",
  " (in 'Hierarchical', line 4, column 2 to column 20)",
  " (in 'Hierarchical', line 5, column 9 to column 10)",
  " (in 'Hierarchical', line 5, column 12 to column 13)",
  " (in 'Hierarchical', line 5, column 2 to column 27)",
  " (in 'Hierarchical', line 6, column 9 to column 10)",
  " (in 'Hierarchical', line 6, column 11 to column 12)",
  " (in 'Hierarchical', line 6, column 2 to column 24)",
  " (in 'Hierarchical', line 9, column 19 to column 20)",
  " (in 'Hierarchical', line 9, column 22 to column 23)",
  " (in 'Hierarchical', line 10, column 9 to column 10)",
  " (in 'Hierarchical', line 10, column 11 to column 12)",
  " (in 'Hierarchical', line 13, column 10 to column 11)",
  " (in 'Hierarchical', line 13, column 13 to column 14)",
  " (in 'Hierarchical', line 30, column 9 to column 10)",
  " (in 'Hierarchical', line 30, column 11 to column 12)"};
#include <stan_meta_header.hpp>
class model_Hierarchical final : public model_base_crtp<model_Hierarchical> {
private:
  int J;
  int K;
  std::vector<int> n;
  Eigen::Matrix<double,-1,-1> source_mean_data__;
  Eigen::Matrix<double,-1,-1> source_sd_data__;
  Eigen::Map<Eigen::Matrix<double,-1,-1>> source_mean{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> source_sd{nullptr, 0, 0};
public:
  ~model_Hierarchical() {}
  model_Hierarchical(stan::io::var_context& context__, unsigned int
                     random_seed__ = 0, std::ostream* pstream__ = nullptr)
      : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_Hierarchical_namespace::model_Hierarchical";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 17;
      context__.validate_dims("data initialization", "J", "int",
        std::vector<size_t>{});
      J = std::numeric_limits<int>::min();
      current_statement__ = 17;
      J = context__.vals_i("J")[(1 - 1)];
      current_statement__ = 17;
      stan::math::check_greater_or_equal(function__, "J", J, 1);
      current_statement__ = 18;
      context__.validate_dims("data initialization", "K", "int",
        std::vector<size_t>{});
      K = std::numeric_limits<int>::min();
      current_statement__ = 18;
      K = context__.vals_i("K")[(1 - 1)];
      current_statement__ = 18;
      stan::math::check_greater_or_equal(function__, "K", K, 1);
      current_statement__ = 19;
      stan::math::validate_non_negative_index("n", "K", K);
      current_statement__ = 20;
      context__.validate_dims("data initialization", "n", "int",
        std::vector<size_t>{static_cast<size_t>(K)});
      n = std::vector<int>(K, std::numeric_limits<int>::min());
      current_statement__ = 20;
      n = context__.vals_i("n");
      current_statement__ = 20;
      stan::math::check_greater_or_equal(function__, "n", n, 1);
      current_statement__ = 21;
      stan::math::validate_non_negative_index("source_mean", "K", K);
      current_statement__ = 22;
      stan::math::validate_non_negative_index("source_mean", "J", J);
      current_statement__ = 23;
      context__.validate_dims("data initialization", "source_mean", "double",
        std::vector<size_t>{static_cast<size_t>(K), static_cast<size_t>(J)});
      source_mean_data__ = Eigen::Matrix<double,-1,-1>::Constant(K, J,
                             std::numeric_limits<double>::quiet_NaN());
      new (&source_mean)
        Eigen::Map<Eigen::Matrix<double,-1,-1>>(source_mean_data__.data(), K,
        J);
      {
        std::vector<local_scalar_t__> source_mean_flat__;
        current_statement__ = 23;
        source_mean_flat__ = context__.vals_r("source_mean");
        current_statement__ = 23;
        pos__ = 1;
        current_statement__ = 23;
        for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
          current_statement__ = 23;
          for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
            current_statement__ = 23;
            stan::model::assign(source_mean, source_mean_flat__[(pos__ - 1)],
              "assigning variable source_mean",
              stan::model::index_uni(sym2__), stan::model::index_uni(sym1__));
            current_statement__ = 23;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 24;
      stan::math::validate_non_negative_index("source_sd", "K", K);
      current_statement__ = 25;
      stan::math::validate_non_negative_index("source_sd", "J", J);
      current_statement__ = 26;
      context__.validate_dims("data initialization", "source_sd", "double",
        std::vector<size_t>{static_cast<size_t>(K), static_cast<size_t>(J)});
      source_sd_data__ = Eigen::Matrix<double,-1,-1>::Constant(K, J,
                           std::numeric_limits<double>::quiet_NaN());
      new (&source_sd)
        Eigen::Map<Eigen::Matrix<double,-1,-1>>(source_sd_data__.data(), K,
        J);
      {
        std::vector<local_scalar_t__> source_sd_flat__;
        current_statement__ = 26;
        source_sd_flat__ = context__.vals_r("source_sd");
        current_statement__ = 26;
        pos__ = 1;
        current_statement__ = 26;
        for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
          current_statement__ = 26;
          for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
            current_statement__ = 26;
            stan::model::assign(source_sd, source_sd_flat__[(pos__ - 1)],
              "assigning variable source_sd", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 26;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 27;
      stan::math::validate_non_negative_index("tmp_X", "K", K);
      current_statement__ = 28;
      stan::math::validate_non_negative_index("tmp_X", "J", J);
      current_statement__ = 29;
      stan::math::validate_non_negative_index("mu", "K", K);
      current_statement__ = 30;
      stan::math::validate_non_negative_index("mu", "J", J);
      current_statement__ = 31;
      stan::math::validate_non_negative_index("src_tau", "K", K);
      current_statement__ = 32;
      stan::math::validate_non_negative_index("src_tau", "J", J);
      current_statement__ = 33;
      stan::math::validate_non_negative_index("sigma", "K", K);
      current_statement__ = 34;
      stan::math::validate_non_negative_index("sigma", "J", J);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = (K * J) + (K * J);
  }
  inline std::string model_name() const final {
    return "model_Hierarchical";
  }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.32.2",
             "stancflags = --allow-undefined"};
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI,
            stan::require_vector_like_t<VecR>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR>
  log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
                pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    static constexpr const char* function__ =
      "model_Hierarchical_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<local_scalar_t__,-1,-1> tmp_X =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, J, DUMMY_VAR__);
      current_statement__ = 1;
      tmp_X = in__.template read_constrain_lb<
                Eigen::Matrix<local_scalar_t__,-1,-1>, jacobian__>(0, lp__,
                K, J);
      Eigen::Matrix<local_scalar_t__,-1,-1> mu =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, J, DUMMY_VAR__);
      current_statement__ = 2;
      mu = in__.template read<Eigen::Matrix<local_scalar_t__,-1,-1>>(K, J);
      Eigen::Matrix<local_scalar_t__,-1,-1> src_tau =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, J, DUMMY_VAR__);
      current_statement__ = 9;
      for (int k = 1; k <= K; ++k) {
        current_statement__ = 7;
        for (int j = 1; j <= J; ++j) {
          current_statement__ = 5;
          stan::model::assign(src_tau,
            (stan::model::rvalue(tmp_X, "tmp_X", stan::model::index_uni(k),
               stan::model::index_uni(j)) /
            (stan::model::rvalue(source_sd, "source_sd",
               stan::model::index_uni(k), stan::model::index_uni(j)) *
            (stan::model::rvalue(n, "n", stan::model::index_uni(k)) - 1))),
            "assigning variable src_tau", stan::model::index_uni(k),
            stan::model::index_uni(j));
        }
      }
      {
        current_statement__ = 16;
        for (int k = 1; k <= K; ++k) {
          current_statement__ = 14;
          for (int j = 1; j <= J; ++j) {
            current_statement__ = 11;
            lp_accum__.add(stan::math::normal_lpdf<propto__>(
                             stan::model::rvalue(mu, "mu",
                               stan::model::index_uni(k),
                               stan::model::index_uni(j)),
                             stan::model::rvalue(source_mean, "source_mean",
                               stan::model::index_uni(k),
                               stan::model::index_uni(j)),
                             (stan::model::rvalue(n, "n",
                                stan::model::index_uni(k)) /
                             stan::model::rvalue(source_sd, "source_sd",
                               stan::model::index_uni(k),
                               stan::model::index_uni(j)))));
            current_statement__ = 12;
            lp_accum__.add(stan::math::chi_square_lpdf<propto__>(
                             stan::model::rvalue(tmp_X, "tmp_X",
                               stan::model::index_uni(k),
                               stan::model::index_uni(j)),
                             stan::model::rvalue(n, "n",
                               stan::model::index_uni(k))));
          }
        }
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
  }
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
            stan::require_vector_like_vt<std::is_floating_point,
            VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
            VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
            VecVar>* = nullptr>
  inline void
  write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
                   VecVar& vars__, const bool
                   emit_transformed_parameters__ = true, const bool
                   emit_generated_quantities__ = true, std::ostream*
                   pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    // suppress unused var warning
    (void) propto__;
    double lp__ = 0.0;
    // suppress unused var warning
    (void) lp__;
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    constexpr bool jacobian__ = false;
    static constexpr const char* function__ =
      "model_Hierarchical_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<double,-1,-1> tmp_X =
        Eigen::Matrix<double,-1,-1>::Constant(K, J,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 1;
      tmp_X = in__.template read_constrain_lb<
                Eigen::Matrix<local_scalar_t__,-1,-1>, jacobian__>(0, lp__,
                K, J);
      Eigen::Matrix<double,-1,-1> mu =
        Eigen::Matrix<double,-1,-1>::Constant(K, J,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 2;
      mu = in__.template read<Eigen::Matrix<local_scalar_t__,-1,-1>>(K, J);
      Eigen::Matrix<double,-1,-1> src_tau =
        Eigen::Matrix<double,-1,-1>::Constant(K, J,
          std::numeric_limits<double>::quiet_NaN());
      out__.write(tmp_X);
      out__.write(mu);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      current_statement__ = 9;
      for (int k = 1; k <= K; ++k) {
        current_statement__ = 7;
        for (int j = 1; j <= J; ++j) {
          current_statement__ = 5;
          stan::model::assign(src_tau,
            (stan::model::rvalue(tmp_X, "tmp_X", stan::model::index_uni(k),
               stan::model::index_uni(j)) /
            (stan::model::rvalue(source_sd, "source_sd",
               stan::model::index_uni(k), stan::model::index_uni(j)) *
            (stan::model::rvalue(n, "n", stan::model::index_uni(k)) - 1))),
            "assigning variable src_tau", stan::model::index_uni(k),
            stan::model::index_uni(j));
        }
      }
      if (emit_transformed_parameters__) {
        out__.write(src_tau);
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
      Eigen::Matrix<double,-1,-1> sigma =
        Eigen::Matrix<double,-1,-1>::Constant(K, J,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 10;
      stan::model::assign(sigma,
        stan::math::divide(1, stan::math::sqrt(src_tau)),
        "assigning variable sigma");
      out__.write(sigma);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, typename VecI,
            stan::require_vector_t<VecVar>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void
  unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
                         VecVar& vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,-1> tmp_X =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, J, DUMMY_VAR__);
      current_statement__ = 1;
      stan::model::assign(tmp_X,
        in__.read<Eigen::Matrix<local_scalar_t__,-1,-1>>(K, J),
        "assigning variable tmp_X");
      out__.write_free_lb(0, tmp_X);
      Eigen::Matrix<local_scalar_t__,-1,-1> mu =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, J, DUMMY_VAR__);
      current_statement__ = 2;
      stan::model::assign(mu,
        in__.read<Eigen::Matrix<local_scalar_t__,-1,-1>>(K, J),
        "assigning variable mu");
      out__.write(mu);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
  inline void
  transform_inits_impl(const stan::io::var_context& context__, VecVar&
                       vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      current_statement__ = 1;
      context__.validate_dims("parameter initialization", "tmp_X", "double",
        std::vector<size_t>{static_cast<size_t>(K), static_cast<size_t>(J)});
      current_statement__ = 2;
      context__.validate_dims("parameter initialization", "mu", "double",
        std::vector<size_t>{static_cast<size_t>(K), static_cast<size_t>(J)});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,-1> tmp_X =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, J, DUMMY_VAR__);
      {
        std::vector<local_scalar_t__> tmp_X_flat__;
        current_statement__ = 1;
        tmp_X_flat__ = context__.vals_r("tmp_X");
        current_statement__ = 1;
        pos__ = 1;
        current_statement__ = 1;
        for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
          current_statement__ = 1;
          for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
            current_statement__ = 1;
            stan::model::assign(tmp_X, tmp_X_flat__[(pos__ - 1)],
              "assigning variable tmp_X", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 1;
            pos__ = (pos__ + 1);
          }
        }
      }
      out__.write_free_lb(0, tmp_X);
      Eigen::Matrix<local_scalar_t__,-1,-1> mu =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, J, DUMMY_VAR__);
      {
        std::vector<local_scalar_t__> mu_flat__;
        current_statement__ = 2;
        mu_flat__ = context__.vals_r("mu");
        current_statement__ = 2;
        pos__ = 1;
        current_statement__ = 2;
        for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
          current_statement__ = 2;
          for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
            current_statement__ = 2;
            stan::model::assign(mu, mu_flat__[(pos__ - 1)],
              "assigning variable mu", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 2;
            pos__ = (pos__ + 1);
          }
        }
      }
      out__.write(mu);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"tmp_X", "mu"};
    if (emit_transformed_parameters__) {
      std::vector<std::string> temp{"src_tau"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {
      std::vector<std::string> temp{"sigma"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{static_cast<
                                                                    size_t>(K),
                                                 static_cast<size_t>(J)},
                std::vector<size_t>{static_cast<size_t>(K),
                  static_cast<size_t>(J)}};
    if (emit_transformed_parameters__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(K),
               static_cast<size_t>(J)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(K),
               static_cast<size_t>(J)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
        param_names__.emplace_back(std::string() + "tmp_X" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
        param_names__.emplace_back(std::string() + "mu" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
          param_names__.emplace_back(std::string() + "src_tau" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
    }
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
          param_names__.emplace_back(std::string() + "sigma" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
    }
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
        param_names__.emplace_back(std::string() + "tmp_X" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
        param_names__.emplace_back(std::string() + "mu" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
          param_names__.emplace_back(std::string() + "src_tau" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
    }
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
          param_names__.emplace_back(std::string() + "sigma" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
    }
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"tmp_X\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(J) + "},\"block\":\"parameters\"},{\"name\":\"mu\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(J) + "},\"block\":\"parameters\"},{\"name\":\"src_tau\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(J) + "},\"block\":\"transformed_parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(J) + "},\"block\":\"generated_quantities\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"tmp_X\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(J) + "},\"block\":\"parameters\"},{\"name\":\"mu\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(J) + "},\"block\":\"parameters\"},{\"name\":\"src_tau\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(J) + "},\"block\":\"transformed_parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(J) + "},\"block\":\"generated_quantities\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = ((K * J) + (K * J));
    const size_t num_transformed = emit_transformed_parameters * ((K * J));
    const size_t num_gen_quantities = emit_generated_quantities * ((K * J));
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    std::vector<int> params_i;
    vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <typename RNG> inline void
  write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
              params_i, std::vector<double>& vars, bool
              emit_transformed_parameters = true, bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = ((K * J) + (K * J));
    const size_t num_transformed = emit_transformed_parameters * ((K * J));
    const size_t num_gen_quantities = emit_generated_quantities * ((K * J));
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    vars = std::vector<double>(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
    Eigen::Matrix<int,-1,1> params_i;
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
           std::ostream* pstream = nullptr) const {
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  inline void
  transform_inits(const stan::io::var_context& context,
                  Eigen::Matrix<double,-1,1>& params_r, std::ostream*
                  pstream = nullptr) const final {
    std::vector<double> params_r_vec(params_r.size());
    std::vector<int> params_i;
    transform_inits(context, params_i, params_r_vec, pstream);
    params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
                 params_r_vec.size());
  }
  inline void
  transform_inits(const stan::io::var_context& context, std::vector<int>&
                  params_i, std::vector<double>& vars, std::ostream*
                  pstream__ = nullptr) const {
    vars.resize(num_params_r__);
    transform_inits_impl(context, vars, pstream__);
  }
  inline void
  unconstrain_array(const std::vector<double>& params_constrained,
                    std::vector<double>& params_unconstrained, std::ostream*
                    pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = std::vector<double>(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
  inline void
  unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
                    Eigen::Matrix<double,-1,1>& params_unconstrained,
                    std::ostream* pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
};
}
using stan_model = model_Hierarchical_namespace::model_Hierarchical;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_Hierarchical_namespace::profiles__;
}
#endif
#endif
