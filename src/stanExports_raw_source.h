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
namespace model_raw_source_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 37> locations_array__ =
  {" (found before start of program)",
  " (in 'raw_source', line 10, column 1 to column 19)",
  " (in 'raw_source', line 11, column 1 to column 25)",
  " (in 'raw_source', line 12, column 1 to column 24)",
  " (in 'raw_source', line 15, column 4 to column 25)",
  " (in 'raw_source', line 17, column 4 to column 65)",
  " (in 'raw_source', line 16, column 20 to line 18, column 3)",
  " (in 'raw_source', line 16, column 5 to line 18, column 3)",
  " (in 'raw_source', line 21, column 3 to column 24)",
  " (in 'raw_source', line 25, column 6 to column 33)",
  " (in 'raw_source', line 24, column 17 to line 26, column 5)",
  " (in 'raw_source', line 24, column 4 to line 26, column 5)",
  " (in 'raw_source', line 23, column 15 to line 27, column 3)",
  " (in 'raw_source', line 23, column 2 to line 27, column 3)",
  " (in 'raw_source', line 29, column 2 to column 33)",
  " (in 'raw_source', line 28, column 13 to line 30, column 1)",
  " (in 'raw_source', line 28, column 0 to line 30, column 1)",
  " (in 'raw_source', line 33, column 6 to column 63)",
  " (in 'raw_source', line 32, column 17 to line 34, column 5)",
  " (in 'raw_source', line 32, column 2 to line 34, column 5)",
  " (in 'raw_source', line 2, column 2 to column 17)",
  " (in 'raw_source', line 3, column 2 to column 17)",
  " (in 'raw_source', line 4, column 2 to column 17)",
  " (in 'raw_source', line 5, column 9 to column 10)",
  " (in 'raw_source', line 5, column 12 to column 13)",
  " (in 'raw_source', line 5, column 2 to column 17)",
  " (in 'raw_source', line 6, column 31 to column 32)",
  " (in 'raw_source', line 6, column 2 to column 34)",
  " (in 'raw_source', line 7, column 2 to column 26)",
  " (in 'raw_source', line 10, column 8 to column 9)",
  " (in 'raw_source', line 10, column 10 to column 11)",
  " (in 'raw_source', line 11, column 22 to column 23)",
  " (in 'raw_source', line 11, column 13 to column 14)",
  " (in 'raw_source', line 12, column 17 to column 18)",
  " (in 'raw_source', line 15, column 22 to column 23)",
  " (in 'raw_source', line 15, column 11 to column 12)",
  " (in 'raw_source', line 15, column 13 to column 14)"};
#include <stan_meta_header.hpp>
class model_raw_source final : public model_base_crtp<model_raw_source> {
private:
  int J;
  int K;
  int N;
  Eigen::Matrix<double,-1,-1> Y_data__;
  std::vector<int> source;
  double shape_sig;
  Eigen::Map<Eigen::Matrix<double,-1,-1>> Y{nullptr, 0, 0};
public:
  ~model_raw_source() {}
  model_raw_source(stan::io::var_context& context__, unsigned int
                   random_seed__ = 0, std::ostream* pstream__ = nullptr)
      : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_raw_source_namespace::model_raw_source";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 20;
      context__.validate_dims("data initialization", "J", "int",
        std::vector<size_t>{});
      J = std::numeric_limits<int>::min();
      current_statement__ = 20;
      J = context__.vals_i("J")[(1 - 1)];
      current_statement__ = 20;
      stan::math::check_greater_or_equal(function__, "J", J, 1);
      current_statement__ = 21;
      context__.validate_dims("data initialization", "K", "int",
        std::vector<size_t>{});
      K = std::numeric_limits<int>::min();
      current_statement__ = 21;
      K = context__.vals_i("K")[(1 - 1)];
      current_statement__ = 21;
      stan::math::check_greater_or_equal(function__, "K", K, 1);
      current_statement__ = 22;
      context__.validate_dims("data initialization", "N", "int",
        std::vector<size_t>{});
      N = std::numeric_limits<int>::min();
      current_statement__ = 22;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 22;
      stan::math::check_greater_or_equal(function__, "N", N, 1);
      current_statement__ = 23;
      stan::math::validate_non_negative_index("Y", "N", N);
      current_statement__ = 24;
      stan::math::validate_non_negative_index("Y", "J", J);
      current_statement__ = 25;
      context__.validate_dims("data initialization", "Y", "double",
        std::vector<size_t>{static_cast<size_t>(N), static_cast<size_t>(J)});
      Y_data__ = Eigen::Matrix<double,-1,-1>::Constant(N, J,
                   std::numeric_limits<double>::quiet_NaN());
      new (&Y) Eigen::Map<Eigen::Matrix<double,-1,-1>>(Y_data__.data(), N, J);
      {
        std::vector<local_scalar_t__> Y_flat__;
        current_statement__ = 25;
        Y_flat__ = context__.vals_r("Y");
        current_statement__ = 25;
        pos__ = 1;
        current_statement__ = 25;
        for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
          current_statement__ = 25;
          for (int sym2__ = 1; sym2__ <= N; ++sym2__) {
            current_statement__ = 25;
            stan::model::assign(Y, Y_flat__[(pos__ - 1)],
              "assigning variable Y", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 25;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 26;
      stan::math::validate_non_negative_index("source", "N", N);
      current_statement__ = 27;
      context__.validate_dims("data initialization", "source", "int",
        std::vector<size_t>{static_cast<size_t>(N)});
      source = std::vector<int>(N, std::numeric_limits<int>::min());
      current_statement__ = 27;
      source = context__.vals_i("source");
      current_statement__ = 27;
      stan::math::check_greater_or_equal(function__, "source", source, 1);
      current_statement__ = 27;
      stan::math::check_less_or_equal(function__, "source", source, K);
      current_statement__ = 28;
      context__.validate_dims("data initialization", "shape_sig", "double",
        std::vector<size_t>{});
      shape_sig = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 28;
      shape_sig = context__.vals_r("shape_sig")[(1 - 1)];
      current_statement__ = 28;
      stan::math::check_greater_or_equal(function__, "shape_sig", shape_sig,
        1);
      current_statement__ = 29;
      stan::math::validate_non_negative_index("mu_jk", "J", J);
      current_statement__ = 30;
      stan::math::validate_non_negative_index("mu_jk", "K", K);
      current_statement__ = 31;
      stan::math::validate_non_negative_index("Omega", "K", K);
      current_statement__ = 32;
      stan::math::validate_non_negative_index("Omega", "J", J);
      current_statement__ = 32;
      stan::math::validate_non_negative_index("Omega", "J", J);
      current_statement__ = 33;
      stan::math::validate_non_negative_index("tau", "K", K);
      current_statement__ = 34;
      stan::math::validate_non_negative_index("Sigma", "K", K);
      current_statement__ = 35;
      stan::math::validate_non_negative_index("Sigma", "J", J);
      current_statement__ = 36;
      stan::math::validate_non_negative_index("Sigma", "J", J);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = (J * K) + (K * ((J * (J - 1)) / 2)) + K;
  }
  inline std::string model_name() const final {
    return "model_raw_source";
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
      "model_raw_source_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<local_scalar_t__,-1,-1> mu_jk =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(J, K, DUMMY_VAR__);
      current_statement__ = 1;
      mu_jk = in__.template read<Eigen::Matrix<local_scalar_t__,-1,-1>>(J, K);
      std::vector<Eigen::Matrix<local_scalar_t__,-1,-1>> Omega =
        std::vector<Eigen::Matrix<local_scalar_t__,-1,-1>>(K,
          Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(J, J, DUMMY_VAR__));
      current_statement__ = 2;
      Omega = in__.template read_constrain_corr_matrix<
                std::vector<Eigen::Matrix<local_scalar_t__,-1,-1>>,
                jacobian__>(lp__, K, J);
      Eigen::Matrix<local_scalar_t__,-1,1> tau =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(K, DUMMY_VAR__);
      current_statement__ = 3;
      tau = in__.template read_constrain_lb<
              Eigen::Matrix<local_scalar_t__,-1,1>, jacobian__>(0, lp__, K);
      std::vector<Eigen::Matrix<local_scalar_t__,-1,-1>> Sigma =
        std::vector<Eigen::Matrix<local_scalar_t__,-1,-1>>(K,
          Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(J, J, DUMMY_VAR__));
      current_statement__ = 7;
      for (int k = 1; k <= K; ++k) {
        current_statement__ = 5;
        stan::model::assign(Sigma,
          stan::math::quad_form_diag(
            stan::model::rvalue(Omega, "Omega", stan::model::index_uni(k)),
            stan::math::multiply(
              stan::model::rvalue(tau, "tau", stan::model::index_uni(k)),
              stan::math::ones_vector(J))), "assigning variable Sigma",
          stan::model::index_uni(k));
      }
      {
        current_statement__ = 8;
        lp_accum__.add(stan::math::cauchy_lpdf<propto__>(tau, 0, 2.5));
        current_statement__ = 13;
        for (int j = 1; j <= J; ++j) {
          current_statement__ = 11;
          for (int k = 1; k <= K; ++k) {
            current_statement__ = 9;
            lp_accum__.add(stan::math::normal_lpdf<propto__>(
                             stan::model::rvalue(mu_jk, "mu_jk",
                               stan::model::index_uni(j),
                               stan::model::index_uni(k)), 0, 100));
          }
        }
        current_statement__ = 16;
        for (int k = 1; k <= K; ++k) {
          current_statement__ = 14;
          lp_accum__.add(stan::math::lkj_corr_lpdf<propto__>(
                           stan::model::rvalue(Omega, "Omega",
                             stan::model::index_uni(k)), shape_sig));
        }
        current_statement__ = 19;
        for (int i = 1; i <= N; ++i) {
          current_statement__ = 17;
          lp_accum__.add(stan::math::multi_normal_lpdf<propto__>(
                           stan::model::rvalue(Y, "Y",
                             stan::model::index_uni(i)),
                           stan::model::rvalue(mu_jk, "mu_jk",
                             stan::model::index_omni(),
                             stan::model::index_uni(
                               stan::model::rvalue(source, "source",
                                 stan::model::index_uni(i)))),
                           stan::model::rvalue(Sigma, "Sigma",
                             stan::model::index_uni(
                               stan::model::rvalue(source, "source",
                                 stan::model::index_uni(i))))));
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
      "model_raw_source_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<double,-1,-1> mu_jk =
        Eigen::Matrix<double,-1,-1>::Constant(J, K,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 1;
      mu_jk = in__.template read<Eigen::Matrix<local_scalar_t__,-1,-1>>(J, K);
      std::vector<Eigen::Matrix<double,-1,-1>> Omega =
        std::vector<Eigen::Matrix<double,-1,-1>>(K,
          Eigen::Matrix<double,-1,-1>::Constant(J, J,
            std::numeric_limits<double>::quiet_NaN()));
      current_statement__ = 2;
      Omega = in__.template read_constrain_corr_matrix<
                std::vector<Eigen::Matrix<local_scalar_t__,-1,-1>>,
                jacobian__>(lp__, K, J);
      Eigen::Matrix<double,-1,1> tau =
        Eigen::Matrix<double,-1,1>::Constant(K,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 3;
      tau = in__.template read_constrain_lb<
              Eigen::Matrix<local_scalar_t__,-1,1>, jacobian__>(0, lp__, K);
      std::vector<Eigen::Matrix<double,-1,-1>> Sigma =
        std::vector<Eigen::Matrix<double,-1,-1>>(K,
          Eigen::Matrix<double,-1,-1>::Constant(J, J,
            std::numeric_limits<double>::quiet_NaN()));
      out__.write(mu_jk);
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= J; ++sym2__) {
          for (int sym3__ = 1; sym3__ <= K; ++sym3__) {
            out__.write(stan::model::rvalue(Omega, "Omega",
                          stan::model::index_uni(sym3__),
                          stan::model::index_uni(sym2__),
                          stan::model::index_uni(sym1__)));
          }
        }
      }
      out__.write(tau);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      current_statement__ = 7;
      for (int k = 1; k <= K; ++k) {
        current_statement__ = 5;
        stan::model::assign(Sigma,
          stan::math::quad_form_diag(
            stan::model::rvalue(Omega, "Omega", stan::model::index_uni(k)),
            stan::math::multiply(
              stan::model::rvalue(tau, "tau", stan::model::index_uni(k)),
              stan::math::ones_vector(J))), "assigning variable Sigma",
          stan::model::index_uni(k));
      }
      if (emit_transformed_parameters__) {
        for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
          for (int sym2__ = 1; sym2__ <= J; ++sym2__) {
            for (int sym3__ = 1; sym3__ <= K; ++sym3__) {
              out__.write(stan::model::rvalue(Sigma, "Sigma",
                            stan::model::index_uni(sym3__),
                            stan::model::index_uni(sym2__),
                            stan::model::index_uni(sym1__)));
            }
          }
        }
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
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
      Eigen::Matrix<local_scalar_t__,-1,-1> mu_jk =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(J, K, DUMMY_VAR__);
      current_statement__ = 1;
      stan::model::assign(mu_jk,
        in__.read<Eigen::Matrix<local_scalar_t__,-1,-1>>(J, K),
        "assigning variable mu_jk");
      out__.write(mu_jk);
      std::vector<Eigen::Matrix<local_scalar_t__,-1,-1>> Omega =
        std::vector<Eigen::Matrix<local_scalar_t__,-1,-1>>(K,
          Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(J, J, DUMMY_VAR__));
      current_statement__ = 2;
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        current_statement__ = 2;
        for (int sym2__ = 1; sym2__ <= J; ++sym2__) {
          current_statement__ = 2;
          for (int sym3__ = 1; sym3__ <= K; ++sym3__) {
            current_statement__ = 2;
            stan::model::assign(Omega, in__.read<local_scalar_t__>(),
              "assigning variable Omega", stan::model::index_uni(sym3__),
              stan::model::index_uni(sym2__), stan::model::index_uni(sym1__));
          }
        }
      }
      out__.write_free_corr_matrix(Omega);
      Eigen::Matrix<local_scalar_t__,-1,1> tau =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(K, DUMMY_VAR__);
      current_statement__ = 3;
      stan::model::assign(tau,
        in__.read<Eigen::Matrix<local_scalar_t__,-1,1>>(K),
        "assigning variable tau");
      out__.write_free_lb(0, tau);
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
      context__.validate_dims("parameter initialization", "mu_jk", "double",
        std::vector<size_t>{static_cast<size_t>(J), static_cast<size_t>(K)});
      current_statement__ = 2;
      context__.validate_dims("parameter initialization", "Omega", "double",
        std::vector<size_t>{static_cast<size_t>(K), static_cast<size_t>(J),
          static_cast<size_t>(J)});
      current_statement__ = 3;
      context__.validate_dims("parameter initialization", "tau", "double",
        std::vector<size_t>{static_cast<size_t>(K)});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,-1> mu_jk =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(J, K, DUMMY_VAR__);
      {
        std::vector<local_scalar_t__> mu_jk_flat__;
        current_statement__ = 1;
        mu_jk_flat__ = context__.vals_r("mu_jk");
        current_statement__ = 1;
        pos__ = 1;
        current_statement__ = 1;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 1;
          for (int sym2__ = 1; sym2__ <= J; ++sym2__) {
            current_statement__ = 1;
            stan::model::assign(mu_jk, mu_jk_flat__[(pos__ - 1)],
              "assigning variable mu_jk", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 1;
            pos__ = (pos__ + 1);
          }
        }
      }
      out__.write(mu_jk);
      std::vector<Eigen::Matrix<local_scalar_t__,-1,-1>> Omega =
        std::vector<Eigen::Matrix<local_scalar_t__,-1,-1>>(K,
          Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(J, J, DUMMY_VAR__));
      {
        std::vector<local_scalar_t__> Omega_flat__;
        current_statement__ = 2;
        Omega_flat__ = context__.vals_r("Omega");
        current_statement__ = 2;
        pos__ = 1;
        current_statement__ = 2;
        for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
          current_statement__ = 2;
          for (int sym2__ = 1; sym2__ <= J; ++sym2__) {
            current_statement__ = 2;
            for (int sym3__ = 1; sym3__ <= K; ++sym3__) {
              current_statement__ = 2;
              stan::model::assign(Omega, Omega_flat__[(pos__ - 1)],
                "assigning variable Omega", stan::model::index_uni(sym3__),
                stan::model::index_uni(sym2__),
                stan::model::index_uni(sym1__));
              current_statement__ = 2;
              pos__ = (pos__ + 1);
            }
          }
        }
      }
      out__.write_free_corr_matrix(Omega);
      Eigen::Matrix<local_scalar_t__,-1,1> tau =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(K, DUMMY_VAR__);
      {
        std::vector<local_scalar_t__> tau_flat__;
        current_statement__ = 3;
        tau_flat__ = context__.vals_r("tau");
        current_statement__ = 3;
        pos__ = 1;
        current_statement__ = 3;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 3;
          stan::model::assign(tau, tau_flat__[(pos__ - 1)],
            "assigning variable tau", stan::model::index_uni(sym1__));
          current_statement__ = 3;
          pos__ = (pos__ + 1);
        }
      }
      out__.write_free_lb(0, tau);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"mu_jk", "Omega", "tau"};
    if (emit_transformed_parameters__) {
      std::vector<std::string> temp{"Sigma"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{static_cast<
                                                                    size_t>(J),
                                                 static_cast<size_t>(K)},
                std::vector<size_t>{static_cast<size_t>(K),
                  static_cast<size_t>(J), static_cast<size_t>(J)},
                std::vector<size_t>{static_cast<size_t>(K)}};
    if (emit_transformed_parameters__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(K),
               static_cast<size_t>(J), static_cast<size_t>(J)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= J; ++sym2__) {
        param_names__.emplace_back(std::string() + "mu_jk" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= J; ++sym2__) {
        for (int sym3__ = 1; sym3__ <= K; ++sym3__) {
          param_names__.emplace_back(std::string() + "Omega" + '.' +
            std::to_string(sym3__) + '.' + std::to_string(sym2__) + '.' +
            std::to_string(sym1__));
        }
      }
    }
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      param_names__.emplace_back(std::string() + "tau" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= J; ++sym2__) {
          for (int sym3__ = 1; sym3__ <= K; ++sym3__) {
            param_names__.emplace_back(std::string() + "Sigma" + '.' +
              std::to_string(sym3__) + '.' + std::to_string(sym2__) + '.' +
              std::to_string(sym1__));
          }
        }
      }
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= J; ++sym2__) {
        param_names__.emplace_back(std::string() + "mu_jk" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    for (int sym1__ = 1; sym1__ <= ((J * (J - 1)) / 2); ++sym1__) {
      for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
        param_names__.emplace_back(std::string() + "Omega" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      param_names__.emplace_back(std::string() + "tau" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= J; ++sym2__) {
          for (int sym3__ = 1; sym3__ <= K; ++sym3__) {
            param_names__.emplace_back(std::string() + "Sigma" + '.' +
              std::to_string(sym3__) + '.' + std::to_string(sym2__) + '.' +
              std::to_string(sym1__));
          }
        }
      }
    }
    if (emit_generated_quantities__) {}
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"mu_jk\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(J) + ",\"cols\":" + std::to_string(K) + "},\"block\":\"parameters\"},{\"name\":\"Omega\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(K) + ",\"element_type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(J) + ",\"cols\":" + std::to_string(J) + "}},\"block\":\"parameters\"},{\"name\":\"tau\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(K) + "},\"block\":\"parameters\"},{\"name\":\"Sigma\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(K) + ",\"element_type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(J) + ",\"cols\":" + std::to_string(J) + "}},\"block\":\"transformed_parameters\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"mu_jk\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(J) + ",\"cols\":" + std::to_string(K) + "},\"block\":\"parameters\"},{\"name\":\"Omega\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(K) + ",\"element_type\":{\"name\":\"vector\",\"length\":" + std::to_string(((J * (J - 1)) /2)) + "}},\"block\":\"parameters\"},{\"name\":\"tau\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(K) + "},\"block\":\"parameters\"},{\"name\":\"Sigma\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(K) + ",\"element_type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(J) + ",\"cols\":" + std::to_string(J) + "}},\"block\":\"transformed_parameters\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = (((J * K) + ((K * J) * J)) + K);
    const size_t num_transformed = emit_transformed_parameters * (((K * J) *
      J));
    const size_t num_gen_quantities = emit_generated_quantities * (0);
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
    const size_t num_params__ = (((J * K) + ((K * J) * J)) + K);
    const size_t num_transformed = emit_transformed_parameters * (((K * J) *
      J));
    const size_t num_gen_quantities = emit_generated_quantities * (0);
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
using stan_model = model_raw_source_namespace::model_raw_source;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_raw_source_namespace::profiles__;
}
#endif
#endif
