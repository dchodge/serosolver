// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/serosolver.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// titre_data_fast_individual_base
void titre_data_fast_individual_base(NumericVector& predicted_titres, const NumericVector& theta, const List& infection_info, const List& vaccination_info, const List& setup_data, const List& indexing, const List& antigenic_maps, const List& other_pars);
RcppExport SEXP _serosolver_titre_data_fast_individual_base(SEXP predicted_titresSEXP, SEXP thetaSEXP, SEXP infection_infoSEXP, SEXP vaccination_infoSEXP, SEXP setup_dataSEXP, SEXP indexingSEXP, SEXP antigenic_mapsSEXP, SEXP other_parsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type infection_info(infection_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type vaccination_info(vaccination_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type setup_data(setup_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type indexing(indexingSEXP);
    Rcpp::traits::input_parameter< const List& >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const List& >::type other_pars(other_parsSEXP);
    titre_data_fast_individual_base(predicted_titres, theta, infection_info, vaccination_info, setup_data, indexing, antigenic_maps, other_pars);
    return R_NilValue;
END_RCPP
}
// titre_data_fast_individual_wane2
void titre_data_fast_individual_wane2(NumericVector& predicted_titres, NumericVector& theta, const List& infection_info, const List& vaccination_info, const List& setup_data, const List& indexing, const List& antigenic_maps, const List& other_pars);
RcppExport SEXP _serosolver_titre_data_fast_individual_wane2(SEXP predicted_titresSEXP, SEXP thetaSEXP, SEXP infection_infoSEXP, SEXP vaccination_infoSEXP, SEXP setup_dataSEXP, SEXP indexingSEXP, SEXP antigenic_mapsSEXP, SEXP other_parsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type infection_info(infection_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type vaccination_info(vaccination_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type setup_data(setup_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type indexing(indexingSEXP);
    Rcpp::traits::input_parameter< const List& >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const List& >::type other_pars(other_parsSEXP);
    titre_data_fast_individual_wane2(predicted_titres, theta, infection_info, vaccination_info, setup_data, indexing, antigenic_maps, other_pars);
    return R_NilValue;
END_RCPP
}
// titre_data_fast_individual_titredep
void titre_data_fast_individual_titredep(NumericVector& predicted_titres, NumericVector& theta, const List& infection_info, const List& vaccination_info, const List& setup_data, const List& indexing, const List& antigenic_maps, const List& other_pars);
RcppExport SEXP _serosolver_titre_data_fast_individual_titredep(SEXP predicted_titresSEXP, SEXP thetaSEXP, SEXP infection_infoSEXP, SEXP vaccination_infoSEXP, SEXP setup_dataSEXP, SEXP indexingSEXP, SEXP antigenic_mapsSEXP, SEXP other_parsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type infection_info(infection_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type vaccination_info(vaccination_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type setup_data(setup_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type indexing(indexingSEXP);
    Rcpp::traits::input_parameter< const List& >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const List& >::type other_pars(other_parsSEXP);
    titre_data_fast_individual_titredep(predicted_titres, theta, infection_info, vaccination_info, setup_data, indexing, antigenic_maps, other_pars);
    return R_NilValue;
END_RCPP
}
// titre_data_fast_individual_strain_dependent
void titre_data_fast_individual_strain_dependent(NumericVector& predicted_titres, NumericVector& theta, const List& infection_info, const List& vaccination_info, const List& setup_data, const List& indexing, const List& antigenic_maps, const List& other_pars);
RcppExport SEXP _serosolver_titre_data_fast_individual_strain_dependent(SEXP predicted_titresSEXP, SEXP thetaSEXP, SEXP infection_infoSEXP, SEXP vaccination_infoSEXP, SEXP setup_dataSEXP, SEXP indexingSEXP, SEXP antigenic_mapsSEXP, SEXP other_parsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type infection_info(infection_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type vaccination_info(vaccination_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type setup_data(setup_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type indexing(indexingSEXP);
    Rcpp::traits::input_parameter< const List& >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const List& >::type other_pars(other_parsSEXP);
    titre_data_fast_individual_strain_dependent(predicted_titres, theta, infection_info, vaccination_info, setup_data, indexing, antigenic_maps, other_pars);
    return R_NilValue;
END_RCPP
}
// fun_cpp
void fun_cpp(int x);
RcppExport SEXP _serosolver_fun_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    fun_cpp(x);
    return R_NilValue;
END_RCPP
}
// timesTwo
void timesTwo(NumericVector x);
RcppExport SEXP _serosolver_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    timesTwo(x);
    return R_NilValue;
END_RCPP
}
// create_ptr
Rcpp::XPtr<FunctionPointer> create_ptr();
RcppExport SEXP _serosolver_create_ptr() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(create_ptr());
    return rcpp_result_gen;
END_RCPP
}
// invokeCallback
void invokeCallback(XPtr<FunctionPointer> callback);
RcppExport SEXP _serosolver_invokeCallback(SEXP callbackSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<FunctionPointer> >::type callback(callbackSEXP);
    invokeCallback(callback);
    return R_NilValue;
END_RCPP
}
// get_callback
InternalFunction get_callback();
RcppExport SEXP _serosolver_get_callback() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(get_callback());
    return rcpp_result_gen;
END_RCPP
}
// subset_nullable_vector
NumericVector subset_nullable_vector(const Nullable<NumericVector>& x, int index1, int index2);
RcppExport SEXP _serosolver_subset_nullable_vector(SEXP xSEXP, SEXP index1SEXP, SEXP index2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type index1(index1SEXP);
    Rcpp::traits::input_parameter< int >::type index2(index2SEXP);
    rcpp_result_gen = Rcpp::wrap(subset_nullable_vector(x, index1, index2));
    return rcpp_result_gen;
END_RCPP
}
// sum_likelihoods
NumericVector sum_likelihoods(NumericVector liks, IntegerVector indices, int n_indivs);
RcppExport SEXP _serosolver_sum_likelihoods(SEXP liksSEXP, SEXP indicesSEXP, SEXP n_indivsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type liks(liksSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< int >::type n_indivs(n_indivsSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_likelihoods(liks, indices, n_indivs));
    return rcpp_result_gen;
END_RCPP
}
// create_cross_reactivity_vector
NumericVector create_cross_reactivity_vector(NumericVector x, double sigma);
RcppExport SEXP _serosolver_create_cross_reactivity_vector(SEXP xSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(create_cross_reactivity_vector(x, sigma));
    return rcpp_result_gen;
END_RCPP
}
// sum_buckets
NumericVector sum_buckets(NumericVector a, NumericVector buckets);
RcppExport SEXP _serosolver_sum_buckets(SEXP aSEXP, SEXP bucketsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type buckets(bucketsSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_buckets(a, buckets));
    return rcpp_result_gen;
END_RCPP
}
// sum_infections_by_group
IntegerMatrix sum_infections_by_group(IntegerMatrix inf_hist, IntegerVector group_ids_vec, int n_groups);
RcppExport SEXP _serosolver_sum_infections_by_group(SEXP inf_histSEXP, SEXP group_ids_vecSEXP, SEXP n_groupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type inf_hist(inf_histSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group_ids_vec(group_ids_vecSEXP);
    Rcpp::traits::input_parameter< int >::type n_groups(n_groupsSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_infections_by_group(inf_hist, group_ids_vec, n_groups));
    return rcpp_result_gen;
END_RCPP
}
// add_measurement_shifts
void add_measurement_shifts(NumericVector& predicted_titres, const NumericVector& to_add, const int& start_index_in_data, const int& end_index_in_data);
RcppExport SEXP _serosolver_add_measurement_shifts(SEXP predicted_titresSEXP, SEXP to_addSEXP, SEXP start_index_in_dataSEXP, SEXP end_index_in_dataSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type to_add(to_addSEXP);
    Rcpp::traits::input_parameter< const int& >::type start_index_in_data(start_index_in_dataSEXP);
    Rcpp::traits::input_parameter< const int& >::type end_index_in_data(end_index_in_dataSEXP);
    add_measurement_shifts(predicted_titres, to_add, start_index_in_data, end_index_in_data);
    return R_NilValue;
END_RCPP
}
// titre_data_fast
NumericVector titre_data_fast(NumericVector theta, const IntegerMatrix infection_history_mat, const List vaccination_hist_info, const Function ab_kin_func, const NumericVector circulation_times, const IntegerVector circulation_times_indices, const NumericVector sample_times, const IntegerVector rows_per_indiv_in_samples, const IntegerVector cum_nrows_per_individual_in_data, const IntegerVector nrows_per_blood_sample, const IntegerVector measurement_strain_indices, const List antigenic_maps, const NumericVector antigenic_distances, const List other_pars, bool boost_before_infection);
RcppExport SEXP _serosolver_titre_data_fast(SEXP thetaSEXP, SEXP infection_history_matSEXP, SEXP vaccination_hist_infoSEXP, SEXP ab_kin_funcSEXP, SEXP circulation_timesSEXP, SEXP circulation_times_indicesSEXP, SEXP sample_timesSEXP, SEXP rows_per_indiv_in_samplesSEXP, SEXP cum_nrows_per_individual_in_dataSEXP, SEXP nrows_per_blood_sampleSEXP, SEXP measurement_strain_indicesSEXP, SEXP antigenic_mapsSEXP, SEXP antigenic_distancesSEXP, SEXP other_parsSEXP, SEXP boost_before_infectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type infection_history_mat(infection_history_matSEXP);
    Rcpp::traits::input_parameter< const List >::type vaccination_hist_info(vaccination_hist_infoSEXP);
    Rcpp::traits::input_parameter< const Function >::type ab_kin_func(ab_kin_funcSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type circulation_times(circulation_timesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type circulation_times_indices(circulation_times_indicesSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type sample_times(sample_timesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type rows_per_indiv_in_samples(rows_per_indiv_in_samplesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type cum_nrows_per_individual_in_data(cum_nrows_per_individual_in_dataSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nrows_per_blood_sample(nrows_per_blood_sampleSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type measurement_strain_indices(measurement_strain_indicesSEXP);
    Rcpp::traits::input_parameter< const List >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type antigenic_distances(antigenic_distancesSEXP);
    Rcpp::traits::input_parameter< const List >::type other_pars(other_parsSEXP);
    Rcpp::traits::input_parameter< bool >::type boost_before_infection(boost_before_infectionSEXP);
    rcpp_result_gen = Rcpp::wrap(titre_data_fast(theta, infection_history_mat, vaccination_hist_info, ab_kin_func, circulation_times, circulation_times_indices, sample_times, rows_per_indiv_in_samples, cum_nrows_per_individual_in_data, nrows_per_blood_sample, measurement_strain_indices, antigenic_maps, antigenic_distances, other_pars, boost_before_infection));
    return rcpp_result_gen;
END_RCPP
}
// inf_mat_prior_cpp
double inf_mat_prior_cpp(const IntegerMatrix& infection_history, const IntegerVector& n_alive, double alpha, double beta);
RcppExport SEXP _serosolver_inf_mat_prior_cpp(SEXP infection_historySEXP, SEXP n_aliveSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type infection_history(infection_historySEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n_alive(n_aliveSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(inf_mat_prior_cpp(infection_history, n_alive, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// inf_mat_prior_cpp_vector
double inf_mat_prior_cpp_vector(const IntegerMatrix& infection_history, const IntegerVector& n_alive, const NumericVector& alphas, const NumericVector& betas);
RcppExport SEXP _serosolver_inf_mat_prior_cpp_vector(SEXP infection_historySEXP, SEXP n_aliveSEXP, SEXP alphasSEXP, SEXP betasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type infection_history(infection_historySEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n_alive(n_aliveSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type betas(betasSEXP);
    rcpp_result_gen = Rcpp::wrap(inf_mat_prior_cpp_vector(infection_history, n_alive, alphas, betas));
    return rcpp_result_gen;
END_RCPP
}
// inf_mat_prior_group_cpp
double inf_mat_prior_group_cpp(const IntegerMatrix& n_infections, const IntegerMatrix& n_alive, double alpha, double beta);
RcppExport SEXP _serosolver_inf_mat_prior_group_cpp(SEXP n_infectionsSEXP, SEXP n_aliveSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type n_infections(n_infectionsSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type n_alive(n_aliveSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(inf_mat_prior_group_cpp(n_infections, n_alive, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// inf_mat_prior_group_cpp_vector
double inf_mat_prior_group_cpp_vector(const IntegerMatrix& n_infections, const IntegerMatrix& n_alive, const NumericVector& alphas, const NumericVector& betas);
RcppExport SEXP _serosolver_inf_mat_prior_group_cpp_vector(SEXP n_infectionsSEXP, SEXP n_aliveSEXP, SEXP alphasSEXP, SEXP betasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type n_infections(n_infectionsSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type n_alive(n_aliveSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type betas(betasSEXP);
    rcpp_result_gen = Rcpp::wrap(inf_mat_prior_group_cpp_vector(n_infections, n_alive, alphas, betas));
    return rcpp_result_gen;
END_RCPP
}
// inf_mat_prior_total_group_cpp
double inf_mat_prior_total_group_cpp(const IntegerVector& n_infections_group, const IntegerVector& n_alive_group, double alpha, double beta);
RcppExport SEXP _serosolver_inf_mat_prior_total_group_cpp(SEXP n_infections_groupSEXP, SEXP n_alive_groupSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type n_infections_group(n_infections_groupSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n_alive_group(n_alive_groupSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(inf_mat_prior_total_group_cpp(n_infections_group, n_alive_group, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// likelihood_func_fast
NumericVector likelihood_func_fast(const NumericVector& theta, const NumericVector& obs, const NumericVector& predicted_titres);
RcppExport SEXP _serosolver_likelihood_func_fast(SEXP thetaSEXP, SEXP obsSEXP, SEXP predicted_titresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type predicted_titres(predicted_titresSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_func_fast(theta, obs, predicted_titres));
    return rcpp_result_gen;
END_RCPP
}
// inf_hist_prop_prior_v3
arma::mat inf_hist_prop_prior_v3(arma::mat infection_history_mat, const IntegerVector& sampled_indivs, const IntegerVector& age_mask, const IntegerVector& strain_mask, const IntegerVector& move_sizes, const IntegerVector& n_infs, double alpha, double beta, const NumericVector& rand_ns, const double& swap_propn);
RcppExport SEXP _serosolver_inf_hist_prop_prior_v3(SEXP infection_history_matSEXP, SEXP sampled_indivsSEXP, SEXP age_maskSEXP, SEXP strain_maskSEXP, SEXP move_sizesSEXP, SEXP n_infsSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP rand_nsSEXP, SEXP swap_propnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type infection_history_mat(infection_history_matSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type sampled_indivs(sampled_indivsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type age_mask(age_maskSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type strain_mask(strain_maskSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type move_sizes(move_sizesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n_infs(n_infsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rand_ns(rand_nsSEXP);
    Rcpp::traits::input_parameter< const double& >::type swap_propn(swap_propnSEXP);
    rcpp_result_gen = Rcpp::wrap(inf_hist_prop_prior_v3(infection_history_mat, sampled_indivs, age_mask, strain_mask, move_sizes, n_infs, alpha, beta, rand_ns, swap_propn));
    return rcpp_result_gen;
END_RCPP
}
// inf_hist_prop_prior_v2_and_v4
List inf_hist_prop_prior_v2_and_v4(NumericVector theta, const IntegerMatrix infection_history_mat, List vaccination_hist_info, const Function ab_kin_func, const NumericVector old_probs_1, const IntegerVector sampled_indivs, const IntegerVector n_years_samp_vec, const IntegerVector age_mask, const IntegerVector strain_mask, const IntegerMatrix n_alive, IntegerMatrix n_infections, IntegerVector n_infected_group, const NumericMatrix prior_lookup, const double swap_propn, const int swap_distance, const bool propose_from_prior, const double alpha, const double beta, const NumericVector circulation_times, const IntegerVector circulation_times_indices, const NumericVector sample_times, const IntegerVector rows_per_indiv_in_samples, const IntegerVector cum_nrows_per_individual_in_data, const IntegerVector cum_nrows_per_individual_in_repeat_data, const IntegerVector nrows_per_blood_sample, const IntegerVector group_id_vec, const IntegerVector measurement_strain_indices, const List antigenic_maps, const NumericVector antigenic_distances, const NumericVector data, const NumericVector repeat_data, const IntegerVector repeat_indices, const NumericVector titre_shifts, IntegerVector proposal_iter, IntegerVector accepted_iter, IntegerVector proposal_swap, IntegerVector accepted_swap, IntegerMatrix overall_swap_proposals, IntegerMatrix overall_add_proposals, const NumericVector time_sample_probs, const List other_pars, const IntegerVector total_alive, const double temp, bool solve_likelihood);
RcppExport SEXP _serosolver_inf_hist_prop_prior_v2_and_v4(SEXP thetaSEXP, SEXP infection_history_matSEXP, SEXP vaccination_hist_infoSEXP, SEXP ab_kin_funcSEXP, SEXP old_probs_1SEXP, SEXP sampled_indivsSEXP, SEXP n_years_samp_vecSEXP, SEXP age_maskSEXP, SEXP strain_maskSEXP, SEXP n_aliveSEXP, SEXP n_infectionsSEXP, SEXP n_infected_groupSEXP, SEXP prior_lookupSEXP, SEXP swap_propnSEXP, SEXP swap_distanceSEXP, SEXP propose_from_priorSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP circulation_timesSEXP, SEXP circulation_times_indicesSEXP, SEXP sample_timesSEXP, SEXP rows_per_indiv_in_samplesSEXP, SEXP cum_nrows_per_individual_in_dataSEXP, SEXP cum_nrows_per_individual_in_repeat_dataSEXP, SEXP nrows_per_blood_sampleSEXP, SEXP group_id_vecSEXP, SEXP measurement_strain_indicesSEXP, SEXP antigenic_mapsSEXP, SEXP antigenic_distancesSEXP, SEXP dataSEXP, SEXP repeat_dataSEXP, SEXP repeat_indicesSEXP, SEXP titre_shiftsSEXP, SEXP proposal_iterSEXP, SEXP accepted_iterSEXP, SEXP proposal_swapSEXP, SEXP accepted_swapSEXP, SEXP overall_swap_proposalsSEXP, SEXP overall_add_proposalsSEXP, SEXP time_sample_probsSEXP, SEXP other_parsSEXP, SEXP total_aliveSEXP, SEXP tempSEXP, SEXP solve_likelihoodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type infection_history_mat(infection_history_matSEXP);
    Rcpp::traits::input_parameter< List >::type vaccination_hist_info(vaccination_hist_infoSEXP);
    Rcpp::traits::input_parameter< const Function >::type ab_kin_func(ab_kin_funcSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type old_probs_1(old_probs_1SEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type sampled_indivs(sampled_indivsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type n_years_samp_vec(n_years_samp_vecSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type age_mask(age_maskSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type strain_mask(strain_maskSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type n_alive(n_aliveSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type n_infections(n_infectionsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_infected_group(n_infected_groupSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type prior_lookup(prior_lookupSEXP);
    Rcpp::traits::input_parameter< const double >::type swap_propn(swap_propnSEXP);
    Rcpp::traits::input_parameter< const int >::type swap_distance(swap_distanceSEXP);
    Rcpp::traits::input_parameter< const bool >::type propose_from_prior(propose_from_priorSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type circulation_times(circulation_timesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type circulation_times_indices(circulation_times_indicesSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type sample_times(sample_timesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type rows_per_indiv_in_samples(rows_per_indiv_in_samplesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type cum_nrows_per_individual_in_data(cum_nrows_per_individual_in_dataSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type cum_nrows_per_individual_in_repeat_data(cum_nrows_per_individual_in_repeat_dataSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type nrows_per_blood_sample(nrows_per_blood_sampleSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type group_id_vec(group_id_vecSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type measurement_strain_indices(measurement_strain_indicesSEXP);
    Rcpp::traits::input_parameter< const List >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type antigenic_distances(antigenic_distancesSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type repeat_data(repeat_dataSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type repeat_indices(repeat_indicesSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type titre_shifts(titre_shiftsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type proposal_iter(proposal_iterSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type accepted_iter(accepted_iterSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type proposal_swap(proposal_swapSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type accepted_swap(accepted_swapSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type overall_swap_proposals(overall_swap_proposalsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type overall_add_proposals(overall_add_proposalsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type time_sample_probs(time_sample_probsSEXP);
    Rcpp::traits::input_parameter< const List >::type other_pars(other_parsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type total_alive(total_aliveSEXP);
    Rcpp::traits::input_parameter< const double >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< bool >::type solve_likelihood(solve_likelihoodSEXP);
    rcpp_result_gen = Rcpp::wrap(inf_hist_prop_prior_v2_and_v4(theta, infection_history_mat, vaccination_hist_info, ab_kin_func, old_probs_1, sampled_indivs, n_years_samp_vec, age_mask, strain_mask, n_alive, n_infections, n_infected_group, prior_lookup, swap_propn, swap_distance, propose_from_prior, alpha, beta, circulation_times, circulation_times_indices, sample_times, rows_per_indiv_in_samples, cum_nrows_per_individual_in_data, cum_nrows_per_individual_in_repeat_data, nrows_per_blood_sample, group_id_vec, measurement_strain_indices, antigenic_maps, antigenic_distances, data, repeat_data, repeat_indices, titre_shifts, proposal_iter, accepted_iter, proposal_swap, accepted_swap, overall_swap_proposals, overall_add_proposals, time_sample_probs, other_pars, total_alive, temp, solve_likelihood));
    return rcpp_result_gen;
END_RCPP
}
// wane_function
double wane_function(NumericVector theta, double time_infected, double wane);
RcppExport SEXP _serosolver_wane_function(SEXP thetaSEXP, SEXP time_infectedSEXP, SEXP waneSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type time_infected(time_infectedSEXP);
    Rcpp::traits::input_parameter< double >::type wane(waneSEXP);
    rcpp_result_gen = Rcpp::wrap(wane_function(theta, time_infected, wane));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_serosolver_titre_data_fast_individual_base", (DL_FUNC) &_serosolver_titre_data_fast_individual_base, 8},
    {"_serosolver_titre_data_fast_individual_wane2", (DL_FUNC) &_serosolver_titre_data_fast_individual_wane2, 8},
    {"_serosolver_titre_data_fast_individual_titredep", (DL_FUNC) &_serosolver_titre_data_fast_individual_titredep, 8},
    {"_serosolver_titre_data_fast_individual_strain_dependent", (DL_FUNC) &_serosolver_titre_data_fast_individual_strain_dependent, 8},
    {"_serosolver_fun_cpp", (DL_FUNC) &_serosolver_fun_cpp, 1},
    {"_serosolver_timesTwo", (DL_FUNC) &_serosolver_timesTwo, 1},
    {"_serosolver_create_ptr", (DL_FUNC) &_serosolver_create_ptr, 0},
    {"_serosolver_invokeCallback", (DL_FUNC) &_serosolver_invokeCallback, 1},
    {"_serosolver_get_callback", (DL_FUNC) &_serosolver_get_callback, 0},
    {"_serosolver_subset_nullable_vector", (DL_FUNC) &_serosolver_subset_nullable_vector, 3},
    {"_serosolver_sum_likelihoods", (DL_FUNC) &_serosolver_sum_likelihoods, 3},
    {"_serosolver_create_cross_reactivity_vector", (DL_FUNC) &_serosolver_create_cross_reactivity_vector, 2},
    {"_serosolver_sum_buckets", (DL_FUNC) &_serosolver_sum_buckets, 2},
    {"_serosolver_sum_infections_by_group", (DL_FUNC) &_serosolver_sum_infections_by_group, 3},
    {"_serosolver_add_measurement_shifts", (DL_FUNC) &_serosolver_add_measurement_shifts, 4},
    {"_serosolver_titre_data_fast", (DL_FUNC) &_serosolver_titre_data_fast, 15},
    {"_serosolver_inf_mat_prior_cpp", (DL_FUNC) &_serosolver_inf_mat_prior_cpp, 4},
    {"_serosolver_inf_mat_prior_cpp_vector", (DL_FUNC) &_serosolver_inf_mat_prior_cpp_vector, 4},
    {"_serosolver_inf_mat_prior_group_cpp", (DL_FUNC) &_serosolver_inf_mat_prior_group_cpp, 4},
    {"_serosolver_inf_mat_prior_group_cpp_vector", (DL_FUNC) &_serosolver_inf_mat_prior_group_cpp_vector, 4},
    {"_serosolver_inf_mat_prior_total_group_cpp", (DL_FUNC) &_serosolver_inf_mat_prior_total_group_cpp, 4},
    {"_serosolver_likelihood_func_fast", (DL_FUNC) &_serosolver_likelihood_func_fast, 3},
    {"_serosolver_inf_hist_prop_prior_v3", (DL_FUNC) &_serosolver_inf_hist_prop_prior_v3, 10},
    {"_serosolver_inf_hist_prop_prior_v2_and_v4", (DL_FUNC) &_serosolver_inf_hist_prop_prior_v2_and_v4, 44},
    {"_serosolver_wane_function", (DL_FUNC) &_serosolver_wane_function, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_serosolver(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
