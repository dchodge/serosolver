// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// infection_model_indiv
NumericVector infection_model_indiv(NumericVector theta, IntegerVector infectionHistory, NumericVector infectionTimes, IntegerVector infectionMapIndices, double samplingTime, IntegerVector measurementMapIndices, NumericVector antigenicMapLong, NumericVector antigenicMapShort, int numberStrains);
RcppExport SEXP _serosolver_infection_model_indiv(SEXP thetaSEXP, SEXP infectionHistorySEXP, SEXP infectionTimesSEXP, SEXP infectionMapIndicesSEXP, SEXP samplingTimeSEXP, SEXP measurementMapIndicesSEXP, SEXP antigenicMapLongSEXP, SEXP antigenicMapShortSEXP, SEXP numberStrainsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type infectionHistory(infectionHistorySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type infectionTimes(infectionTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type infectionMapIndices(infectionMapIndicesSEXP);
    Rcpp::traits::input_parameter< double >::type samplingTime(samplingTimeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type measurementMapIndices(measurementMapIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapLong(antigenicMapLongSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapShort(antigenicMapShortSEXP);
    Rcpp::traits::input_parameter< int >::type numberStrains(numberStrainsSEXP);
    rcpp_result_gen = Rcpp::wrap(infection_model_indiv(theta, infectionHistory, infectionTimes, infectionMapIndices, samplingTime, measurementMapIndices, antigenicMapLong, antigenicMapShort, numberStrains));
    return rcpp_result_gen;
END_RCPP
}
// titre_data_individual
NumericVector titre_data_individual(NumericVector theta, IntegerVector infectionHistory, NumericVector circulationTimes, IntegerVector circulationMapIndices, NumericVector samplingTimes, IntegerVector dataIndices, IntegerVector measuredMapIndices, NumericVector antigenicMapLong, NumericVector antigenicMapShort, int numberStrains);
RcppExport SEXP _serosolver_titre_data_individual(SEXP thetaSEXP, SEXP infectionHistorySEXP, SEXP circulationTimesSEXP, SEXP circulationMapIndicesSEXP, SEXP samplingTimesSEXP, SEXP dataIndicesSEXP, SEXP measuredMapIndicesSEXP, SEXP antigenicMapLongSEXP, SEXP antigenicMapShortSEXP, SEXP numberStrainsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type infectionHistory(infectionHistorySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type circulationTimes(circulationTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type circulationMapIndices(circulationMapIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samplingTimes(samplingTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dataIndices(dataIndicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type measuredMapIndices(measuredMapIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapLong(antigenicMapLongSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapShort(antigenicMapShortSEXP);
    Rcpp::traits::input_parameter< int >::type numberStrains(numberStrainsSEXP);
    rcpp_result_gen = Rcpp::wrap(titre_data_individual(theta, infectionHistory, circulationTimes, circulationMapIndices, samplingTimes, dataIndices, measuredMapIndices, antigenicMapLong, antigenicMapShort, numberStrains));
    return rcpp_result_gen;
END_RCPP
}
// titre_data_group
NumericVector titre_data_group(NumericVector theta, IntegerMatrix infectionHistories, NumericVector circulationTimes, IntegerVector circulationMapIndices, NumericVector samplingTimes, IntegerVector indicesTitreDataSample, IntegerVector indicesTitreDataOverall, IntegerVector indicesSamples, IntegerVector measuredMapIndices, NumericVector antigenicMapLong, NumericVector antigenicMapShort);
RcppExport SEXP _serosolver_titre_data_group(SEXP thetaSEXP, SEXP infectionHistoriesSEXP, SEXP circulationTimesSEXP, SEXP circulationMapIndicesSEXP, SEXP samplingTimesSEXP, SEXP indicesTitreDataSampleSEXP, SEXP indicesTitreDataOverallSEXP, SEXP indicesSamplesSEXP, SEXP measuredMapIndicesSEXP, SEXP antigenicMapLongSEXP, SEXP antigenicMapShortSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type infectionHistories(infectionHistoriesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type circulationTimes(circulationTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type circulationMapIndices(circulationMapIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samplingTimes(samplingTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesTitreDataSample(indicesTitreDataSampleSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesTitreDataOverall(indicesTitreDataOverallSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesSamples(indicesSamplesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type measuredMapIndices(measuredMapIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapLong(antigenicMapLongSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapShort(antigenicMapShortSEXP);
    rcpp_result_gen = Rcpp::wrap(titre_data_group(theta, infectionHistories, circulationTimes, circulationMapIndices, samplingTimes, indicesTitreDataSample, indicesTitreDataOverall, indicesSamples, measuredMapIndices, antigenicMapLong, antigenicMapShort));
    return rcpp_result_gen;
END_RCPP
}
// likelihood_titre
double likelihood_titre(NumericVector expected, NumericVector data, NumericVector theta);
RcppExport SEXP _serosolver_likelihood_titre(SEXP expectedSEXP, SEXP dataSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type expected(expectedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_titre(expected, data, theta));
    return rcpp_result_gen;
END_RCPP
}
// likelihood_data_individual
double likelihood_data_individual(NumericVector theta, IntegerVector infectionHistory, NumericVector circulationTimes, IntegerVector circulationMapIndices, NumericVector samplingTimes, IntegerVector dataIndices, IntegerVector measuredMapIndices, NumericVector antigenicMapLong, NumericVector antigenicMapShort, int numberStrains, NumericVector data);
RcppExport SEXP _serosolver_likelihood_data_individual(SEXP thetaSEXP, SEXP infectionHistorySEXP, SEXP circulationTimesSEXP, SEXP circulationMapIndicesSEXP, SEXP samplingTimesSEXP, SEXP dataIndicesSEXP, SEXP measuredMapIndicesSEXP, SEXP antigenicMapLongSEXP, SEXP antigenicMapShortSEXP, SEXP numberStrainsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type infectionHistory(infectionHistorySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type circulationTimes(circulationTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type circulationMapIndices(circulationMapIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samplingTimes(samplingTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dataIndices(dataIndicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type measuredMapIndices(measuredMapIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapLong(antigenicMapLongSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapShort(antigenicMapShortSEXP);
    Rcpp::traits::input_parameter< int >::type numberStrains(numberStrainsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_data_individual(theta, infectionHistory, circulationTimes, circulationMapIndices, samplingTimes, dataIndices, measuredMapIndices, antigenicMapLong, antigenicMapShort, numberStrains, data));
    return rcpp_result_gen;
END_RCPP
}
// individual_likelihood
double individual_likelihood(NumericVector theta, NumericVector infectionHistory, NumericVector samplingTimes, IntegerVector indivIndices, NumericVector strainIsolationTimes, IntegerVector indivVirusIndices, NumericVector antigenicMapLong, NumericVector antigenicMapShort, NumericVector titres, int numberStrains);
RcppExport SEXP _serosolver_individual_likelihood(SEXP thetaSEXP, SEXP infectionHistorySEXP, SEXP samplingTimesSEXP, SEXP indivIndicesSEXP, SEXP strainIsolationTimesSEXP, SEXP indivVirusIndicesSEXP, SEXP antigenicMapLongSEXP, SEXP antigenicMapShortSEXP, SEXP titresSEXP, SEXP numberStrainsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type infectionHistory(infectionHistorySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samplingTimes(samplingTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indivIndices(indivIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type strainIsolationTimes(strainIsolationTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indivVirusIndices(indivVirusIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapLong(antigenicMapLongSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapShort(antigenicMapShortSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type titres(titresSEXP);
    Rcpp::traits::input_parameter< int >::type numberStrains(numberStrainsSEXP);
    rcpp_result_gen = Rcpp::wrap(individual_likelihood(theta, infectionHistory, samplingTimes, indivIndices, strainIsolationTimes, indivVirusIndices, antigenicMapLong, antigenicMapShort, titres, numberStrains));
    return rcpp_result_gen;
END_RCPP
}
// group_likelihood_vector
NumericVector group_likelihood_vector(NumericVector theta, NumericMatrix infectionHistories, IntegerVector indicesSamples, IntegerVector indicesData, IntegerVector indicesDataOverall, NumericVector samplingTimes, NumericVector strainIsolationTimes, IntegerVector virusIndices, NumericVector antigenicMapLong, NumericVector antigenicMapShort, NumericVector titres);
RcppExport SEXP _serosolver_group_likelihood_vector(SEXP thetaSEXP, SEXP infectionHistoriesSEXP, SEXP indicesSamplesSEXP, SEXP indicesDataSEXP, SEXP indicesDataOverallSEXP, SEXP samplingTimesSEXP, SEXP strainIsolationTimesSEXP, SEXP virusIndicesSEXP, SEXP antigenicMapLongSEXP, SEXP antigenicMapShortSEXP, SEXP titresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type infectionHistories(infectionHistoriesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesSamples(indicesSamplesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesData(indicesDataSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesDataOverall(indicesDataOverallSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samplingTimes(samplingTimesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type strainIsolationTimes(strainIsolationTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type virusIndices(virusIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapLong(antigenicMapLongSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapShort(antigenicMapShortSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type titres(titresSEXP);
    rcpp_result_gen = Rcpp::wrap(group_likelihood_vector(theta, infectionHistories, indicesSamples, indicesData, indicesDataOverall, samplingTimes, strainIsolationTimes, virusIndices, antigenicMapLong, antigenicMapShort, titres));
    return rcpp_result_gen;
END_RCPP
}
// group_likelihood_total
double group_likelihood_total(NumericVector theta, NumericMatrix infectionHistories, IntegerVector indicesSamples, IntegerVector indicesData, IntegerVector indicesDataOverall, NumericVector samplingTimes, NumericVector strainIsolationTimes, IntegerVector virusIndices, NumericVector antigenicMapLong, NumericVector antigenicMapShort, NumericVector titres);
RcppExport SEXP _serosolver_group_likelihood_total(SEXP thetaSEXP, SEXP infectionHistoriesSEXP, SEXP indicesSamplesSEXP, SEXP indicesDataSEXP, SEXP indicesDataOverallSEXP, SEXP samplingTimesSEXP, SEXP strainIsolationTimesSEXP, SEXP virusIndicesSEXP, SEXP antigenicMapLongSEXP, SEXP antigenicMapShortSEXP, SEXP titresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type infectionHistories(infectionHistoriesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesSamples(indicesSamplesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesData(indicesDataSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesDataOverall(indicesDataOverallSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samplingTimes(samplingTimesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type strainIsolationTimes(strainIsolationTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type virusIndices(virusIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapLong(antigenicMapLongSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapShort(antigenicMapShortSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type titres(titresSEXP);
    rcpp_result_gen = Rcpp::wrap(group_likelihood_total(theta, infectionHistories, indicesSamples, indicesData, indicesDataOverall, samplingTimes, strainIsolationTimes, virusIndices, antigenicMapLong, antigenicMapShort, titres));
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
// c_model_original
NumericVector c_model_original(int n, int nsamp, IntegerVector x, NumericVector theta, NumericVector dd, NumericVector dd2, int t_sample);
RcppExport SEXP _serosolver_c_model_original(SEXP nSEXP, SEXP nsampSEXP, SEXP xSEXP, SEXP thetaSEXP, SEXP ddSEXP, SEXP dd2SEXP, SEXP t_sampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type nsamp(nsampSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dd(ddSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dd2(dd2SEXP);
    Rcpp::traits::input_parameter< int >::type t_sample(t_sampleSEXP);
    rcpp_result_gen = Rcpp::wrap(c_model_original(n, nsamp, x, theta, dd, dd2, t_sample));
    return rcpp_result_gen;
END_RCPP
}
// infection_model_indiv_mus
NumericVector infection_model_indiv_mus(NumericVector theta, NumericVector mus, IntegerVector infectionHistory, NumericVector infectionTimes, IntegerVector boostingVecIndices, IntegerVector infectionMapIndices, double samplingTime, IntegerVector measurementMapIndices, NumericVector antigenicMapLong, NumericVector antigenicMapShort, int numberStrains);
RcppExport SEXP _serosolver_infection_model_indiv_mus(SEXP thetaSEXP, SEXP musSEXP, SEXP infectionHistorySEXP, SEXP infectionTimesSEXP, SEXP boostingVecIndicesSEXP, SEXP infectionMapIndicesSEXP, SEXP samplingTimeSEXP, SEXP measurementMapIndicesSEXP, SEXP antigenicMapLongSEXP, SEXP antigenicMapShortSEXP, SEXP numberStrainsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mus(musSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type infectionHistory(infectionHistorySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type infectionTimes(infectionTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type boostingVecIndices(boostingVecIndicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type infectionMapIndices(infectionMapIndicesSEXP);
    Rcpp::traits::input_parameter< double >::type samplingTime(samplingTimeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type measurementMapIndices(measurementMapIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapLong(antigenicMapLongSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapShort(antigenicMapShortSEXP);
    Rcpp::traits::input_parameter< int >::type numberStrains(numberStrainsSEXP);
    rcpp_result_gen = Rcpp::wrap(infection_model_indiv_mus(theta, mus, infectionHistory, infectionTimes, boostingVecIndices, infectionMapIndices, samplingTime, measurementMapIndices, antigenicMapLong, antigenicMapShort, numberStrains));
    return rcpp_result_gen;
END_RCPP
}
// titre_data_individual_mus
NumericVector titre_data_individual_mus(NumericVector theta, NumericVector mus, IntegerVector infectionHistory, NumericVector circulationTimes, IntegerVector circulationMapIndices, IntegerVector musIndices, NumericVector samplingTimes, IntegerVector dataIndices, IntegerVector measuredMapIndices, NumericVector antigenicMapLong, NumericVector antigenicMapShort, int numberStrains);
RcppExport SEXP _serosolver_titre_data_individual_mus(SEXP thetaSEXP, SEXP musSEXP, SEXP infectionHistorySEXP, SEXP circulationTimesSEXP, SEXP circulationMapIndicesSEXP, SEXP musIndicesSEXP, SEXP samplingTimesSEXP, SEXP dataIndicesSEXP, SEXP measuredMapIndicesSEXP, SEXP antigenicMapLongSEXP, SEXP antigenicMapShortSEXP, SEXP numberStrainsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mus(musSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type infectionHistory(infectionHistorySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type circulationTimes(circulationTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type circulationMapIndices(circulationMapIndicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type musIndices(musIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samplingTimes(samplingTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dataIndices(dataIndicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type measuredMapIndices(measuredMapIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapLong(antigenicMapLongSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapShort(antigenicMapShortSEXP);
    Rcpp::traits::input_parameter< int >::type numberStrains(numberStrainsSEXP);
    rcpp_result_gen = Rcpp::wrap(titre_data_individual_mus(theta, mus, infectionHistory, circulationTimes, circulationMapIndices, musIndices, samplingTimes, dataIndices, measuredMapIndices, antigenicMapLong, antigenicMapShort, numberStrains));
    return rcpp_result_gen;
END_RCPP
}
// titre_data_group_mus
NumericVector titre_data_group_mus(NumericVector theta, NumericVector mus, IntegerMatrix infectionHistories, NumericVector circulationTimes, IntegerVector circulationMapIndices, IntegerVector musIndices, NumericVector samplingTimes, IntegerVector indicesTitreDataSample, IntegerVector indicesTitreDataOverall, IntegerVector indicesSamples, IntegerVector measuredMapIndices, NumericVector antigenicMapLong, NumericVector antigenicMapShort);
RcppExport SEXP _serosolver_titre_data_group_mus(SEXP thetaSEXP, SEXP musSEXP, SEXP infectionHistoriesSEXP, SEXP circulationTimesSEXP, SEXP circulationMapIndicesSEXP, SEXP musIndicesSEXP, SEXP samplingTimesSEXP, SEXP indicesTitreDataSampleSEXP, SEXP indicesTitreDataOverallSEXP, SEXP indicesSamplesSEXP, SEXP measuredMapIndicesSEXP, SEXP antigenicMapLongSEXP, SEXP antigenicMapShortSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mus(musSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type infectionHistories(infectionHistoriesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type circulationTimes(circulationTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type circulationMapIndices(circulationMapIndicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type musIndices(musIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samplingTimes(samplingTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesTitreDataSample(indicesTitreDataSampleSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesTitreDataOverall(indicesTitreDataOverallSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesSamples(indicesSamplesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type measuredMapIndices(measuredMapIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapLong(antigenicMapLongSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapShort(antigenicMapShortSEXP);
    rcpp_result_gen = Rcpp::wrap(titre_data_group_mus(theta, mus, infectionHistories, circulationTimes, circulationMapIndices, musIndices, samplingTimes, indicesTitreDataSample, indicesTitreDataOverall, indicesSamples, measuredMapIndices, antigenicMapLong, antigenicMapShort));
    return rcpp_result_gen;
END_RCPP
}
// infection_history_proposal_gibbs
IntegerMatrix infection_history_proposal_gibbs(NumericVector pars, IntegerMatrix infHist, double indivSampPropn, int n_years_samp, IntegerVector ageMask, IntegerVector n_alive, double swapPropn, int swapDistance, double alpha, double beta, NumericVector circulationTimes, IntegerVector circulationMapIndices, NumericVector samplingTimes, IntegerVector indicesTitreDataSample, IntegerVector indicesTitreDataOverall, IntegerVector indicesSamples, IntegerVector measuredMapIndices, NumericVector antigenicMapLong, NumericVector antigenicMapShort, NumericVector data);
RcppExport SEXP _serosolver_infection_history_proposal_gibbs(SEXP parsSEXP, SEXP infHistSEXP, SEXP indivSampPropnSEXP, SEXP n_years_sampSEXP, SEXP ageMaskSEXP, SEXP n_aliveSEXP, SEXP swapPropnSEXP, SEXP swapDistanceSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP circulationTimesSEXP, SEXP circulationMapIndicesSEXP, SEXP samplingTimesSEXP, SEXP indicesTitreDataSampleSEXP, SEXP indicesTitreDataOverallSEXP, SEXP indicesSamplesSEXP, SEXP measuredMapIndicesSEXP, SEXP antigenicMapLongSEXP, SEXP antigenicMapShortSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type infHist(infHistSEXP);
    Rcpp::traits::input_parameter< double >::type indivSampPropn(indivSampPropnSEXP);
    Rcpp::traits::input_parameter< int >::type n_years_samp(n_years_sampSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ageMask(ageMaskSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_alive(n_aliveSEXP);
    Rcpp::traits::input_parameter< double >::type swapPropn(swapPropnSEXP);
    Rcpp::traits::input_parameter< int >::type swapDistance(swapDistanceSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type circulationTimes(circulationTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type circulationMapIndices(circulationMapIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samplingTimes(samplingTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesTitreDataSample(indicesTitreDataSampleSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesTitreDataOverall(indicesTitreDataOverallSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indicesSamples(indicesSamplesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type measuredMapIndices(measuredMapIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapLong(antigenicMapLongSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapShort(antigenicMapShortSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(infection_history_proposal_gibbs(pars, infHist, indivSampPropn, n_years_samp, ageMask, n_alive, swapPropn, swapDistance, alpha, beta, circulationTimes, circulationMapIndices, samplingTimes, indicesTitreDataSample, indicesTitreDataOverall, indicesSamples, measuredMapIndices, antigenicMapLong, antigenicMapShort, data));
    return rcpp_result_gen;
END_RCPP
}
// inf_hist_prop_cpp
arma::mat inf_hist_prop_cpp(arma::mat infHist, IntegerVector sampledIndivs, IntegerVector ageMask, IntegerVector moveSizes, IntegerVector nInfs, double alpha, double beta, NumericVector randNs);
RcppExport SEXP _serosolver_inf_hist_prop_cpp(SEXP infHistSEXP, SEXP sampledIndivsSEXP, SEXP ageMaskSEXP, SEXP moveSizesSEXP, SEXP nInfsSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP randNsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type infHist(infHistSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sampledIndivs(sampledIndivsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ageMask(ageMaskSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type moveSizes(moveSizesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nInfs(nInfsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type randNs(randNsSEXP);
    rcpp_result_gen = Rcpp::wrap(inf_hist_prop_cpp(infHist, sampledIndivs, ageMask, moveSizes, nInfs, alpha, beta, randNs));
    return rcpp_result_gen;
END_RCPP
}
// subset_test
NumericVector subset_test(NumericVector x, IntegerVector y);
RcppExport SEXP _serosolver_subset_test(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(subset_test(x, y));
    return rcpp_result_gen;
END_RCPP
}
// subset_test1
NumericVector subset_test1(NumericVector x, LogicalVector y);
RcppExport SEXP _serosolver_subset_test1(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(subset_test1(x, y));
    return rcpp_result_gen;
END_RCPP
}
// infection_model_indiv_OLD
NumericVector infection_model_indiv_OLD(NumericVector theta, NumericVector infectionHistory, double samplingTime, NumericVector strainIsolationTimes, IntegerVector virusIndices, NumericVector antigenicMapLong, NumericVector antigenicMapShort, int numberStrains);
RcppExport SEXP _serosolver_infection_model_indiv_OLD(SEXP thetaSEXP, SEXP infectionHistorySEXP, SEXP samplingTimeSEXP, SEXP strainIsolationTimesSEXP, SEXP virusIndicesSEXP, SEXP antigenicMapLongSEXP, SEXP antigenicMapShortSEXP, SEXP numberStrainsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type infectionHistory(infectionHistorySEXP);
    Rcpp::traits::input_parameter< double >::type samplingTime(samplingTimeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type strainIsolationTimes(strainIsolationTimesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type virusIndices(virusIndicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapLong(antigenicMapLongSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type antigenicMapShort(antigenicMapShortSEXP);
    Rcpp::traits::input_parameter< int >::type numberStrains(numberStrainsSEXP);
    rcpp_result_gen = Rcpp::wrap(infection_model_indiv_OLD(theta, infectionHistory, samplingTime, strainIsolationTimes, virusIndices, antigenicMapLong, antigenicMapShort, numberStrains));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_serosolver_infection_model_indiv", (DL_FUNC) &_serosolver_infection_model_indiv, 9},
    {"_serosolver_titre_data_individual", (DL_FUNC) &_serosolver_titre_data_individual, 10},
    {"_serosolver_titre_data_group", (DL_FUNC) &_serosolver_titre_data_group, 11},
    {"_serosolver_likelihood_titre", (DL_FUNC) &_serosolver_likelihood_titre, 3},
    {"_serosolver_likelihood_data_individual", (DL_FUNC) &_serosolver_likelihood_data_individual, 11},
    {"_serosolver_individual_likelihood", (DL_FUNC) &_serosolver_individual_likelihood, 10},
    {"_serosolver_group_likelihood_vector", (DL_FUNC) &_serosolver_group_likelihood_vector, 11},
    {"_serosolver_group_likelihood_total", (DL_FUNC) &_serosolver_group_likelihood_total, 11},
    {"_serosolver_sum_buckets", (DL_FUNC) &_serosolver_sum_buckets, 2},
    {"_serosolver_create_cross_reactivity_vector", (DL_FUNC) &_serosolver_create_cross_reactivity_vector, 2},
    {"_serosolver_c_model_original", (DL_FUNC) &_serosolver_c_model_original, 7},
    {"_serosolver_infection_model_indiv_mus", (DL_FUNC) &_serosolver_infection_model_indiv_mus, 11},
    {"_serosolver_titre_data_individual_mus", (DL_FUNC) &_serosolver_titre_data_individual_mus, 12},
    {"_serosolver_titre_data_group_mus", (DL_FUNC) &_serosolver_titre_data_group_mus, 13},
    {"_serosolver_infection_history_proposal_gibbs", (DL_FUNC) &_serosolver_infection_history_proposal_gibbs, 20},
    {"_serosolver_inf_hist_prop_cpp", (DL_FUNC) &_serosolver_inf_hist_prop_cpp, 8},
    {"_serosolver_subset_test", (DL_FUNC) &_serosolver_subset_test, 2},
    {"_serosolver_subset_test1", (DL_FUNC) &_serosolver_subset_test1, 2},
    {"_serosolver_infection_model_indiv_OLD", (DL_FUNC) &_serosolver_infection_model_indiv_OLD, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_serosolver(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
