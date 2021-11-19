#include <cmath>
#include "wane_function.h"
#include "boosting_functions_fast.h"
#include "helpers.h"

//' Overall model function, fast implementation
//'
//' @param theta NumericVector, the named vector of model parameters
//' @param infection_history_mat IntegerMatrix, the matrix of 1s and 0s showing presence/absence of infection for each possible time for each individual. 
//' @param circulation_times NumericVector, the actual times of circulation that the infection history vector corresponds to
//' @param circulation_times_indices IntegerVector, which entry in the melted antigenic map that these infection times correspond to
//' @param sample_times NumericVector, the times that each blood sample was taken
//' @param rows_per_indiv_in_samples IntegerVector, one entry for each individual. Each entry dictates how many indices through sample_times to iterate per individual (ie. how many sample times does each individual have?)
//' @param cum_nrows_per_individual_in_data IntegerVector, How many cumulative rows in the titre data correspond to each individual?
//' @param nrows_per_blood_sample IntegerVector, one entry per sample taken. Dictates how many entries to iterate through cum_nrows_per_individual_in_data for each sampling time considered
//' @param measurement_strain_indices IntegerVector, the indices of all measured strains in the melted antigenic map, with one entry per measured titre
//' @param antigenic_map_long NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1 see \code{\link{create_cross_reactivity_vector}}
//' @param antigenic_map_short NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2, see \code{\link{create_cross_reactivity_vector}}
//' @param antigenic_distances NumericVector, the collapsed cross reactivity map giving euclidean antigenic distances, see \code{\link{create_cross_reactivity_vector}}
//' @param mus NumericVector, if length is greater than one, assumes that strain-specific boosting is used rather than a single boosting parameter
//' @param boosting_vec_indices IntegerVector, same length as circulation_times, giving the index in the vector \code{mus} that each entry should use as its boosting parameter.
//' @param boost_before_infection bool to indicate if calculated titre for that time should be before the infection has occurred, used to calculate titre-mediated immunity
//' @return NumericVector of predicted titres for each entry in measurement_strain_indices
//' @export
//' @family titre_model
// [[Rcpp::export(rng = false)]]
NumericVector titre_data_fast(NumericVector theta,
			      const IntegerMatrix infection_history_mat,
            const List vaccination_hist_info,
				    const Function abkinetics_model,
			      const NumericVector circulation_times,
			      const IntegerVector circulation_times_indices,
			      const NumericVector sample_times,
			      const IntegerVector rows_per_indiv_in_samples, // How many rows in titre data correspond to each individual, sample and repeat?
			      const IntegerVector cum_nrows_per_individual_in_data, // How many rows in the titre data correspond to each individual?
			      const IntegerVector nrows_per_blood_sample, // Split the sample times and runs for each individual
			      const IntegerVector measurement_strain_indices, // For each titre measurement, corresponding entry in antigenic map
			      const List antigenic_maps, // add in antigenic_map_long_vac ?
			      const NumericVector antigenic_distances,	// Currently not doing anything, but has uses for model extensions		      
			      const NumericVector mus,
			      const IntegerVector boosting_vec_indices,
			      bool boost_before_infection = false
			      ){

  // Dimensions of structures
  int n = infection_history_mat.nrow();
  int number_strains = infection_history_mat.ncol();
  int total_titres = measurement_strain_indices.size();

  // To track how far through the larger vectors we move for each individual
  int index_in_samples;
  int end_index_in_samples;
  int start_index_in_data;
  
  // Only use the infections that actually happened
  IntegerVector infection_history(number_strains);
  LogicalVector indices;
  NumericVector infection_times;
  IntegerVector infection_strain_indices_tmp;

  // Only use the vacciation that actually happened
  NumericMatrix vac_history_matrix;
  NumericVector vac_history_strains;
  NumericVector vac_history_strains_indices;
  bool vac_null_ind = false;
  
  if (vaccination_hist_info["vac_history_matrix"] == R_NilValue) {
    vac_null_ind = true;
  } else {
    NumericMatrix vac_history_matrix = vaccination_hist_info["vac_history_matrix"];
    NumericVector vac_history_strains = vaccination_hist_info["vac_history_strains"];
    NumericVector vac_history_strains_indices = vaccination_hist_info["vac_history_strains_indices"];
  }

  IntegerVector vaccination_history;
  LogicalVector indices_vac;
  NumericVector vaccination_times;
  IntegerVector vaccination_strain_indices_tmp;

  bool base_function = true;
  bool alternative_wane_func = false;
  bool titre_dependent_boosting = false;
  bool strain_dep_boost = false;
  if (mus.size() > 1) {
    strain_dep_boost = false;    
  }

  double min_titre = 0;
  NumericVector predicted_titres(total_titres, min_titre);

  // Create the Lists for cycling
  List setup_dat = List::create(
    _("measurement_strain_indices") = measurement_strain_indices, 
    _["sample_times"] = sample_times,
    _["nrows_per_blood_sample"] = nrows_per_blood_sample,
    _["number_strains"] = number_strains
  );

  List indexing = List::create(
      _("index_in_samples") = R_NilValue, 
      _["end_index_in_samples"] = R_NilValue,
      _["start_index_in_data"] = R_NilValue
  );

  List infection_info = List::create( 
      _["inf_times"] = R_NilValue,
      _["inf_indices"] = R_NilValue
  );

  List vaccination_info_i = List::create( 
      _["vac_times"] = R_NilValue,
      _["vac_indices"] = R_NilValue,
      _["vac_null_ind"] = vac_null_ind
  );

  for (int i = 1; i <= n; ++i) {
    // Find infection times for individual i
    infection_history = infection_history_mat(i-1, _);
    //Rcpp::Rcout << "infection_history: " << infection_history << ". Number: "<< i << std::endl;
    indices = infection_history > 0;
    infection_times = circulation_times[indices];

    // Find vaccination history for individual i
    if (!vac_null_ind) {
      vaccination_history = vac_history_matrix(i-1, _);
      indices_vac = vaccination_history > 0;
      vaccination_times = vac_history_strains[indices_vac];
    }

    // Only run if either an infection or vaccination exists in person i's history
    if (infection_times.size() > 0 || vaccination_times.size() > 0) {
      // Infection info fill for individual i
      infection_strain_indices_tmp = circulation_times_indices[indices];
      infection_info["inf_times"] = infection_times;
      infection_info["inf_indices"] = infection_strain_indices_tmp;

      // Vaccination fill infro for individual i 
      if (!vac_null_ind) {
        vaccination_strain_indices_tmp = vac_history_strains_indices[indices_vac];
        vaccination_info_i["vac_times"] = vaccination_times;
        vaccination_info_i["vac_indices"] = vaccination_strain_indices_tmp;
      }

      // Indexing info fill
      indexing["index_in_samples"] = rows_per_indiv_in_samples[i-1];
      indexing["end_index_in_samples"] = rows_per_indiv_in_samples[i] - 1;
      indexing["start_index_in_data"] = cum_nrows_per_individual_in_data[i-1];

      // ====================================================== //
      // =============== CHOOSE MODEL TO SOLVE =============== //
      // ====================================================== //
      // Go to sub function - this is where we have options for different models
      // Note, these are in "boosting_functions.cpp"
      if (base_function) {

      //  titre_data_fast_individual_base(
            abkinetics_model(
              predicted_titres, 
              theta,
              infection_info,
              vaccination_info_i,
              setup_dat, 
              indexing,
              antigenic_maps);
      } else if (titre_dependent_boosting) {
	          titre_data_fast_individual_titredep(
              predicted_titres, 
              theta, 
					    infection_info,
              vaccination_info_i,
              setup_dat, 
              indexing,
					    antigenic_maps);	
      } else if (strain_dep_boost) {
            titre_data_fast_individual_strain_dependent(
              predicted_titres, 
              mus, boosting_vec_indices, 
              theta,
              infection_info,
              vaccination_info_i,
              setup_dat, 
              indexing,
              antigenic_maps);
      } else if(alternative_wane_func) {
            titre_data_fast_individual_wane2(
              predicted_titres, 
              theta,
              infection_info,
              vaccination_info_i,
              setup_dat, 
              indexing,
              antigenic_maps);
      } else {
	         //titre_data_fast_individual_base(
            abkinetics_model(
              predicted_titres, 
              theta,
              infection_info,
              vaccination_info_i,
              setup_dat, 
              indexing,
              antigenic_maps);
      }
    }
  }
  return(predicted_titres);
}
