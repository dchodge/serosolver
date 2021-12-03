#include "../inst/include/serosolver.h"
#include "boosting_functions_fast.h"

#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b)) // define MAX function for use later
#endif

//' Base boosting fast
//' 
//' A fast implementation of the basic boosting function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
// [[Rcpp::export]]
void titre_data_fast_individual_base(
            NumericVector &predicted_titres,
            const NumericVector &theta,
            const List &infection_info,
            const List &vaccination_info,
            const List &setup_data,
            const List &indexing,
            const List &antigenic_maps,
            const List &other_pars
          ){

  int max_vaccinations;
  bool vac_flag = !vaccination_info["vac_null_ind"];
  NumericVector vaccination_times;
  if (vac_flag) {
    NumericVector vaccination_times = vaccination_info["vac_times"];
    IntegerVector vaccination_strain_indices_tmp = vaccination_info["vac_indices"];
    max_vaccinations = vaccination_times.size();
  }

  NumericVector infection_times = infection_info["inf_times"];
  IntegerVector infection_strain_indices_tmp = infection_info["inf_indices"];

	bool boost_before_infection = false;

  int index_in_samples = indexing["index_in_samples"];
  int end_index_in_samples = indexing["end_index_in_samples"];
  int start_index_in_data = indexing["start_index_in_data"];

  NumericVector sample_times = setup_data["sample_times"];
  IntegerVector measurement_strain_indices = setup_data["measurement_strain_indices"];
  IntegerVector nrows_per_blood_sample = setup_data["nrows_per_blood_sample"];
  int number_strains = setup_data["number_strains"];

  double sampling_time;
  double time;
  double n_inf;
  double n_vac;
  int x_inf;
  int x_vac;

  double wane_amount, wane_amount_vac;
  double seniority;
  double vac_suppress;

  int n_titres;

  int end_index_in_data;
  int tmp_titre_index;
  int inf_map_index;
  int vac_map_index;
  int index;

  int max_infections = infection_times.size();

  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double kappa = theta["kappa"];
  
  NumericVector antigenic_map_long = antigenic_maps["long"];
  NumericVector antigenic_map_short = antigenic_maps["short"];
 // NumericVector antigenic_map_long_vac = antigenic_maps["long_vac"];
  //NumericVector antigenic_map_short_vac = antigenic_maps["short_vac"];

  NumericVector inf_vac_times;
  if (vac_flag) {
    if (max_vaccinations > 0  & max_infections > 0) {
      inf_vac_times = union_(infection_times, vaccination_times);
    } else if(max_vaccinations == 0  & max_infections > 0) {
      inf_vac_times = infection_times;
    } else if (max_vaccinations > 0  & max_infections == 0) {
      inf_vac_times = vaccination_times;
    } 
  }
  else {
    inf_vac_times = infection_times;
  }
  std::sort(inf_vac_times.begin(), inf_vac_times.end());

  int max_inf_vac_times = inf_vac_times.size();

  LogicalVector indicatior_inf(max_inf_vac_times);
  LogicalVector indicatior_vac(max_inf_vac_times);

  for (int i = 0; i < inf_vac_times.size(); i++) {
    if (vac_flag) {
      for (int k = 0; k < vaccination_times.size(); k++) {
        indicatior_vac[i] = inf_vac_times[i] == vaccination_times[k];
        if (indicatior_vac[i])
          break;
      }
    }
    for (int j = 0; j < infection_times.size(); j++) {
      indicatior_inf[i] = inf_vac_times[i] == infection_times[j];
      if (indicatior_inf[i])
        break;
    }
  }
 // Rcpp::Rcout <<  "infection_times: " <<  infection_times << std::endl;
 // Rcpp::Rcout <<  "inf_vac_times: " <<  inf_vac_times << std::endl;
 // Rcpp::Rcout <<  "indicatior_inf: " <<  indicatior_inf << std::endl;

 // Rcpp::Rcout <<  "infection_strain_indices_tmp: " <<  infection_strain_indices_tmp << std::endl;


  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 0.0;
    n_vac = 0.0;
    x_inf = 0.0;
    x_vac = 0.0;

    // Find number of titres in the predicted_titres vector that correspond to this sample
    n_titres = nrows_per_blood_sample[j];
    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;

    // Sum all infections that would contribute towards observed titres at this time
     //  Rcpp::Rcout << "max_inf_vac_times: " << max_inf_vac_times << std::endl;
    for (int x = 0; x < max_inf_vac_times; ++x){
     //    Rcpp::Rcout << "indicatior_inf[x]: " << indicatior_inf[x] << std::endl;
      if (indicatior_inf[x]) {
        ++n_inf;
        // sampling through each predicted infection
        // Only go further if this sample happened after the infection
        if((boost_before_infection && sampling_time > infection_times[x_inf]) ||
          (!boost_before_infection && sampling_time >= infection_times[x_inf])){
     //       Rcpp::Rcout << "Inside real bit" << std::endl;
            time = sampling_time - infection_times[x_inf]; // Time between sample and infection
            wane_amount = MAX(0, 1.0 - (wane*time)); // Basic waning function
            seniority = MAX(0, 1.0 - tau*(n_inf + n_vac - 1.0)); // Antigenic seniority
            
            inf_map_index = infection_strain_indices_tmp[x_inf]; // Index of this infecting strain in antigenic map
            for(int k = 0; k < n_titres; ++k){
              index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
              predicted_titres[tmp_titre_index + k] += seniority *
                ((mu*antigenic_map_long[index]) + (mu_short*antigenic_map_short[index])*wane_amount);
            }
        }
        ++x_inf;
      }

          /*  if (indicatior_vac[x] & (vac_flag)) {
          ++n_vac;
          if((boost_before_infection && sampling_time > vaccination_times[x_vac]) ||
          (!boost_before_infection && sampling_time >= vaccination_times[x_vac])){

            double mu_vac_t = 1.0 / (1.0 + exp(-mu_vac));

            time = sampling_time - vaccination_times[x_vac]; // Time er vaccination
            wane_amount_vac = MAX(0, 1.0 - (wane_vac*time)); // Basic waning function
            seniority = MAX(0, 1.0 - tau_prev_vac*tau*(n_inf + n_vac - 1.0)); // Antigenic seniority
            vac_map_index = vaccination_strain_indices_tmp[x_vac]; // Index of this vaccinating strain in antigenic map

            for(int k = 0; k < n_titres; ++k){
              index = measurement_strain_indices[tmp_titre_index + k]*number_strains + vac_map_index;
              predicted_titres[tmp_titre_index + k] += (seniority) * ((mu*mu_vac_t*antigenic_map_long_vac[index]) + (mu_short_vac*antigenic_map_short_vac[index])*wane_amount_vac);
            }
          }
          ++x_vac;
        }*/
    }
    start_index_in_data = end_index_in_data;
  }
}

 // abfunc titre_data_fast_individual_base_ptr{ &titre_data_fast_individual_base };

//' Alternative waning fast
//' 
//' A fast implementation of the alternative waning function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
// [[Rcpp::export]]
void titre_data_fast_individual_wane2(
            NumericVector &predicted_titres,
            NumericVector &theta,
            const List &infection_info,
            const List &vaccination_info,
            const List &setup_data,
				    const List &indexing,
            const List &antigenic_maps,
            const List &other_pars
				  ){

	bool boost_before_infection = false;

  NumericVector vaccination_times = vaccination_info["vac_times"];
  IntegerVector vaccination_strain_indices_tmp = vaccination_info["vac_indices"];

  NumericVector infection_times = infection_info["inf_times"];
  IntegerVector infection_strain_indices_tmp = infection_info["inf_indices"];

  NumericVector antigenic_map_long = antigenic_maps["long"];
  NumericVector antigenic_map_short = antigenic_maps["short"];

  int index_in_samples = indexing["index_in_samples"];
  int end_index_in_samples = indexing["end_index_in_samples"];
  int start_index_in_data = indexing["start_index_in_data"];

  NumericVector sample_times = setup_data["sample_times"];
  IntegerVector measurement_strain_indices = setup_data["measurement_strain_indices"];
  IntegerVector nrows_per_blood_sample = setup_data["nrows_per_blood_sample"];
  int number_strains = setup_data["number_strains"];

  double sampling_time;
  double time;
  double n_inf;
  double wane_amount;
  double seniority;

  int n_titres;
  int max_infections = infection_times.size();
  int end_index_in_data;
  int tmp_titre_index;
  int inf_map_index;
  int index;

  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double kappa = theta["kappa"];
  double t_change = theta["t_change"];
  double wane_2 = -kappa*wane;
  double wane_2_val; // Interaction term
  // For each sample this individual has
  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 1.0;

    // Find number of titres in the predicted_titres vector that correspond to this sample
    n_titres = nrows_per_blood_sample[j];

    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;

    // Sum all infections that would contribute towards observed titres at this time
    for(int x = 0; x < max_infections; ++x){
      // Only go further if this sample happened after the infection
      //if(sampling_time >= infection_times[x]){
      if((boost_before_infection && sampling_time > infection_times[x]) ||
	 (!boost_before_infection && sampling_time >= infection_times[x])){
	time = sampling_time - infection_times[x]; // Time between sample and infection

	/////////////////////////////////
	// Amanda's alternative waning function
	if(time > t_change){
	  wane_2_val = wane_2*(time - t_change); 
	}else{
	  wane_2_val = 0;
	}
	wane_amount = MAX(0, 1.0 - (wane*time+wane_2_val));
	/////////////////////////////////

	seniority = MAX(0, 1.0 - tau*(n_inf - 1.0)); // Antigenic seniority
	inf_map_index = infection_strain_indices_tmp[x]; // Index of this infecting strain in antigenic map

	// Find contribution to each measured titre from this infection
	for(int k = 0; k < n_titres; ++k){
	  index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
	  predicted_titres[tmp_titre_index + k] += seniority *
	    ((mu*antigenic_map_long[index]) + (mu_short*antigenic_map_short[index])*wane_amount);
	}
	++n_inf;
      }
    }
    start_index_in_data = end_index_in_data;
  }
}


//' Titre dependent boosting fast
//' 
//' A fast implementation of the titre dependent boosting function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
// [[Rcpp::export]]
void titre_data_fast_individual_titredep(
            NumericVector &predicted_titres,
            NumericVector &theta, 
            const List &infection_info,
            const List &vaccination_info,
				  	const List &setup_data,
				    const List &indexing,
            const List &antigenic_maps,
            const List &other_pars
					 ){
	bool boost_before_infection = false;

  NumericVector vaccination_times = vaccination_info["vac_times"];
  IntegerVector vaccination_strain_indices_tmp = vaccination_info["vac_indices"];

  NumericVector infection_times = infection_info["inf_times"];
  IntegerVector infection_strain_indices_tmp = infection_info["inf_indices"];

  NumericVector antigenic_map_long = antigenic_maps["long"];
  NumericVector antigenic_map_short = antigenic_maps["short"];
  double index_in_samples = indexing["index_in_samples"];
  double end_index_in_samples = indexing["end_index_in_samples"];
  double start_index_in_data = indexing["start_index_in_data"];

  NumericVector sample_times = setup_data["sample_times"];
  IntegerVector measurement_strain_indices = setup_data["measurement_strain_indices"];
  IntegerVector nrows_per_blood_sample = setup_data["nrows_per_blood_sample"];
  int number_strains = setup_data["number_strains"];

  double sampling_time;
  double time;
  double n_inf;
  double wane_amount;
  double seniority;
  double infection_time;

  double boost = 0;
  double long_boost=0;
  double monitored_titre=0;
  double short_boost=0;

  int n_titres;
  int max_infections = infection_times.size();
  int end_index_in_data;
  int tmp_titre_index;
  int inf_map_index;
  int inf_map_index_tmp;
  int index;

  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double gradient = theta["gradient"];
  double boost_limit = theta["boost_limit"];
  double titre_suppression = MAX(0, 1.0 - gradient*boost_limit);


  NumericVector monitored_titres(max_infections);

  // For each sample this individual has
  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 0.0;

    // Find number of titres in the predicted_titres vector that correspond to this sample
    n_titres = nrows_per_blood_sample[j];

    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;

    // Sum all infections that would contribute towards observed titres at this time
    for(int x = 0; x < max_infections; ++x){
      // Only go further if this sample happened after the infection
      if((boost_before_infection && sampling_time > infection_times[x]) ||
	          (!boost_before_infection && sampling_time >= infection_times[x])){
        monitored_titre = 0;
        infection_time = infection_times[x];
        time = sampling_time - infection_time; // Time between sample and infection
        inf_map_index = infection_strain_indices_tmp[x]; // Index of this infecting strain in antigenic map

        // Add up contribution of all previous infections to titre that
        // would be observed at this infection time
        for(int ii = x - 1; ii >= 0; --ii){
          inf_map_index_tmp = inf_map_index * number_strains + infection_strain_indices_tmp[ii];
          seniority = MAX(0, 1.0 - tau * ii);
          wane_amount = MAX(0, 1.0 - wane * (infection_time - infection_times[ii]));

          long_boost = seniority * mu * antigenic_map_long[inf_map_index_tmp];
          short_boost = seniority * mu_short * antigenic_map_short[inf_map_index_tmp];
          if(monitored_titres[ii] >= boost_limit){
            long_boost *= titre_suppression;
            short_boost *= titre_suppression;
          } else {
            long_boost *= MAX(0, 1.0 - gradient*monitored_titres[ii]);
            short_boost *= MAX(0, 1.0 - gradient*monitored_titres[ii]);	    
          }
          long_boost = MAX(0, long_boost);
          short_boost = MAX(0, short_boost);
          boost = long_boost + short_boost * wane_amount;
          monitored_titre += boost;
        }
        monitored_titres[x] = monitored_titre;

        wane_amount= MAX(0, 1.0 - (wane*time)); // Basic waning function
        seniority = MAX(0, 1.0 - tau*n_inf); // Antigenic seniority
        
        // Find contribution to each measured titre from this infection
        for(int k = 0; k < n_titres; ++k){
          index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;

          long_boost = seniority * mu * antigenic_map_long[index];
          short_boost = seniority * mu_short * antigenic_map_short[index];
          
          // Titre dependent boosting - at ceiling
          if(monitored_titres[x] >= boost_limit){
            long_boost *= titre_suppression;
            short_boost *= titre_suppression;
          // Titre dependent boosting - below ceiling
          } else {
            long_boost = long_boost * (1 - gradient * monitored_titres[x]); 
            short_boost = short_boost * (1 - gradient * monitored_titres[x]); // Titre dependent boosting - below ceiling
          }
          long_boost = MAX(0, long_boost);
          short_boost = MAX(0, short_boost);
          boost = long_boost + short_boost * wane_amount;
          predicted_titres[tmp_titre_index + k] += boost;
        }
      ++n_inf;
      }
    }
    start_index_in_data = end_index_in_data;
  }
}



//' Base boosting fast
//' 
//' A fast implementation of the basic boosting function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
// [[Rcpp::export]]
void titre_data_fast_individual_strain_dependent(
            NumericVector &predicted_titres,
            NumericVector &theta,
            const List &infection_info,
            const List &vaccination_info,
						const List &setup_data,
				    const List &indexing,
            const List &antigenic_maps,
            const List &other_pars
						 ){

  bool boost_before_infection = false;

  NumericVector infection_times = infection_info["inf_times"];
  IntegerVector infection_strain_indices_tmp = infection_info["inf_indices"];

  NumericVector antigenic_map_long = antigenic_maps["long"];
  NumericVector antigenic_map_short = antigenic_maps["short"];

  double index_in_samples = indexing["index_in_samples"];
  double end_index_in_samples = indexing["end_index_in_samples"];
  double start_index_in_data = indexing["start_index_in_data"];

  NumericVector sample_times = setup_data["sample_times"];
  IntegerVector measurement_strain_indices = setup_data["measurement_strain_indices"];
  IntegerVector nrows_per_blood_sample = setup_data["nrows_per_blood_sample"];
  int number_strains = setup_data["number_strains"];

	NumericVector mus = other_pars["mu"];
  IntegerVector boosting_vec_indices = other_pars["boosting_vec_indices"];

  double sampling_time;
  double time;
  double n_inf;
  double wane_amount;
  double seniority;

  double mu;

  int n_titres;
  int max_infections = infection_times.size();
  int end_index_in_data;
  int tmp_titre_index;
  int inf_map_index;
  int index;

  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];

  // For each sample this individual has
  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 1.0;

    // Find number of titres in the predicted_titres vector that correspond to this sample
    n_titres = nrows_per_blood_sample[j];

    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;

    // Sum all infections that would contribute towards observed titres at this time
    for(int x = 0; x < max_infections; ++x){
      // Only go further if this sample happened after the infection
        if((boost_before_infection && sampling_time > infection_times[x]) ||
	   (!boost_before_infection && sampling_time >= infection_times[x])){
	  //if(sampling_time >= infection_times[x]){
	time = sampling_time - infection_times[x]; // Time between sample and infection
	wane_amount= MAX(0, 1.0 - (wane*time)); // Basic waning function
	seniority = MAX(0, 1.0 - tau*(n_inf - 1.0)); // Antigenic seniority
	inf_map_index = infection_strain_indices_tmp[x]; // Index of this infecting strain in antigenic map
	mu = mus[boosting_vec_indices[inf_map_index]];
	// Find contribution to each measured titre from this infection
	for(int k = 0; k < n_titres; ++k){
	  index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
	  predicted_titres[tmp_titre_index + k] += seniority *
	    ((mu*antigenic_map_long[index]) + (mu_short*antigenic_map_short[index])*wane_amount);
	}
	++n_inf;
      }
    }
    start_index_in_data = end_index_in_data;
  }
}



void titre_data_fast_custom(
            NumericVector &predicted_titres,
            NumericVector &theta,
            const List &infection_info,
            const List &vaccination_info,
            const List &setup_data,
            const List &indexing,
            const List &antigenic_maps,
            const List &other_pars
          ) 
{
  ////////////////////////////////////
  // Get all the values from Lists //
  ////////////////////////////////////
  // vacciantion_info  - get the vaccinations times and strain indicies for individual i
  // - only run if provided
  int max_vaccinations;
  bool vac_flag = !vaccination_info["vac_null_ind"];
  NumericVector vaccination_times;
  if (vac_flag) {
    NumericVector vaccination_times = vaccination_info["vac_times"];
    IntegerVector vaccination_strain_indices_tmp = vaccination_info["vac_indices"];
    max_vaccinations = vaccination_times.size();
  }

  // info_info - get the infection times and strain indicies for individual i
  NumericVector infection_times = infection_info["inf_times"];
  IntegerVector infection_strain_indices_tmp = infection_info["inf_indices"];

  // Get some index values assocaited with individual i (Derived from titre data)
  double index_in_samples = indexing["index_in_samples"];
  double end_index_in_samples = indexing["end_index_in_samples"];
  double start_index_in_data = indexing["start_index_in_data"];

  // Get some infomation from titre_data
  NumericVector sample_times = setup_data["sample_times"];
  IntegerVector measurement_strain_indices = setup_data["measurement_strain_indices"];
  IntegerVector nrows_per_blood_sample = setup_data["nrows_per_blood_sample"];
  int number_strains = setup_data["number_strains"];

  // info_info - get the antigen maps defined in function
  NumericVector antigenic_map_long = antigenic_maps["long"];
  NumericVector antigenic_map_short = antigenic_maps["short"];

  // Add function here 

}