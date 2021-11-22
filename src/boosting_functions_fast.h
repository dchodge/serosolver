#include <Rcpp.h>
using namespace Rcpp;

#ifndef TITRE_DATA_FAST_INDIVIDUAL_BASE_H
#define TITRE_DATA_FAST_INDIVIDUAL_BASE_H
void titre_data_fast_individual_base(
					NumericVector &predicted_titres,
					const NumericVector &theta,
                    const List & infection_info,
                    const List &vaccination_info,
					const List &setup_dat,
					const List &indexing,
				    const List &antigenic_maps,
					const List &other_pars
				    );
#endif


#ifndef TITRE_DATA_FAST_INDIVIDUAL_WANE2_H
#define TITRE_DATA_FAST_INDIVIDUAL_WANE2_H
void titre_data_fast_individual_wane2(
					  NumericVector &predicted_titres,
					  NumericVector &theta,
                      const List & infection_info,
					   const List &vaccination_info,
					  const List &setup_dat,
					  const List &indexing,
				    	const List &antigenic_maps,
						const List &other_pars
				      );
#endif

#ifndef TITRE_DATA_FAST_INDIVIDUAL_TITREDEP_H
#define TITRE_DATA_FAST_INDIVIDUAL_TITREDEP_H
void titre_data_fast_individual_titredep(
					NumericVector &predicted_titres,
					NumericVector &theta,
                    const List & infection_info,					   
					const List &vaccination_info,
					const List &setup_dat,
					const List &indexing,
				    const List &antigenic_maps,
					const List &other_pars
				);
#endif

#ifndef TITRE_DATA_FAST_INDIVIDUAL_STRAIN_DEPENDENT_H
#define TITRE_DATA_FAST_INDIVIDUAL_STRAIN_DEPENDENT_H
void titre_data_fast_individual_strain_dependent(NumericVector &predicted_titres,
						 NumericVector &theta,
                         const List & infection_info,
						 const List &vaccination_info,
						 const List &setup_dat,
					     const List &indexing,
				    const List &antigenic_maps,
					const List &other_pars
						 );
#endif
