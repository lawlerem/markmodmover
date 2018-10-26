#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
	//
	// Throughout, the last dimension of arrays should be indexing the observations for ease of coding
	//

	DATA_VECTOR(theta);
	DATA_VECTOR(st_dist);
	DATA_VECTOR(st_dist_starts);
	DATA_VECTOR(grouping);

	DATA_INTEGER(step_distribution); // 1 = Gamma; 2 = Log-Normal
	DATA_INTEGER(step_zero_inflation); // 0 = No; 1 = Yes

	DATA_INTEGER(angle_zero_inflation); // 0 = No; 1 = Yes

	// Parameters for transition probability matrix
	PARAMETER_ARRAY(tpm_working_pars_matrix);	// Row is the current state. Each next state has one parameter, although we won't use the parameter for staying in the same state (on the diagonals) for a size (n.state,n.state)

	// Parameters for theta distribution (wrapped Cauchy distribution)
	PARAMETER_ARRAY(theta_working_pars);	// Row indexes state, column 1 is mu and column 2 is logitRho

	// Parameters for st_dist distribution (Normal distribution)
	PARAMETER_ARRAY(dist_working_pars);	// Row indexes state, column 1 is mu and column 2 is logSigma
	PARAMETER_VECTOR(acf_working_pars);

	// Parameters for zero inflation
	PARAMETER_VECTOR(logit_step_zero_probs);
	PARAMETER_VECTOR(logit_angle_zero_probs);
	PARAMETER_ARRAY(angle_zero_working_pars);	// Row indexes state, column 1 is mu and column 2 is logitRho

	PARAMETER(dummy);	// Allows us to hold all the other parameters constant when running a second time to permute states


	Type pi = atan(1)*4;
	int n_state = theta_working_pars.rows();
	int n_obs = theta.size();
	int grp;





	//
	//
	//	Transition Probability Matrix, Initial Distribution, State Simulation
	//
	//

	// Transition Probability Matrix from beta parameters
	array<Type> tpm(n_state, n_state);
	Type row_sum = 0.0;
	for(int i=0; i<n_state; i++) {    // Need to fill in by row-first to do row sums easily
		for(int j=0; j<n_state; j++) {
			if( i==j ) {
				tpm(i,j) = 1.0;	// Set diagonal entries equal to one (not free to vary)
			} else {
				tpm(i,j) = 1.0 / (
						1.0 + exp(tpm_working_pars_matrix(i,j))
					); // Off-diagonals use parameters (free to vary
			}
			row_sum += tpm(i,j);	// Adds entry to row sum
		}

		// Done filling in columns for a particular row
		for(int j=0; j<n_state; j++) {
			tpm(i,j) = tpm(i,j) / row_sum; // Enforces row sum equal to 1
		}

		// Done looping through columns, move on to next row (with new row_sum = 0)
		row_sum = 0;
	}
	//
	REPORT(tpm_working_pars_matrix);
	REPORT(tpm);
	//



	// System of equations for finding stationary distribution:	'I - Gamma + 1'
	matrix<Type> system(n_state, n_state);
	for(int i=0; i<n_state; i++) {
		for(int j=0; j<n_state; j++) {
			if( i ==j ) {
				system(i,j) = 1 - tpm(j,i) + 1;
			} else {
				system(i,j) = 0 - tpm(j,i) + 1;
			}
		}
	}

	// vector of ones for finding stationary distribution
	vector<Type> ones(n_state);
	for(int i=0; i<n_state; i++) {
		ones(i) = 1.0;
	}

	matrix<Type> sys_inverse = system.inverse();
	vector<Type> start_probs = sys_inverse*ones;

	matrix<Type> start_probs_row(1,n_state);
	for(int i=0; i<n_state; i++) {
		start_probs_row(0,i) = start_probs(i);
	}
	//
	REPORT(start_probs);
	//



	// Simulation for state process
	vector<int> sim_state_path(n_obs);
	SIMULATE {
		vector<Type> rand(n_obs);

		vector<Type> zero(n_obs);
		vector<Type> one(n_obs);
		zero = 0.0;
		one = 1.0;

		rand = runif(zero,one);
		Type check = 0.0;
		int sim_state = 0;

		while( rand(0) > check ) {
			check += start_probs(sim_state);
			sim_state++;
		}
		sim_state_path(0) = sim_state;  // sim_state will be between (inclusive) 1 and n_state

		for(int k=1; k<n_obs; k++) {
			int last = sim_state_path(k-1) - 1;	// Subtract 1 since we need to start indices from 0

			check = 0.0;
			sim_state = 0;

			while( rand(k) > check ) {
				check += tpm(last,sim_state);
				sim_state++;
			}
			sim_state_path(k) = sim_state;

		}
		//
		REPORT(sim_state_path);
		//
	}





	//
	//
	//	Natural parameters, likelihood, cdf, and simulation for theta (Wrapped Cauchy)
	//
	//
	array<Type> theta_pars = theta_working_pars;
	theta_pars.col(1) = 1/(
			1+exp(theta_working_pars.col(1))
		);	// concentration between 0 and 1

	vector<Type> angle_zero_probs(n_state);
	array<Type> angle_zero_pars = angle_zero_working_pars;

	for(int i=0; i<n_state; i++) {
		angle_zero_probs(i) = exp(logit_angle_zero_probs(i))/(1 + exp(logit_angle_zero_probs(i)));		//  Between 0 and 1
	}

	angle_zero_pars.col(0) = 0.1 - 0.1*exp(angle_zero_working_pars.col(0))/(1 + exp(angle_zero_working_pars.col(0)));	// Between -0.01 and 0.01
	angle_zero_pars.col(1) = 0.95 + 0.05*exp(angle_zero_working_pars.col(1))/(1 + exp(angle_zero_working_pars.col(1)));	// Between 0.95 and 1.0
//	angle_zero_pars.col(1) = exp(angle_zero_working_pars.col(1))/(1 + exp(angle_zero_working_pars.col(1)));	// Between 0 and 1.0
	//
	REPORT(theta_working_pars);
	REPORT(theta_pars);
	ADREPORT(theta_pars);

	REPORT(logit_angle_zero_probs);
	ADREPORT(angle_zero_probs);
	REPORT(angle_zero_probs);

	REPORT(angle_zero_working_pars);
	ADREPORT(angle_zero_pars);
	REPORT(angle_zero_pars);
	//

	// Probability Density Matrices for theta
	array<Type> theta_array(n_state, n_state, n_obs);
	for(int k=0; k<n_obs; k++) {	// k indexes time
		for(int i=0; i<n_state; i++) {
			for(int j=0; j<n_state; j++) {
				if( i==j ) {
					theta_array(i,j,k) = 1.0/(2*pi) * (
							 (1-theta_pars(i,1)*theta_pars(i,1))
							 /
							 (1+theta_pars(i,1)*theta_pars(i,1) - 2*theta_pars(i,1)*cos(theta(k) - theta_pars(i,0)))
						);
				} else {
					theta_array(i,j,k) = 0.0;
				}
			}
		}
	}

	if( angle_zero_inflation == 1.0 ) {
		for(int k=0; k<n_obs; k++) {
			for(int i=0; i<n_state; i++) {
				Type zero_dist = 1.0/(2*pi) * (
						(1-pow(angle_zero_pars(i,1),2))
						/
						(1+pow(angle_zero_pars(i,1),2) - 2*angle_zero_pars(i,1)*cos(theta(k) - angle_zero_pars(i,0) ))
					);
				theta_array(i,i,k) = angle_zero_probs(i)*zero_dist + (1-angle_zero_probs(i))*theta_array(i,i,k);
			}
		}
	}
	//
	REPORT(theta_array);
	//

	// CDF Matrices for theta. Need to make this mod 1 in R.
	array<Type> theta_cdf_array(n_state, n_state, n_obs);
	for(int k=0; k<n_obs; k++) {
		for(int i=0; i<n_state; i++) {
			for(int j=0; j<n_state; j++) {
				if( i==j ) {
                    if( theta_pars(i,0) == -pi ){
                        theta_cdf_array(i,j,k) = (-1/pi) * (
                        		atan(
                        			(theta_pars(i,1) + 1)
                        			/
                        			(theta_pars(i,1) - 1)
                        			*
                        			tan(
                        				0.5*(theta(k) - theta_pars(i,0))
                        			)
                        		) - atan(
                        			(theta_pars(i,1) + 1)
                        			/
                        			(theta_pars(i,1) - 1)
                        			*
                        			tan(
                        				0.5*(-pi+0.00000001 - theta_pars(i,0))
                        			)
                        		)
                        	);
                    } else {
                    theta_cdf_array(i,j,k) = (-1/pi) * (
                    		atan(
                    			(theta_pars(i,1) + 1)
                    			/
                    			(theta_pars(i,1) - 1)
                    			*
                    			tan(
                    				0.5*(theta(k) - theta_pars(i,0))
                    			)
                    		)- atan(
                    			(theta_pars(i,1) + 1)
                    			/
                    			(theta_pars(i,1) - 1)
                    			*
                    			tan(
                    				0.5*(-pi - theta_pars(i,0))
                    			)
                    		)
                    	);
                    }
                } else {
                    theta_cdf_array(i,j,k) = 0.0;
                }
			}
		}
	}
	//
	REPORT(theta_cdf_array);
	//




	// Simulation for theta
	vector<Type> sim_theta_series(n_obs); // Make this between -180 and 180 in R
	SIMULATE {
		for(int k=0; k<n_obs; k++) {
			int state = sim_state_path(k) - 1;

			vector<Type> zero(n_obs);
			vector<Type> one(n_obs);
			zero = 0.0;
			one = 1.0;
			vector<Type> rand = runif(zero,one);

			Type peak = theta_pars(state,0);
			Type conc = theta_pars(state,1);
			Type q0 = atan( (conc+1)/(conc-1) * tan( (-pi - peak) / 2));

			sim_theta_series(k) = 2*atan(
					(conc-1)
					/
					(conc+1)
					*
					tan(
						-pi*rand(k)+q0
					)
				) + peak;
			sim_theta_series(k) = 180/pi * sim_theta_series(k);
		}
		//
		REPORT(sim_theta_series);
		//
	}













	//
	//
 	//	Natural parameters, likelihood, cdf, and simulation for st_dist (1 == Autocorrelated Gamma; 2 == Autocorrelated Log-Normal)
	//
	//
	array<Type> dist_pars = dist_working_pars;
	if( step_distribution == 1 ) {
		dist_pars.col(0) = exp(dist_working_pars.col(0));
		dist_pars.col(1) = exp(dist_working_pars.col(1))+0.01;
	} else if( step_distribution == 2 ) {
		dist_pars.col(1) = exp(dist_working_pars.col(1))+0.01;
	}

	vector<Type> acf_pars = acf_working_pars;
	if( step_distribution == 1 ) {
		for(int i=0; i<n_state; i++) {
			acf_pars(i) = exp(acf_working_pars(i));
		}
	} else if( step_distribution == 2 ) {
		for(int i=0; i<n_state; i++) {
			acf_pars(i) = acf_working_pars(i);
		}
	}
	vector<Type> step_zero_probs = exp(logit_step_zero_probs)/(1+exp(logit_step_zero_probs));
	//
	REPORT(dist_working_pars);
	REPORT(dist_pars);
	ADREPORT(dist_pars);
	REPORT(acf_working_pars);
	ADREPORT(acf_pars);
	REPORT(acf_pars);
	REPORT(logit_step_zero_probs);
	REPORT(step_zero_probs);
	//




	// Probability Density Matrices for st_dist (1 == Autocorrelated Gamma; 2 == Autocorrelated Log-Normal)



	// Initial Contribution
	array<Type> dist_array(n_state, n_state, n_obs);

	grp = 0;



	// Series Contribution
	if( step_distribution == 1 ) {
		for(int k=0; k<n_obs; k++) {	// k indexes time
			for(int i=0; i<n_state; i++) {
				for(int j=0; j<n_state; j++) {
					if( i==j ) {
						if( st_dist(k) == 0.0 ) {
							dist_array(i,j,k) = 0.0;
						} else if( k == grouping(grp) ) {		// If k is in grouping, then st_dist(k) is part of a new contiguous block
							Type mean = acf_pars(i)*st_dist_starts(grp) + dist_pars(i,0);
						//Type mean = acf_pars(i) + dist_pars(i,0);
							Type sd = dist_pars(i,1);

							dist_array(i,j,k) =
								dgamma(st_dist(k),
										pow(mean/sd,2),
										pow(sd,2)/mean);
						} else {
							Type mean = acf_pars(i)*st_dist(k-1) + dist_pars(i,0);
							Type sd = dist_pars(i,1);

							dist_array(i,j,k) =
								dgamma(st_dist(k),
										pow(mean/sd,2),
										pow(sd,2)/mean);
						}
					} else {
						dist_array(i,j,k) = 0.0;
					}
				}
			}
			 if( k == grouping(grp) ) grp++;
		}
	} else if( step_distribution == 2 ) {
		for(int k=0; k<n_obs; k++) {	// k indexes time
			for(int i=0; i<n_state; i++) {
				for(int j=0; j<n_state; j++) {
					if( i==j ) {
						if( st_dist(k) == 0.0 ) {
							dist_array(i,j,k) = 0.0;
						} else if( k == grouping(grp) ) {
							dist_array(i,j,k) =
								dnorm(log(st_dist(k)),
									acf_pars(i)*log(st_dist_starts(grp)) + dist_pars(i,0),
									dist_pars(i,1),
									false);
						} else {
							dist_array(i,j,k) =
								dnorm(log(st_dist(k)),
									acf_pars(i)*log(st_dist(k-1)) + dist_pars(i,0),
									dist_pars(i,1),
									false);
						}
					} else {
						dist_array(i,j,k) = 0.0;
					}
				}
			}
			if( k == grouping(grp) ) grp++;
		}
	}
	grp = 0;




	// Zero-inflation for step length
	if( step_zero_inflation == 1 ){
		for(int k=0; k<n_obs; k++) {
			for(int i=0; i<n_state; i++) {
				Type is_zero;
				if( st_dist(k) == 0.0 ){
					is_zero = 1.0;
				} else {
					is_zero = 0.0;
				}
				dist_array(i,i,k) = step_zero_probs(i)*is_zero + (1.0-step_zero_probs(i))*dist_array(i,i,k);
			}
		}
	}
	//
	REPORT(dist_array);
	//




	// CDF Matrix for st_dist
	array<Type> dist_cdf_array(n_state, n_state, n_obs);



	// Series contribution
	if( step_distribution == 1 ){
		for(int k=0; k<n_obs; k++) {
			for(int i=0; i<n_state; i++) {
				for(int j=0; j<n_state; j++) {
					if( i==j ) {
						if( st_dist(k) == 0.0 ) {
							dist_cdf_array(i,j,k) = 0.0;
						} else if( k == grouping(grp) ) {
							Type mean = acf_pars(i)*st_dist_starts(grp) + dist_pars(i,0);
							Type sd = dist_pars(i,1);

							dist_cdf_array(i,j,k) =
								pgamma(st_dist(k),
									pow(mean/sd,2),
									pow(sd,2)/mean);
						} else {
							Type mean = acf_pars(i)*st_dist(k-1) + dist_pars(i,0);
							Type sd = dist_pars(i,1);

							dist_cdf_array(i,j,k) =
								pgamma(st_dist(k),
									pow(mean/sd,2),
									pow(sd,2)/mean);
						}
					} else {
						dist_cdf_array(i,j,k) = 0.0;
					}
				}
			}
			if( k == grouping(grp) ) grp++;
		}
	} else if( step_distribution == 2 ){
		for(int k=0; k<n_obs; k++) {
			for(int i=0; i<n_state; i++) {
				for(int j=0; j<n_state; j++) {
					if( i==j ) {
						if( st_dist(k) == 0.0 ) {
							dist_cdf_array(i,j,k) = 0.0;
						} else if( k == grouping(grp) ) {
							dist_cdf_array(i,j,k) =
								pnorm(log(st_dist(k)),
									acf_pars(i)*log(st_dist_starts(grp)) + dist_pars(i,0),
									dist_pars(i,1));
						} else {
							dist_cdf_array(i,j,k) =
								pnorm(log(st_dist(k)),
									acf_pars(i)*log(st_dist(k-1)) + dist_pars(i,0),
									dist_pars(i,1));
						}
					} else {
						dist_cdf_array(i,j,k) = 0.0;
					}
				}
			}
			if( k == grouping(grp) ) grp++;
		}
	}
	grp = 0;


	// Zero inflation
	if( step_zero_inflation == 1 ){
		for(int k=0; k<n_obs; k++) {
			for(int i=0; i<n_state; i++) {
				dist_cdf_array(i,i,k) = step_zero_probs(i) + (1-step_zero_probs(i))*dist_cdf_array(i,i,k);
			}
		}
	}
	//
	REPORT(dist_cdf_array);
	//



	// Simulation for st_dist
	SIMULATE {
		vector<Type> sim_dist_series(n_obs+1);
		vector<Type> rand(n_obs);

		if( step_zero_inflation == 1 ){
			vector<Type> zero(n_obs);
			vector<Type> one(n_obs);
			zero = 0.0;
			one = 1.0;
			rand = runif(zero,one);
		} else {
			rand = 1.0;
		}

		if( step_distribution == 1 ){
			sim_dist_series(0) = st_dist_starts(0);

			for(int k=0; k<n_obs; k++) {
				if( step_zero_inflation == 1 && rand(k) < step_zero_probs(sim_state_path(k)) ) {
					sim_dist_series(k+1) = 0.0;
				} else {
					Type mean = acf_pars(sim_state_path(k)-1) * sim_dist_series(k-1+1) + dist_pars(sim_state_path(k)-1,0);
					Type sd = dist_pars(sim_state_path(k)-1,1);

					sim_dist_series(k+1) = rgamma(pow(mean/sd,2),
								      pow(sd,2)/mean);
				}
			}
		} else if( step_distribution == 2 ){
			vector<Type> mu(n_obs);	// need a vector to get n_obs random normal variables
			vector<Type> sd(n_obs);
			for(int k=0; k<n_obs; k++) {
				int state = sim_state_path(k) - 1; // need minus 1 since index starts at 0

				mu(k) = dist_pars(state,0);
				sd(k) = dist_pars(state,1);
			}

			vector<Type> sim_error_series = rnorm(mu,sd); // Call rnorm here to get better run time

			sim_dist_series(0) = log(st_dist_starts(0));

			for(int k=0; k<n_obs; k++) {
				int state = sim_state_path(k) - 1;

				sim_dist_series(k+1) = acf_pars(state)*sim_dist_series(k-1+1) + sim_error_series(k);
			}

			for(int k=0; k<n_obs+1; k++) {
				sim_dist_series(k) = exp(sim_dist_series(k)); // transform back to ratio scale
			}
		}
		//
		REPORT(sim_dist_series);
		//
	}






	//
	//
	// Viterbi algorithm for global decoding of states
	//
	//
	array<Type> viterbi(n_state, n_obs);	// Columns index time, rows index state
	array<Type> forward_max(n_state);	// Vector of most likely states
	array<Type> l_array = theta_array * dist_array;
	vector<int> viterbi_path(n_obs);


	vector<int> int_grouping(grouping.size());
	for(int k=0; k<=n_obs; k++){
		if( k == grouping(grp) ){
			int_grouping(grp) = k;
			grp++;
		}
	}
	grp = 0;

	for(grp=0; grp<(grouping.size()-1); grp++) {
		int min = int_grouping(grp);
		int max = int_grouping(grp+1);


		//starting state likelihoods
		for(int i=0; i<n_state; i++) {
			viterbi(i,min) = log(start_probs(i)*l_array(i,i,min));
		}

		for(int k=(min+1); k<max; k++) {
			for(int i=0; i<n_state; i++) {
				for(int j=0; j<n_state; j++) {
					forward_max(j) = viterbi(j,k-1) + log(tpm(i,j)); // transition probabilities
				}
				viterbi(i,k) = forward_max.maxCoeff() + log(l_array(i,i,k));	// Choose most likely transition
			}
		}

		Eigen::Index max_row;
		Eigen::Index max_col;
		Type foo;

		// Get the last column of viterbi matrix for backwards probabilities
		matrix<Type> column = viterbi.col(max-1);  // column matrix
		foo = column.maxCoeff(&max_row, &max_col);	// Sends index of largest coefficient to max_row (don't care about max_col)
		viterbi_path(max-1) = max_row + 1;	// Eigen::Index starts at 0, want it to start at 1

		// Backwards recursion
		for(int k=max-2; k>=min; k--) {
			for(int i=0; i<n_state; i++) {
				column(i,0) = viterbi(i,k) + log(tpm(i,max_row));	// backwards Viterbi coefficient
			}
			foo = column.maxCoeff(&max_row, &max_col);	// sends index of largest coefficient to max_row
			viterbi_path(k) = max_row + 1;	// Eigen::Index starts at 0, want it to start at 1
		}
	}
	grp = 0;
	//
	REPORT(viterbi_path);
	//




	//
	//
	// Likelihood computation, OSA forecast predictions
	//
	//
	//array<Type> l_array = theta_array * dist_array; // declared earlier for Viterbi algorithm
	Type ll = 0.0;

	// Forward Probabilities as row vector
	matrix<Type> alpha(1,n_state);

	// Store forward probabilities
	matrix<Type> forecast(n_state,n_obs);



	// Likelihood for first observation
//	alpha = start_probs_row * l_array.col(0).matrix();
//	ll += log(alpha.sum());	//add log of forward probablity to recursive likelihood
//	forecast.col(0) = ( alpha*(tpm.matrix()) ).transpose() / alpha.sum();	// OSA forecast prediction



	for(int k=0; k<n_obs; k++) {
		if( k == grouping(grp) ) {
			forecast.col(k) = start_probs;
			alpha = start_probs_row * l_array.col(k).matrix();
		} else {
			forecast.col(k) = ( alpha*tpm.matrix() ).transpose() / alpha.sum();   // OSA forecast prediction
			alpha = alpha* ( tpm.matrix() ) * ( l_array.col(k).matrix() ); // Add k'th observation to forward probability
		}
		ll += log(alpha.sum());  // add log of forward probability to recursive likelihood
		alpha = alpha/alpha.sum();	// rescale forward probabilities to sum to 1
		if( k == grouping(grp) ) grp++;
	}
	//
	REPORT(forecast);
	return -ll + dummy*dummy;
	//
}
