#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "mtwister.h"


//---------------------------------------------------------------------------------------------------
// VARIABLES
#define _fp 0
#define _fy 1
#define _sow 2
#define _mp 3
#define _my 4
#define _boar 5

#define _m 0
#define _s 1
#define _e 2
#define _i 3
#define _r 4

#define col_idx 0
#define col_sow 1
#define col_nonsow 2
#define col_fid 3

#define wpy 52 // 52 weeks per year
#define w_update 43 // A week to update the weeks for breeding (weaning), farrowing, and dispersal for every cell

#define mu_breed 50.5 // Breeding period: Normal distribution
#define sd_breed 3 // Breeding period: Normal distribution
#define breed_min 44 // Earliest week for breeding period
#define breed_max 5 // Latest week for breeding period

#define gest 17 // Gestation period
#define disp 35.0 // Indicating 8 months after weaning: Dispersal of weaned WBs can be anytime from weaning ~ +8 months 

#define col_breed 0
#define col_farrow 1
#define col_disp 2

#define w_mab 12
#define w_latent 1
#define w_infectious 1

#define p_death 0.99

#define _na 99999
//#define n_grid_total 36186
#define n_grid 5480
#define n_adj 8
#define max_dist 15.0
#define max_reach 145
#define n_particle 2
#define n_chain 5000
#define n_parameter 7
#define Old 0
#define New 1
#define spatial_lim 0.5

#define fwk_1 5
#define fwk_2 6
#define fwk_3 7
#define fwk_4 8
#define fwk_5 11
#define fwk_6 13
#define fwk_7 14
#define fwk_8 18
#define fwk_9 19
#define fwk_10 22

#define week_index 39
#define cid_00 262
#define cid_01 167
#define cid_02 168
#define cid_03 578
#define cid_04 849
#define cid_05 1220
#define cid_06 579
#define cid_07 1716
#define cid_08 7
#define first_wk 0 // One week previous to the index case which is on Week 1
#define second_wk 1
#define third_wk 2
#define fourth_wk 3
#define fifth_wk 27

#define d_fence 1

const int g_demo= 6; // 0: Female piglet, 1: Female yearling, 2: Sow, 3: Male piglet, 4: Male yearling, 5: Male boar
const int g_dz= 5; // 0: Maternal antibody, 1: Susceptible, 2: Exposed, 3: Infectious, 4: Recovered
const int w_SS= 40;
const int n_row_road_cross= 171588;
const int n_row_river_cross= 58915;

char grid_info_file[]= "D:/ASF Kor/grid_info.csv";
char grid_adj_file[]= "D:/ASF Kor/adj_cid.csv";
char grid_xy_file[]= "D:/ASF Kor/grid_xy.csv";
char SS_n_file[]= "D:/ASF Kor/SS_n.csv";
char SS_k_file[]= "D:/ASF Kor/SS_k.csv";
char SS_n_h_file[]= "D:/ASF Kor/SS_n_h.csv";
char SS_k_h_file[]= "D:/ASF Kor/SS_k_h.csv";
char SS_grid_file[]= "D:/ASF Kor/SS_grid.csv";
char road_cross_file[]= "D:/ASF Kor/road_cross.csv";
char town_cross_file[]= "D:/ASF Kor/town_cross.csv";
char river_cross_file[]= "D:/ASF Kor/river_cross.csv";
char para_file[]= "D:/ASF Kor/Result.csv";


//---------------------------------------------------------------------------------------------------
// FUNCTIONS
double _Binomial_function(int k, int n, double p);
//double _Relative_binfunc(int k, int n, double p);
void _Random_order(int *seq, int objects);
int _Random_binomial(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, int n, double p);
double _Random_normal(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double mu, double sd);
double _Random_normal_simple(double mu, double sd);
double _PDF_normal(double x, double mu, double sd);
double _CDF_normal(double x, double mu, double sd);
//double _Random_exp(double lambda);
double _Random_uniform(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double a, double b);
double _Cross_cell(int **grid_info, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid, int adj_cid);
int _Random_grid(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, 
	int **grid_info, int **adj_cid, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid);
int _Natal_migration(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, 
	int **grid_info, int **adj_cid, int **n_wb, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid, int sex);
//void _Initial_wb(int **grid_info, int **n_wb);
double _Dist_grid(double **grid_xy, int i, int j);
void _Introduce_ASF(int week_cont, int **n_wb);
void _LL_estimation(int **n_wb, int **carcass, int **SS_n, int **SS_k, int **SS_n_h, int **SS_k_h, int week_cont, int current_cid, double *log_lik);	
void _Simulation(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
	int **grid_info, int **adj_cid, double **grid_xy, int **reach_cid, double **reach_dist, int **n_wb, int **n_wb_tmp, int **carcass, int **carcass_tmp, int **w_demo, 
	int **SS_n, int **SS_k, int **SS_grid, int **SS_n_h, int **SS_k_h, int **road_cross, int **town_cross, int **river_cross,
	int week_cont, double *beta_wb, double *beta_car, double *pi_between, double *p_road, double *p_town, double *p_river, double *p_fence, double *alpha, double *log_lik);
void _Read_init_wb(char file_name[], int **grid_info, int **n_wb, int **carcass);
void _Read_SS_grid(char file_name[], int **SS_grid);
double _Run_the_model(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
	int **grid_info, int **adj_cid, double **grid_xy, int **reach_cid, double **reach_dist, int **SS_n, int **SS_k, int **SS_n_h, int **SS_k_h, 
	int **road_cross, int **town_cross, int **river_cross, 
	double *beta_wb, double *beta_car, double *pi_between, double *p_road, double *p_town, double *p_river, double *p_fence, double *alpha);
void _Propose_parameters(double **_parameters, int para_id);
double _Selection_prob(double **_parameters, double *loglikelihood, int para_id);
void _Read_grid_info(char file_name[], int **grid_info);
void _Read_SS(char file_name[], int **SS);
void _Read_SS_grid(char file_name[], int **SS_grid);
void _Write_init_wb(char file_name[], int **n_wb, int **carcass);
void _Read_init_wb(char file_name[], int **grid_info, int **n_wb, int **carcass);
void _Read_grid_adj(char file_name[], int **adj_cid);
void _Read_grid_xy(char file_name[], double **grid_xy);
void _Read_grid_cross(char file_name[], int **grid_cross, int n_row);
void _Write_parameters(char file_name[], double **parameters);


//---------------------------------------------------------------------------------------------------
// MAIN FUNCTION
int main(void) {
	int l, m, n;
	double dist, prob;

	// Declare arrays for the simulation
	int **adj_cid;
	int **grid_info; // Grid information
	double **grid_xy; // Grid centroid coordinates
	int **reach_cid;
	double **reach_dist;
	int **SS_n; // Summary statistics (number of tested per cell)
	int **SS_k; // Summary statistics (number of test positives per cell)
	int **SS_n_h; // Summary statistics (number of tested per cell)
	int **SS_k_h; // Summary statistics (number of test positives per cell)
	int **road_cross;
	int **town_cross;
	int **river_cross;
	adj_cid= calloc(n_grid, sizeof(int *));
	grid_info= calloc(n_grid, sizeof(int *));
	grid_xy= calloc(n_grid, sizeof(double *));
	reach_cid= calloc(n_grid, sizeof(int *));
	reach_dist= calloc(n_grid, sizeof(double *));
	SS_n= calloc(n_grid, sizeof(int *));
	SS_k= calloc(n_grid, sizeof(int *));
	SS_n_h= calloc(n_grid, sizeof(int *));
	SS_k_h= calloc(n_grid, sizeof(int *));
	road_cross= calloc(n_grid, sizeof(int *));
	town_cross= calloc(n_grid, sizeof(int *));
	river_cross= calloc(n_grid, sizeof(int *));
	for (l= 0; l < n_grid; ++l) {
		adj_cid[l]= calloc(n_adj, sizeof(int));
		grid_info[l]= calloc(4, sizeof(int)); // 0: max number of sows, 1: max number of non-sows, 2: inside fence
		grid_xy[l]= calloc(2, sizeof(double)); // 0: x, 1: y
		reach_cid[l]= calloc(max_reach, sizeof(int));
		reach_dist[l]= calloc(max_reach, sizeof(double));
		SS_n[l]= calloc(w_SS, sizeof(int));
		SS_k[l]= calloc(w_SS, sizeof(int));
		SS_n_h[l]= calloc(w_SS, sizeof(int));
		SS_k_h[l]= calloc(w_SS, sizeof(int));
		road_cross[l]= calloc(n_grid, sizeof(int));
		town_cross[l]= calloc(n_grid, sizeof(int));
		river_cross[l]= calloc(n_grid, sizeof(int));
	}
	for (l= 0; l < n_grid; ++l) {
		for (m= 0; m < n_grid; ++m) {
			road_cross[l][m]= 0;
			town_cross[l][m]= 0;
			river_cross[l][m]= 0;
		}
	}

	// Import a file for grid information
	_Read_grid_adj(grid_adj_file, adj_cid);
	_Read_grid_info(grid_info_file, grid_info);
	_Read_grid_xy(grid_xy_file, grid_xy);
	_Read_SS(SS_n_file, SS_n);
	_Read_SS(SS_k_file, SS_k);
	_Read_SS(SS_n_h_file, SS_n_h);
	_Read_SS(SS_k_h_file, SS_k_h);
	_Read_grid_cross(road_cross_file, road_cross, n_row_road_cross);
//	_Read_grid_cross(town_cross_file, town_cross);
	_Read_grid_cross(river_cross_file, river_cross, n_row_river_cross);

	// Store the cids and its relative distance 
	for (l= 0; l < n_grid; ++l) {
		for (m= 0; m < max_reach; ++m) {
			reach_cid[l][m]= _na;
		}
	}
	
	for (l= 0; l < n_grid; ++l) {
		n= 0;
		for (m= 0; m < n_grid; ++m) {
			dist= _Dist_grid(grid_xy, l, m);
			if (dist <= max_dist && l != m) {
				reach_cid[l][n]= m;
				reach_dist[l][n]= dist;
				n += 1;
			}
		}
	}

	double *loglikelihood;
	loglikelihood= calloc(n_particle, sizeof(double));

	// Parameters
	// To store
	double **parameters;
	parameters= calloc(n_chain * n_parameter + 1, sizeof(double *));
	for (l= 0; l < n_chain * n_parameter + 1; ++l) {
		parameters[l]= calloc(n_parameter, sizeof(double)); // Number of parameters to estimate
	}
	// Temperory
	double **_parameters;
	_parameters= calloc(n_particle, sizeof(double *));
	for (l= 0; l < n_particle; ++l) {
		_parameters[l]= calloc(n_parameter, sizeof(double));
	}
	
	// Set the initial values of the parameters

	parameters[0][0]= 0.01;//0.01 : Beta_wb
	parameters[0][1]= 0.01;//0.01 : Beta_car
	parameters[0][2]= 0.5;//0.5 : alpha
	parameters[0][3]= 0.1;//0.1 : pi_between
	parameters[0][4]= 0.1;//0.1 : p_road
	parameters[0][5]= 0.1;//0.1 : p_river
	parameters[0][6]= 0.1;//0.1 : p_fence
/*
	parameters[0][0]= 0.1;//0.01 : Beta_wb
	parameters[0][1]= 0.1;//0.01 : Beta_car
	parameters[0][2]= 2.0;//0.5 : alpha
	parameters[0][3]= 0.9;//0.1 : pi_between
	parameters[0][4]= 0.9;//0.1 : p_road
	parameters[0][5]= 0.9;//0.1 : p_river
	parameters[0][6]= 0.9;//0.1 : p_fence
*/
	// Start MH algorithm
	for (l= 0; l < n_chain; ++l) {
		for (m= 0; m < n_parameter; ++m) {
			n= n_parameter * l + m + 1; // Number of row in the "parameters"

			// Bring-up the previous parameter values
			_parameters[Old][0]= parameters[n-1][0];
			_parameters[Old][1]= parameters[n-1][1];
			_parameters[Old][2]= parameters[n-1][2];
			_parameters[Old][3]= parameters[n-1][3];
			_parameters[Old][4]= parameters[n-1][4];
			_parameters[Old][5]= parameters[n-1][5];
			_parameters[Old][6]= parameters[n-1][6];

			// Propose a new parameter
			_Propose_parameters(_parameters, m);
//			printf("%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n", _parameters[New][0], _parameters[New][1], _parameters[New][2], _parameters[New][3], _parameters[New][4], _parameters[New][5], _parameters[New][6]);
		
			#pragma omp parallel num_threads(n_particle)
			{
				int val= omp_get_thread_num();
				double ll;
				srand(n);

				// Mersenne Twister variables for each thread
				int _MTcnt, *MTcnt; _MTcnt= 0; MTcnt= &_MTcnt;
				int _MTleft, *MTleft; _MTleft= 1; MTleft= &_MTleft; 
				int _MTinitf, *MTinitf; _MTinitf= 0; MTinitf= &_MTinitf;
				unsigned long MTstate[MTn], MTnext[MTn];
				_MTinit(MTstate, MTleft, MTinitf, n);

				// Parameters
				double *beta_wb;
				beta_wb= &_parameters[val][0];
				double *beta_car;
				beta_car= &_parameters[val][1];
				double *alpha;
				alpha= &_parameters[val][2];
				double *pi_between;
				pi_between= &_parameters[val][3];
				double *p_road;
				p_road= &_parameters[val][4];
				double *p_river;
				p_river= &_parameters[val][5];
				double *p_fence;
				p_fence= &_parameters[val][6];
				double *p_town, _p_town;
				_p_town= 1.0;
				p_town= &_p_town;

//				printf("IDX %i: %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n", val, _parameters[val][0], _parameters[val][1], _parameters[val][2], _parameters[val][3], _parameters[val][4], _parameters[val][5], _parameters[val][6]);
				ll= _Run_the_model(MTstate, MTnext, MTleft, MTcnt, 
					grid_info, adj_cid, grid_xy, reach_cid, reach_dist, SS_n, SS_k, SS_n_h, SS_k_h, road_cross, town_cross, river_cross, beta_wb, beta_car, pi_between, p_road, p_town, p_river, p_fence, alpha);
//				printf("IDX %i: %4.3f\n", val, ll);
				loglikelihood[val]= ll;
			}
			printf("%4.3f, %4.3f\n", loglikelihood[0], loglikelihood[1]);

			// Update the parameter to "parameters"
			// Choose to select the proposed parameter = Posterior likelihood + Prior adjusting for the probability due to Hastings algorithm
			prob= _Selection_prob(_parameters, loglikelihood, m);
//			printf("%4.3f, %4.3f by %4.3f\n", loglikelihood[0], loglikelihood[1], prob);
			//prob= exp(loglikelihood[New] - loglikelihood[Old] + log((1.0/_parameters[Old][m])*_PDF_normal(log(_parameters[Old][m]),_parameters[New][m])) - log((1.0/_parameters[New][m])*_PDF_normal(log(_parameters[New][m]),_parameters[Old][m])));
			if ((double)rand()/(double)RAND_MAX < prob && loglikelihood[1] != 0.0) {
				// Keep the proposed parameter
				parameters[n][0]= _parameters[New][0];
				parameters[n][1]= _parameters[New][1];
				parameters[n][2]= _parameters[New][2];
				parameters[n][3]= _parameters[New][3];
				parameters[n][4]= _parameters[New][4];
				parameters[n][5]= _parameters[New][5];
				parameters[n][6]= _parameters[New][6];
				printf("Update!\n");
			} else {
				// Discard the proposed parameter
				parameters[n][0]= _parameters[Old][0];
				parameters[n][1]= _parameters[Old][1];
				parameters[n][2]= _parameters[Old][2];
				parameters[n][3]= _parameters[Old][3];
				parameters[n][4]= _parameters[Old][4];
				parameters[n][5]= _parameters[Old][5];
				parameters[n][6]= _parameters[Old][6];
				printf("Discard...\n");
			}
	
			printf("%5i iterations is over.\n", n);
		}
	}
	
	// Export the result
	_Write_parameters(para_file, parameters);

	// Free the memory
	for (l= 0; l < n_chain * n_parameter + 1; ++l) {
		free(parameters[l]);
	}	
	free(parameters);

	for (l= 0; l < n_particle; ++l) {
		free(_parameters[l]);
	}	
	free(_parameters);
	
	for (l= 0; l < n_grid; ++l) {
		free(adj_cid[l]);
		free(grid_info[l]);
		free(grid_xy[l]);
		free(reach_cid[l]);
		free(reach_dist[l]);
		free(SS_n[l]);
		free(SS_k[l]);
		free(SS_n_h[l]);
		free(SS_k_h[l]);
		free(road_cross[l]);
		free(town_cross[l]);
		free(river_cross[l]);
	}

	free(adj_cid);
	free(grid_info);
	free(grid_xy);
	free(reach_cid);
	free(reach_dist);
	free(SS_n);
	free(SS_k);
	free(SS_n_h);
	free(SS_k_h);
	free(road_cross);
	free(town_cross);
	free(river_cross);

	free(loglikelihood);
	
	return 0;
}


//---------------------------------------------------------------------------------------------------
// FUNCTIONS IN DETAIL
// Binomial function
double _Binomial_function(int k, int n, double p) {
	int i, x;
	double answer= 0.0;
	
	if (k > n) {
		printf("Warning! k is bigger than n in _Binomial_function!\n");
	}
	if (k < 0 || n <= 0) {
		printf("Warning! k < 0 or n < 0 in _Binomial_function!\n");
	}
	
	// When n is 1, the PMF equals to "p" when k is 1, "1-p" when k is 0.
	if (n == 1) {
		if (k == 0) {
			answer= 1 - p;
		} else {
			answer= p;
		}
	} else {
		// Factorial calculation is requied - n! / (k! x (n-k)!) - so define the small one between n or n-k (= x).
		if (k >= (n / 2)) {
			x= (n - k);
		} else {
			x= k;
		}
		
		// If "p" is zero, the only non-zero PMF value is when k is zero.
		if (p == 0.0) {
			if (k == 0) {
				answer= 1.0;
			} else {
				answer= 0.0;
			}
		// If "p" is one, the only non-zero PMF value is when k is n.
		} else if (p == 1.0) {
			if (k == n) {
				answer= 1.0;
			} else {
				answer= 0.0;
			}
		} else {
			// If k equals to n, the right side of p^k * (1-p)^(n-k) is ignored.
			if (k == n) {
				answer= pow(p, k);
			// If k equals to zero, the left side of p^k * (1-p)^(n-k) is ignored.
			} else if (k == 0) {
				answer= pow((1.0 - p), (n - k));
			} else {
				// Calculate the factorial - n! / (k! x (n-k)!) - as a log scale.
				for (i= 0; i < x; ++i) {
					answer += (log(n - i) - log(i + 1)); // This is the same as "(log(n - i) - log(x - i))", but the denominator is just inverse order.
				}
				answer= exp(answer);
				answer *= pow(p, k);
				answer *= pow((1.0 - p), (n - k));
			}
		}
	}
	
	return answer;
}
/*
// Relative % of binomial function value compared to the maximum PMF value
double _Relative_binfunc(int k, int n, double p) {
	double _Binomial_function(int k, int n, double p);
	
	int i;
	double val, max= 0.0;
	
	// Search for the maximum PMF of binomial function
	i= 0;
	while (n >= i) {
		val= _Binomial_function(i, n, p);
		
		if (val >= max) {
			max= val;
			i += 1;
		} else {
			break;	
		}		
	}
	
	if (max == 0) {
		printf("Warning! The denominator of _Relative_binfunc is ZERO!\n");
	}
	val= _Binomial_function(k, n, p);
	val /= max;
	
	return val;
}
*/
// Random order of simulation: WARNING! THIS CODE IS BASED ON 'rand()' WHICH HAS A LIMITATION OF 32767 IN MINGW. 
// MERSENNE TWISTER SHOULD BE USED IN CASE MORE THAN 32767 OBJECTS ARE ORDERED.
void _Random_order(int *seq, int objects) {
	int i, j, x;
	
	for (i= 0; i < objects; ++i) {
		seq[i]= i;
	}

	for (i= objects; i > 0; --i) {
		j= rand() % (i + 1);
		x= seq[i-1];
		seq[i-1]= seq[j];
		seq[j]= x;
	}
}

// Random binomial
int _Random_binomial(unsigned long MTstate[], unsigned long MTnext[], int *MTleft , int *MTcnt, int n, double p) {
    int index;
	int n_pos= 0;
	int i= 0;
	if (p > 1) {p= 1.0;}
	if (p < 0) {p= 0.0;}

	double x;

	if (p == 0.0 || n == 0) {
		return 0;
	} else if (p == 1.0) {
		return n;
	} else {
	    do {
	    	if (p < 1.0/(double)RAND_MAX) {
			    do {
					x= _Rand(MTstate, MTnext, MTleft, MTcnt);
				} while (x > 1.0 || x < 0.0);
//				x= _Rand(MTstate, MTnext, MTleft, MTcnt);
			} else {
				x= (double)rand()/(double)RAND_MAX; // THIS DOES NOT GENEREATE ANY VALUE LESS THAN 1/32767 UNDER MINGW.
			}
			if(p >= x) {index= 1;} else {index= 0;}
			n_pos += index;
			i++;
	    } while (i < n);
	
		return n_pos;
	}
}

// Random number from normal distribution: It retunrs a random number according to a normal distribution ("mu", "sigma").
double _Random_normal(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double mu, double sd) {
	double two_pi= 2.0 * acos(-1.0); // acos(-1) is PI
	double u, v, z;

	do {
		u= _Rand(MTstate, MTnext, MTleft, MTcnt);
		v= _Rand(MTstate, MTnext, MTleft, MTcnt);
		//u= (double)rand()/(double)RAND_MAX; // THIS DOES NOT GENEREATE ANY VALUE LESS THAN 1/32767 UNDER MINGW.
		//v= (double)rand()/(double)RAND_MAX; // THIS DOES NOT GENEREATE ANY VALUE LESS THAN 1/32767 UNDER MINGW.
	} while (u == 0.0);

	z= sqrt(-2.0 * log(u)) * cos(two_pi * v);

	return (mu + z * sd);
}
double _Random_normal_simple(double mu, double sd) {
	double two_pi= 2.0 * acos(-1.0); // acos(-1) is PI
	double u, v, z;

	do {
		u= (double)rand()/(double)RAND_MAX;
		v= (double)rand()/(double)RAND_MAX;
	} while (u == 0.0);

	z= sqrt(-2.0 * log(u)) * cos(two_pi * v);

	return (mu + z * sd);
}

// Probability density function of normal distribution
double _PDF_normal(double x, double mu, double sd) {
	double val;
	val= (1.0 / (sd * pow(2.0 * acos(-1.0), 0.5))) * exp(-0.5 * pow(((x - mu) / sd), 2.0));
	return val;
}

// Cumulative probability density function of normal distribution
double _CDF_normal(double x, double mu, double sd) {
	double val;
	val= 0.5 * (1 + erf((x - mu) / (sd * pow(2.0, 0.5))));
	return val;
}

// Random number between two numbers: It retunrs a random number according to a uniform distribution ("a", "b").
double _Random_uniform(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double a, double b) {
	double small;
	double diff= fabs(b - a);
	double x;

	if (a > b) {small= b;} else if (a < b) {small= a;} else {return a;}

	if (diff > (double)RAND_MAX) {
		x= diff * _Rand(MTstate, MTnext, MTleft, MTcnt);
	} else {
		x= diff * (double)rand()/(double)RAND_MAX; // THIS DOES NOT GENEREATE ANY VALUE LESS THAN 1/32767 UNDER MINGW.
	}
	small= small + x;

	return small;
}
/*
// Random number from exponential distribution with the rate of "lambda"
double _Random_exp(double lambda) {
    double u;
    u= (double)rand() / (double)(RAND_MAX + 1.0);
    return (-log(1 - u) / lambda);
}
*/
// Function to designate the probability of movement between cells depending on the fence status
double _Cross_cell(int **grid_info, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid, int adj_cid) {
	int fid_org, fid_dst;
	double p_export;
	fid_org= grid_info[current_cid][col_fid];
	fid_dst= grid_info[adj_cid][col_fid];
	
	if (week_cont < fwk_1) {
		fid_org= 0;
		fid_dst= 0;
	} else if (week_cont < fwk_2) {
		if (fid_org > 2) {fid_org= 0;}
		if (fid_dst > 2) {fid_dst= 0;}
	} else if (week_cont < fwk_3) {
		if (fid_org > 4) {fid_org= 0;}
		if (fid_dst > 4) {fid_dst= 0;}
	} else if (week_cont < fwk_4) {
		if (fid_org > 6) {fid_org= 0;}
		if (fid_dst > 6) {fid_dst= 0;}
	} else if (week_cont < fwk_5) {
		if (fid_org > 7) {fid_org= 0;}
		if (fid_dst > 7) {fid_dst= 0;}
	} else if (week_cont < fwk_6) {
		if (fid_org == 16) {fid_org= 0;} else if (fid_org > 9) {fid_org= 8;}
		if (fid_dst == 16) {fid_dst= 0;} else if (fid_dst > 9) {fid_dst= 8;}
	} else if (week_cont < fwk_7) {
		if (fid_org == 16) {fid_org= 0;} else if (fid_org > 11) {fid_org= 8;}
		if (fid_dst == 16) {fid_dst= 0;} else if (fid_dst > 11) {fid_dst= 8;}
	} else if (week_cont < fwk_8) {
		if (fid_org == 16) {fid_org= 0;} else if (fid_org > 12) {fid_org= 8;}
		if (fid_dst == 16) {fid_dst= 0;} else if (fid_dst > 12) {fid_dst= 8;}
	} else if (week_cont < fwk_9) {
		if (fid_org == 16) {fid_org= 0;} else if (fid_org > 13) {fid_org= 8;}
		if (fid_dst == 16) {fid_dst= 0;} else if (fid_dst > 13) {fid_dst= 8;}
	} else if (week_cont < fwk_10) {
		if (fid_org == 16) {fid_org= 0;} else if (fid_org > 15) {fid_org= 8;}
		if (fid_dst == 16) {fid_dst= 0;} else if (fid_dst > 15) {fid_dst= 8;}
	} else {
		if (fid_org > 16) {fid_org= 0;}
		if (fid_dst > 16) {fid_dst= 0;}
	}

	// Fence ON
	if (d_fence == 1) {
		if (fid_org != fid_dst) {p_export= (*p_fence);} else {p_export= 1.0;}
	} else {
		p_export= 1.0;
	}

	if (road_cross[current_cid][adj_cid] == 1) {
		p_export *= (*p_road);
	}
	
	if (town_cross[current_cid][adj_cid] == 1) {
		p_export *= (*p_town);
	}

	if (river_cross[current_cid][adj_cid] == 1) {
		p_export *= (*p_river);
	}
	
	return p_export;
}

// Random choose
int _Random_grid(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, 
	int **grid_info, int **adj_cid, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid) {
	
	void _Random_order(int *seq, int objects);
	double _Cross_cell(int **grid_info, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid, int adj_cid);
	
	int i, j, dst, *seq;
	double r, prob;
	seq= calloc(n_adj, sizeof(int));
	_Random_order(seq, n_adj);

	dst= _na;
	for (i= 0; i < n_adj; ++i) {
		j= adj_cid[current_cid][seq[i]];
		if (j != _na) {
			prob= _Cross_cell(grid_info, road_cross, town_cross, river_cross, p_road, p_town, p_river, p_fence, week_cont, current_cid, j);
		} else {
			prob= 0.0;
		}

		do {
			r= _Rand(MTstate, MTnext, MTleft, MTcnt);
			//r= (double)rand()/(double)RAND_MAX; // THIS DOES NOT GENEREATE ANY VALUE LESS THAN 1/32767 UNDER MINGW.
		} while (r == 0.0 || r == 1.0);
		
		if (r <= prob) {
			dst= j;
			break;
		}
	}
	free(seq);
	
	if (dst == _na) {
		return current_cid;
	} else {
		return dst;
	}
}
/*
// Initialise the wild boar population
void _Initial_wb(int **grid_info, int **n_wb) {
	int i;
	// For each cell in the grid
	for (i= 0; i < n_grid; ++i) {
		// Is the cell on land? Run only those with 1 (= on land)
		if (grid_info[i][0] == 1) {
			n_wb[i][g_dz * _s + _fy]= grid_info[i][1];
			n_wb[i][g_dz * _s + _sow]= grid_info[i][1];
			n_wb[i][g_dz * _s + _boar]= grid_info[i][1];
		}
	}
}
*/
// Natal migration
int _Natal_migration(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, 
	int **grid_info, int **adj_cid, int **n_wb, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid, int sex) {

	void _Random_order(int *seq, int objects);
	double _Cross_cell(int **grid_info, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid, int adj_cid);
	int _Random_grid(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, 
		int **grid_info, int **adj_cid, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid);
	
	int i, d_lim, step, tmp_cid, tmp_total, done;
	
	if (sex == 0) {
		d_lim= 12; // Up to 25 km
	} else {
		d_lim= 24; // Up to 50 km
	}
	
	step= 0; tmp_cid= current_cid;
	
	// Group splitting of female yearlings
	if (sex == 0) {
		while (step < d_lim) {
			tmp_cid= _Random_grid(MTstate, MTnext, MTleft, MTcnt, grid_info, adj_cid, road_cross, town_cross, river_cross, p_road, p_town, p_river, p_fence, week_cont, tmp_cid);

			tmp_total= 0;
			for (i= 0; i < g_dz; ++i) {
				tmp_total += (n_wb[tmp_cid][g_demo * i + _fp] + n_wb[tmp_cid][g_demo * i + _fy] + n_wb[tmp_cid][g_demo * i + _sow]);
			}
			
			if (tmp_total == 0 && grid_info[tmp_cid][col_sow] > 0 && tmp_cid != current_cid) {
				step= d_lim;
				done= 1;	
			} else {
				step += 1;
				done= 0;
			}

			if (done == 1) {break;}
		}
	// Migration of male yearlings
	} else {
		tmp_total= 0;
		for (i= 0; i < g_dz; ++i) {
			tmp_total += (n_wb[tmp_cid][g_demo * i + _fp] + n_wb[tmp_cid][g_demo * i + _fy] + n_wb[tmp_cid][g_demo * i + _sow] + n_wb[tmp_cid][g_demo * i + _my] + n_wb[tmp_cid][g_demo * i + _mp] + n_wb[tmp_cid][g_demo * i + _boar]);
		}
		
		if (tmp_total >= (grid_info[tmp_cid][col_sow] + grid_info[tmp_cid][col_nonsow])) {
			while (step < d_lim) {
				tmp_cid= _Random_grid(MTstate, MTnext, MTleft, MTcnt, grid_info, adj_cid, road_cross, town_cross, river_cross, p_road, p_town, p_river, p_fence, week_cont, tmp_cid);
				
				// End the process in case there is no adjacent cells
				if (step == 0 && tmp_cid == current_cid) {break;}
				
				tmp_total= 0;
				for (i= 0; i < g_dz; ++i) {
					tmp_total += (n_wb[tmp_cid][g_demo * i + _fp] + n_wb[tmp_cid][g_demo * i + _fy] + n_wb[tmp_cid][g_demo * i + _sow] + n_wb[tmp_cid][g_demo * i + _my] + n_wb[tmp_cid][g_demo * i + _mp] + n_wb[tmp_cid][g_demo * i + _boar]);
				}

				if (tmp_total < (grid_info[tmp_cid][col_sow] + grid_info[tmp_cid][col_nonsow]) && tmp_cid != current_cid) {
					step= d_lim;
					done= 1;
				} else {
					step += 1;
					done= 0;
				}
				
				if (done == 1) {break;}
			}
		}
	}
		
	if (sex == 0 && done == 0) {
		return current_cid;
	} else {
		return tmp_cid;
	}
}

// Distance between grid cells
double _Dist_grid(double **grid_xy, int i, int j) {
	double val;
	val= pow(pow(grid_xy[i][0] - grid_xy[j][0], 2.0) + pow(grid_xy[i][1] - grid_xy[j][1], 2.0), 0.5);
	val /= 1000.0;
	
	return val;
}

void _Introduce_ASF(int week_cont, int **n_wb) {
	int n_int= 5;
	int cid;
	// Introduce ASF
	if (week_cont == first_wk) {
		cid= cid_00;
		if (n_wb[cid][g_demo * _s + _boar] > n_int) {
			n_wb[cid][g_demo * _i + _boar] += n_int; n_wb[cid][g_demo * _s + _boar] -= n_int;
		} else {
			n_wb[cid][g_demo * _i + _boar] += n_wb[cid][g_demo * _s + _boar]; n_wb[cid][g_demo * _s + _boar]= 0;
		}
	} else if (week_cont == second_wk) {
		cid= cid_01;
		if (n_wb[cid][g_demo * _s + _boar] > n_int) {
			n_wb[cid][g_demo * _i + _boar] += n_int; n_wb[cid][g_demo * _s + _boar] -= n_int;
		} else {
			n_wb[cid][g_demo * _i + _boar] += n_wb[cid][g_demo * _s + _boar]; n_wb[cid][g_demo * _s + _boar]= 0;
		}
	} else if (week_cont == third_wk) {
		cid= cid_02;
		if (n_wb[cid][g_demo * _s + _boar] > n_int) {
			n_wb[cid][g_demo * _i + _boar] += n_int; n_wb[cid][g_demo * _s + _boar] -= n_int;
		} else {
			n_wb[cid][g_demo * _i + _boar] += n_wb[cid][g_demo * _s + _boar]; n_wb[cid][g_demo * _s + _boar]= 0;
		}
		cid= cid_03;
		if (n_wb[cid][g_demo * _s + _boar] > n_int) {
			n_wb[cid][g_demo * _i + _boar] += n_int; n_wb[cid][g_demo * _s + _boar] -= n_int;
		} else {
			n_wb[cid][g_demo * _i + _boar] += n_wb[cid][g_demo * _s + _boar]; n_wb[cid][g_demo * _s + _boar]= 0;
		}
		cid= cid_04;
		if (n_wb[cid][g_demo * _s + _boar] > n_int) {
			n_wb[cid][g_demo * _i + _boar] += n_int; n_wb[cid][g_demo * _s + _boar] -= n_int;
		} else {
			n_wb[cid][g_demo * _i + _boar] += n_wb[cid][g_demo * _s + _boar]; n_wb[cid][g_demo * _s + _boar]= 0;
		}
		cid= cid_05;
		if (n_wb[cid][g_demo * _s + _boar] > n_int) {
			n_wb[cid][g_demo * _i + _boar] += n_int; n_wb[cid][g_demo * _s + _boar] -= n_int;
		} else {
			n_wb[cid][g_demo * _i + _boar] += n_wb[cid][g_demo * _s + _boar]; n_wb[cid][g_demo * _s + _boar]= 0;
		}
	} else if (week_cont == fourth_wk) {
		cid= cid_06;
		if (n_wb[cid][g_demo * _s + _boar] > n_int) {
			n_wb[cid][g_demo * _i + _boar] += n_int; n_wb[cid][g_demo * _s + _boar] -= n_int;
		} else {
			n_wb[cid][g_demo * _i + _boar] += n_wb[cid][g_demo * _s + _boar]; n_wb[cid][g_demo * _s + _boar]= 0;
		}
		cid= cid_07;
		if (n_wb[cid][g_demo * _s + _boar] > n_int) {
			n_wb[cid][g_demo * _i + _boar] += n_int; n_wb[cid][g_demo * _s + _boar] -= n_int;
		} else {
			n_wb[cid][g_demo * _i + _boar] += n_wb[cid][g_demo * _s + _boar]; n_wb[cid][g_demo * _s + _boar]= 0;
		}
	} else if (week_cont == fifth_wk) {
		cid= cid_08;
		if (n_wb[cid][g_demo * _s + _boar] > n_int) {
			n_wb[cid][g_demo * _i + _boar] += n_int; n_wb[cid][g_demo * _s + _boar] -= n_int;
		} else {
			n_wb[cid][g_demo * _i + _boar] += n_wb[cid][g_demo * _s + _boar]; n_wb[cid][g_demo * _s + _boar]= 0;
		}
	}
}

// ASF simulation
void _Simulation(unsigned long *MTstate, unsigned long *MTnext, int *MTleft, int *MTcnt,
	int **grid_info, int **adj_cid, double **grid_xy, int **reach_cid, double **reach_dist, int **n_wb, int **n_wb_tmp, int **carcass, int **carcass_tmp, int **w_demo, 
	int **SS_n, int **SS_k, int **SS_grid, int **SS_n_h, int **SS_k_h, int **road_cross, int **town_cross, int **river_cross,
	int week_cont, double *beta_wb, double *beta_car, double *pi_between, double *p_road, double *p_town, double *p_river, double *p_fence, double *alpha, double *log_lik) {

	double _Binomial_function(int k, int n, double p);
	void _Random_order(int *seq, int objects);
	int _Random_binomial(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, int n, double p);
	double _Random_normal(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double mu, double sd);
	double _Random_normal_simple(double mu, double sd);
	double _Random_uniform(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double a, double b);
	double _Cross_cell(int **grid_info, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid, int adj_cid);
	int _Random_grid(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, 
		int **grid_info, int **adj_cid, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid);
	int _Natal_migration(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, 
		int **grid_info, int **adj_cid, int **n_wb, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid, int sex);
	double _Dist_grid(double **grid_xy, int i, int j);
	void _Introduce_ASF(int week_cont, int **n_wb);
	void _LL_estimation(int **n_wb, int **carcass, int **SS_n, int **SS_k, int **SS_n_h, int **SS_k_h, int week_cont, int current_cid, double *log_lik);

	int i, j, k;
	int val;
	int n_fp, n_fy, n_sow, n_mp, n_my, n_boar, new_cid, n_move, new_born;
	int	inf_sounder, inf_car;
	double val_d;
	double adj_sounder, adj_car, pi;
	double p_mort, p_move, p_repro, lambda;
	int dM, dS, dE, dI, dR, dC;
//	int safeguard= 0;

	// Calendar week for simulation
	int week;
	if (d_fence == 0) {
		week= week_cont % 52;
	} else {
		week= (week_cont + week_index) % 52;
	}
	
	// Decay rate of carcass
	double w_carcass, w_asf_carcass;
	w_carcass= 22.5 * cos(acos(-1.0) * (week - 1) / 26) + 29.5;
	w_carcass= -1.0 * log(1.0 - 0.99) / w_carcass;
	w_asf_carcass= 25.5 * cos(acos(-1.0) * (week - 1) / 26) + 26.5;
	w_asf_carcass= -1.0 * log(1.0 - 0.99) / w_asf_carcass;

	// For each grid
	if (week == w_update || week_cont == 0) {
		for (i= 0; i < n_grid; ++i) {
			// Determine the farrowing and weaning date at the week 43 which is before the breeding period
			do {
				val= (int)round(_Random_normal_simple(50.5, 3.0));
				val %= 52;
			} while (val < breed_min && val > breed_max);
			w_demo[i][col_breed]= val; // Week of breeding & weaning
			w_demo[i][col_farrow]= (w_demo[i][col_breed] + gest) % 52; // Week of farrowing = (week of breeding + gestation period) % 52
	
			do {
				val= (w_demo[i][0] + (int)round(_Random_uniform(MTstate, MTnext, MTleft, MTcnt, 0.0, disp))) % 52;
			} while (val == w_update);
			w_demo[i][col_disp]= val; // Week of dispersal
		}
	}
	
	// ASF introduction
	_Introduce_ASF(week_cont, n_wb);

	// For each cell in the grid
	for (j= 0; j < n_grid; ++j) {
		// Demographic
		n_fp= n_wb[j][g_demo * _m + _fp] + n_wb[j][g_demo * _s + _fp] + n_wb[j][g_demo * _e + _fp] + n_wb[j][g_demo * _i + _fp] + n_wb[j][g_demo * _r + _fp];
		n_fy= n_wb[j][g_demo * _m + _fy] + n_wb[j][g_demo * _s + _fy] + n_wb[j][g_demo * _e + _fy] + n_wb[j][g_demo * _i + _fy] + n_wb[j][g_demo * _r + _fy];
		n_sow= n_wb[j][g_demo * _m + _sow] + n_wb[j][g_demo * _s + _sow] + n_wb[j][g_demo * _e + _sow] + n_wb[j][g_demo * _i + _sow] + n_wb[j][g_demo * _r + _sow];
		n_mp= n_wb[j][g_demo * _m + _mp] + n_wb[j][g_demo * _s + _mp] + n_wb[j][g_demo * _e + _mp] + n_wb[j][g_demo * _i + _mp] + n_wb[j][g_demo * _r + _mp];
		n_my= n_wb[j][g_demo * _m + _my] + n_wb[j][g_demo * _s + _my] + n_wb[j][g_demo * _e + _my] + n_wb[j][g_demo * _i + _my] + n_wb[j][g_demo * _r + _my];
		n_boar= n_wb[j][g_demo * _m + _boar] + n_wb[j][g_demo * _s + _boar] + n_wb[j][g_demo * _e + _boar] + n_wb[j][g_demo * _i + _boar] + n_wb[j][g_demo * _r + _boar];

		// DEMO: For each ASF disease group
		for (k= 0; k < g_dz; ++k) {
			// Carcass degradation
			val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, carcass[j][k], w_carcass);
			carcass_tmp[j][k]= carcass[j][k] - val;

			if (k == _i) {
				val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, carcass_tmp[j][k], w_asf_carcass);
				carcass_tmp[j][_i]-= val;
				carcass_tmp[j][_r]+= val;
			}

			// Mortality
			p_mort= 0.014735;
			val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb[j][g_demo * k + _fp], p_mort);
			n_wb_tmp[j][g_demo * k + _fp]= n_wb[j][g_demo * k + _fp] - val; n_fp -= val; carcass_tmp[j][k] += val;

			val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb[j][g_demo * k + _mp], p_mort);
			n_wb_tmp[j][g_demo * k + _mp]= n_wb[j][g_demo * k + _mp] - val; n_mp -= val; carcass_tmp[j][k] += val;

			p_mort= 0.008358;
			val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb[j][g_demo * k + _fy], p_mort);
			n_wb_tmp[j][g_demo * k + _fy]= n_wb[j][g_demo * k + _fy] - val; n_fy -= val; carcass_tmp[j][k] += val;

			val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb[j][g_demo * k + _sow], p_mort);
			n_wb_tmp[j][g_demo * k + _sow]= n_wb[j][g_demo * k + _sow] - val; n_sow -= val; carcass_tmp[j][k] += val;

			val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb[j][g_demo * k + _my], p_mort);
			n_wb_tmp[j][g_demo * k + _my]= n_wb[j][g_demo * k + _my] - val; n_my -= val; carcass_tmp[j][k] += val;

			val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb[j][g_demo * k + _boar], p_mort);
			n_wb_tmp[j][g_demo * k + _boar]= n_wb[j][g_demo * k + _boar] - val; n_boar -= val; carcass_tmp[j][k] += val;

			// Hunting
			if (d_fence == 1) {
				if (grid_info[j][col_fid] > 0) {
					p_mort= 0.017710;
				} else {
					p_mort= 0.005296;
				}
				val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * k + _fp], p_mort);
				n_wb_tmp[j][g_demo * k + _fp] -= val; n_fp -= val;

				val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * k + _fy], p_mort);
				n_wb_tmp[j][g_demo * k + _fy] -= val; n_fy -= val;

				val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * k + _sow], p_mort);
				n_wb_tmp[j][g_demo * k + _sow] -= val; n_sow -= val;

				val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * k + _mp], p_mort);
				n_wb_tmp[j][g_demo * k + _mp] -= val; n_mp -= val;

				val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * k + _my], p_mort);
				n_wb_tmp[j][g_demo * k + _my] -= val; n_my -= val;

				val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * k + _boar], p_mort);
				n_wb_tmp[j][g_demo * k + _boar] -= val; n_boar -= val;
			} else {
				p_mort= 0.002466;
				val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * k + _fy], p_mort);
				n_wb_tmp[j][g_demo * k + _fy] -= val; n_fy -= val;

				val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * k + _sow], p_mort);
				n_wb_tmp[j][g_demo * k + _sow] -= val; n_sow -= val;

				val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * k + _my], p_mort);
				n_wb_tmp[j][g_demo * k + _my] -= val; n_my -= val;

				val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * k + _boar], p_mort);
				n_wb_tmp[j][g_demo * k + _boar] -= val; n_boar -= val;
			}

			// Weaning
			if (week == w_demo[j][col_breed]) {
				// Female yearling -> Sow
				n_wb_tmp[j][g_demo * k + _sow] += n_wb_tmp[j][g_demo * k + _fy]; n_sow += n_wb_tmp[j][g_demo * k + _fy]; 
				n_fy -= n_wb_tmp[j][g_demo * k + _fy]; n_wb_tmp[j][g_demo * k + _fy]= 0;

				// Female piglet -> Female yearling
				n_wb_tmp[j][g_demo * k + _fy] += n_wb_tmp[j][g_demo * k + _fp]; n_fy += n_wb_tmp[j][g_demo * k + _fp]; 
				n_fp -= n_wb_tmp[j][g_demo * k + _fp]; n_wb_tmp[j][g_demo * k + _fp]= 0;

				// Male yearling -> Male boar
				n_wb_tmp[j][g_demo * k + _boar] += n_wb_tmp[j][g_demo * k + _my]; n_boar += n_wb_tmp[j][g_demo * k + _my]; 
				n_my -= n_wb_tmp[j][g_demo * k + _my]; n_wb_tmp[j][g_demo * k + _my]= 0;

				// Male piglet -> Male yearling
				n_wb_tmp[j][g_demo * k + _my] += n_wb_tmp[j][g_demo * k + _mp]; n_my += n_wb_tmp[j][g_demo * k + _mp]; 
				n_mp -= n_wb_tmp[j][g_demo * k + _mp]; n_wb_tmp[j][g_demo * k + _mp]= 0;
			} 

			// Dispersal
			if (week == w_demo[j][col_disp]) {
				// Female group splitting
				new_cid= _Natal_migration(MTstate, MTnext, MTleft, MTcnt, grid_info, adj_cid, n_wb, road_cross, town_cross, river_cross, p_road, p_town, p_river, p_fence, week_cont, j, 0);
				// Only if the number of sows is more than its breeding capacity AND the available cell is not the current cell
				if (new_cid != j) {
					if (n_sow > grid_info[j][col_sow] && n_fp != 0) {
						p_move= (double)n_wb_tmp[j][g_demo * k + _fp] / (double)n_fp;

						if (n_fp > grid_info[new_cid][col_sow]) {
							n_move= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, grid_info[new_cid][col_sow], p_move);
							n_wb_tmp[new_cid][g_demo * k + _fy] += n_move; n_wb_tmp[j][g_demo * k + _fp] -= n_move; n_fp -= n_move;
						} else {
							n_move= n_wb_tmp[j][g_demo * k + _fp];
							n_wb_tmp[new_cid][g_demo * k + _fy] += n_move; n_wb_tmp[j][g_demo * k + _fp] -= n_move; n_fp -= n_move;
						}
					}
				}
				// Male dispersal
				while (n_wb_tmp[j][g_demo * k + _my] > 0) {
					new_cid= _Natal_migration(MTstate, MTnext, MTleft, MTcnt, grid_info, adj_cid, n_wb, road_cross, town_cross, river_cross, p_road, p_town, p_river, p_fence, week_cont, j, 1);
					n_wb_tmp[new_cid][g_demo * k + _boar] += 1;
					n_wb_tmp[j][g_demo * k + _my] -= 1; n_my -= 1;
				}
			}
		
			// Farrowing
			if (week == w_demo[j][col_farrow]) {
				p_repro= (double)n_sow / (double)grid_info[j][col_sow];

				if (k == _s || k == _r) {
					if (p_repro < 1.0) {
						p_repro= 1.0;
					} else {
						p_repro= (double)grid_info[j][col_sow] / (double)n_sow;
					}
				} else {
					if (p_repro < 1.0) {
						p_repro= 1.0 * 0.625; // Impaired reproduction (Lange et al 2017)
					} else {
						p_repro= 0.625 * (double)grid_info[j][col_sow] / (double)n_sow; // Impaired reproduction (Lange et al 2017)
					}
				}

				new_born= (int)((int)_Random_normal_simple(5.0, 1.0) * _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * k + _sow], p_repro)) / 2; // Male/female newborn piglet
				if (k == _r) {
					n_wb_tmp[j][g_demo * _m + _fp] += new_born; n_fp += new_born;
					n_wb_tmp[j][g_demo * _m + _mp] += new_born; n_mp += new_born;
				} else {
					n_wb_tmp[j][g_demo * _s + _fp] += new_born; n_fp += new_born;
					n_wb_tmp[j][g_demo * _s + _mp] += new_born; n_mp += new_born;
				}
			}
		} // END DEMO: For each ASF disease group

		// ASF
		// Infectious pressure
		inf_sounder= n_wb[j][g_demo * _i + _fp] + n_wb[j][g_demo * _i + _fy] + n_wb[j][g_demo * _i + _sow] + n_wb[j][g_demo * _i + _mp] + n_wb[j][g_demo * _i + _my] + n_wb[j][g_demo * _i + _boar];
		inf_car= carcass[j][_i];
		adj_sounder= 0.0;
		adj_car= 0.0;

		k= 0;
	
		while (reach_cid[j][k] != _na) {
			pi= exp(-1.0 * (*alpha) * reach_dist[j][k]) * _Cross_cell(grid_info, road_cross, town_cross, river_cross, p_road, p_town, p_river, p_fence, week_cont, j, reach_cid[j][k]);;
			adj_sounder += ((double)(n_wb[reach_cid[j][k]][g_demo * _i + _fp] + n_wb[reach_cid[j][k]][g_demo * _i + _fy] + n_wb[reach_cid[j][k]][g_demo * _i + _sow] + n_wb[reach_cid[j][k]][g_demo * _i + _mp] + n_wb[reach_cid[j][k]][g_demo * _i + _my] + n_wb[reach_cid[j][k]][g_demo * _i + _boar]) * pi);
			adj_car += ((double)carcass[reach_cid[j][k]][_i] * pi);
			k += 1;
		}

		lambda= 1.0 - exp(-1.0 * (((*beta_wb) * (double)inf_sounder) + ((*beta_wb) * (*pi_between) * adj_sounder) + ((*beta_car) * (double)inf_car) + ((*beta_car) * (*pi_between) * adj_car)));
		
		for (k= 0; k < g_demo; ++k) {
			// M -> S
			val_d= -1.0 * log(1.0 - 0.99) / w_mab;
			dM= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * _m + k], val_d);
			// S -> E
			dS= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, n_wb_tmp[j][g_demo * _s + k], lambda);
			// E -> I
			dE= n_wb_tmp[j][g_demo * _e + k];
			// I -> R or death
			dI= n_wb_tmp[j][g_demo * _i + k];
			dC= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, dI, p_death); dR= dI - dC;

			n_wb_tmp[j][g_demo * _m + k] -= dM;
			n_wb_tmp[j][g_demo * _s + k] += dM; n_wb_tmp[j][g_demo * _s + k] -= dS;
			n_wb_tmp[j][g_demo * _e + k] += dS; n_wb_tmp[j][g_demo * _e + k] -= dE;
			n_wb_tmp[j][g_demo * _i + k] += dE; n_wb_tmp[j][g_demo * _i + k] -= dI;
			n_wb_tmp[j][g_demo * _r + k] += dR; carcass_tmp[j][_i] += dC;
	
			if (dC > 0) {
				SS_grid[j][1]= 1;
			}
		}

		// Add the likelihood
		_LL_estimation(n_wb_tmp, carcass_tmp, SS_n, SS_k, SS_n_h, SS_k_h, week_cont, j, log_lik);
	} // END For each cell in the grid

	// Null the log likelihood if there is no infectious animal at the end week of simulation
//	if (week_cont == w_SS - 1 && safeguard <= 25) {
//		(*log_lik)= 0.0;
//	}

	// Update the WB population and carcasses
	for (j= 0; j < n_grid; ++j) {
		for (k= 0; k < (g_demo * g_dz); ++k) {
			n_wb[j][k]= n_wb_tmp[j][k];
		}
		for (k= 0; k < g_dz; ++k) {
			carcass[j][k]= carcass_tmp[j][k];
		}
	}
}

// Estimate the log likelihood
void _LL_estimation(int **n_wb, int **carcass, int **SS_n, int **SS_k, int **SS_n_h, int **SS_k_h, int week_cont, int current_cid, double *log_lik) {
	double _Binomial_function(int k, int n, double p);
	
	int j, k, l, val;
	double val_d;
	
	j= current_cid;
	k= week_cont;

	// Carcasses
	if (SS_n[j][k] > 0) {
	// Among the carcasses
		val= carcass[j][_m] + carcass[j][_s] + carcass[j][_e] + carcass[j][_i] + carcass[j][_r];
		if (val > 0) {
			val_d= (double)carcass[j][_i] / (double)val;
		} else {
			val_d= 0.0;
		}

		if (val == 0) {
//			printf("Warning. No carcasses, while the observation is opposite!! Sim= %3i, Obs= %3i, Pos= %3i\n", val, SS_n[j][k], SS_k[j][k]);
		} else {
			(*log_lik) += _Binomial_function(SS_k[j][k], SS_n[j][k], val_d);
		}
	} else {
//		printf("No carcasses reported in CID %3i at Week %3i\n", current_cid, week_cont);
	}
	
	// Hunting bags
	if (SS_n_h[j][k] > 0) {
		val= 0;
		for (l= 0; l < g_dz; ++l) {
			val += (n_wb[j][g_demo * l + _fp]+n_wb[j][g_demo * l + _fy]+n_wb[j][g_demo * l + _sow]+n_wb[j][g_demo * l + _mp]+n_wb[j][g_demo * l + _my]+n_wb[j][g_demo * l + _boar]);
			if (l == _i) {
				val_d= (double)(n_wb[j][g_demo * l + _fp]+n_wb[j][g_demo * l + _fy]+n_wb[j][g_demo * l + _sow]+n_wb[j][g_demo * l + _mp]+n_wb[j][g_demo * l + _my]+n_wb[j][g_demo * l + _boar]);
			}
		}
		if (val > 0) {
			val_d /= (double)val;
		} else {
			val_d= 0.0;
		}

		if (val == 0) {
//			printf("Warning. No carcasses, while the observation is opposite!! Sim= %3i, Obs= %3i, Pos= %3i\n", val, SS_n[j][k], SS_k[j][k]);
		} else {
			(*log_lik) += _Binomial_function(SS_k_h[j][k], SS_n_h[j][k], val_d);
		}
	} else {
//		printf("No carcasses reported in CID %3i at Week %3i\n", current_cid, week_cont);
	}
}

// Simulate ASF spread and return the log-likelihood value
double _Run_the_model(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
	int **grid_info, int **adj_cid, double **grid_xy, int **reach_cid, double **reach_dist, int **SS_n, int **SS_k, int **SS_n_h, int **SS_k_h, 
	int **road_cross, int **town_cross, int **river_cross, 
	double *beta_wb, double *beta_car, double *pi_between, double *p_road, double *p_town, double *p_river, double *p_fence, double *alpha) {
	
	double _Binomial_function(int k, int n, double p);
	void _Random_order(int *seq, int objects);
	int _Random_binomial(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, int n, double p);
	double _Random_normal(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double mu, double sd);
	double _Random_normal_simple(double mu, double sd);
	double _Random_uniform(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double a, double b);
	double _Cross_cell(int **grid_info, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid, int adj_cid);
	int _Random_grid(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, 
		int **grid_info, int **adj_cid, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid);
	int _Natal_migration(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
		int **grid_info, int **adj_cid, int **n_wb, int **road_cross, int **town_cross, int **river_cross, double *p_road, double *p_town, double *p_river, double *p_fence, int week_cont, int current_cid, int sex);
	double _Dist_grid(double **grid_xy, int i, int j);
	void _Introduce_ASF(int week_cont, int **n_wb);
	void _LL_estimation(int **n_wb, int **carcass, int **SS_n, int **SS_k, int **SS_n_h, int **SS_k_h, int week_cont, int current_cid, double *log_lik);	
	void _Simulation(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
		int **grid_info, int **adj_cid, double **grid_xy, int **reach_cid, double **reach_dist, int **n_wb, int **n_wb_tmp, int **carcass, int **carcass_tmp, int **w_demo, 
		int **SS_n, int **SS_k, int **SS_grid, int **SS_n_h, int **SS_k_h, int **road_cross, int **town_cross, int **river_cross,
		int week_cont, double *beta_wb, double *beta_car, double *pi_between, double *p_road, double *p_town, double *p_river, double *p_fence, double *alpha, double *log_lik);
	void _Read_init_wb(char file_name[], int **grid_info, int **n_wb, int **carcass);
	void _Read_SS_grid(char file_name[], int **SS_grid);

	int i;
	
	int **n_wb; // Number of live WBs per cell
	int **n_wb_tmp; // Number of live WBs per cell
	int **carcass; // Number of wild boar carcasses per cell
	int **carcass_tmp; // Number of wild boar carcasses per cell
	int **w_demo; // Week of breeding and farrowing of WBs for each cell (or sounder)
	int **SS_grid; // Summary statistics (number of test positive cell)

	n_wb= calloc(n_grid, sizeof(int *));
	n_wb_tmp= calloc(n_grid, sizeof(int *));
	carcass= calloc(n_grid, sizeof(int *));
	carcass_tmp= calloc(n_grid, sizeof(int *));
	w_demo= calloc(n_grid, sizeof(int *));
	SS_grid= calloc(n_grid, sizeof(int *));

	for (i= 0; i < n_grid; ++i) {
		n_wb[i]= calloc(g_demo * g_dz, sizeof(int));
		n_wb_tmp[i]= calloc(g_demo * g_dz, sizeof(int));
		carcass[i]= calloc(g_dz, sizeof(int));
		carcass_tmp[i]= calloc(g_dz, sizeof(int));
		w_demo[i]= calloc(3, sizeof(int)); // 0: Breeding week, 1: Farrowing week, 2: Dispersal week - Updated at week 43 (w_update) which is one week before November 1st
		SS_grid[i]= calloc(2, sizeof(int));
	}
	
	// Read the initial population
	char input_file[100];
	snprintf(input_file, 100, "D:/ASF Kor/n_wbs.csv");
	_Read_init_wb(input_file, grid_info, n_wb, carcass);
	_Read_SS_grid(SS_grid_file, SS_grid);

	// Parameters
	double _log_lik, *log_lik;
	_log_lik= 0.0;
	log_lik= &_log_lik;
/*
	int k, n_m= 0, n_s= 0, n_e= 0, n_i= 0, n_r= 0, carca= 0;
	for (k= 0; k < n_grid; ++k) {
		n_m += n_wb[k][g_demo * _m + _fp] + n_wb[k][g_demo * _m + _fy] + n_wb[k][g_demo * _m + _sow] + n_wb[k][g_demo * _m + _mp] + n_wb[k][g_demo * _m + _my] + n_wb[k][g_demo * _m + _boar];
		n_s += n_wb[k][g_demo * _s + _fp] + n_wb[k][g_demo * _s + _fy] + n_wb[k][g_demo * _s + _sow] + n_wb[k][g_demo * _s + _mp] + n_wb[k][g_demo * _s + _my] + n_wb[k][g_demo * _s + _boar];
		n_e += n_wb[k][g_demo * _e + _fp] + n_wb[k][g_demo * _e + _fy] + n_wb[k][g_demo * _e + _sow] + n_wb[k][g_demo * _e + _mp] + n_wb[k][g_demo * _e + _my] + n_wb[k][g_demo * _e + _boar];
		n_i += n_wb[k][g_demo * _i + _fp] + n_wb[k][g_demo * _i + _fy] + n_wb[k][g_demo * _i + _sow] + n_wb[k][g_demo * _i + _mp] + n_wb[k][g_demo * _i + _my] + n_wb[k][g_demo * _i + _boar];
		n_r += n_wb[k][g_demo * _r + _fp] + n_wb[k][g_demo * _r + _fy] + n_wb[k][g_demo * _r + _sow] + n_wb[k][g_demo * _r + _mp] + n_wb[k][g_demo * _r + _my] + n_wb[k][g_demo * _r + _boar];
		carca += carcass[k][_i];
	}
	printf("Week  0: %6i, %6i, %6i, %6i, %6i, %6i.\n", n_m, n_s, n_e, n_i, n_r, carca);
*/
	// Run the simulation
	for (i= 0; i < w_SS; ++i) {
		_Simulation(MTstate, MTnext, MTleft, MTcnt, grid_info, adj_cid, grid_xy, reach_cid, reach_dist, n_wb, n_wb_tmp, carcass, carcass_tmp, w_demo, SS_n, SS_k, SS_grid, SS_n_h, SS_k_h, road_cross, town_cross, river_cross, i, 
			beta_wb, beta_car, pi_between, p_road, p_town, p_river, p_fence, alpha, log_lik);
/*
		n_m= 0, n_s= 0, n_e= 0, n_i= 0, n_r= 0, carca= 0;
		for (k= 0; k < n_grid; ++k) {
			n_m += n_wb[k][g_demo * _m + _fp] + n_wb[k][g_demo * _m + _fy] + n_wb[k][g_demo * _m + _sow] + n_wb[k][g_demo * _m + _mp] + n_wb[k][g_demo * _m + _my] + n_wb[k][g_demo * _m + _boar];
			n_s += n_wb[k][g_demo * _s + _fp] + n_wb[k][g_demo * _s + _fy] + n_wb[k][g_demo * _s + _sow] + n_wb[k][g_demo * _s + _mp] + n_wb[k][g_demo * _s + _my] + n_wb[k][g_demo * _s + _boar];
			n_e += n_wb[k][g_demo * _e + _fp] + n_wb[k][g_demo * _e + _fy] + n_wb[k][g_demo * _e + _sow] + n_wb[k][g_demo * _e + _mp] + n_wb[k][g_demo * _e + _my] + n_wb[k][g_demo * _e + _boar];
			n_i += n_wb[k][g_demo * _i + _fp] + n_wb[k][g_demo * _i + _fy] + n_wb[k][g_demo * _i + _sow] + n_wb[k][g_demo * _i + _mp] + n_wb[k][g_demo * _i + _my] + n_wb[k][g_demo * _i + _boar];
			n_r += n_wb[k][g_demo * _r + _fp] + n_wb[k][g_demo * _r + _fy] + n_wb[k][g_demo * _r + _sow] + n_wb[k][g_demo * _r + _mp] + n_wb[k][g_demo * _r + _my] + n_wb[k][g_demo * _r + _boar];
			carca += carcass[k][_i];
		}
		printf("Week %2i: %6i, %6i, %6i, %6i, %6i, %6i.\n", i, n_m, n_s, n_e, n_i, n_r, carca);
*/
	}

	// Safeguard
/*
	int j= 0, k= 0;
	double val;
	// Null the log likelihood if cells with ASF+ carcasses do not match with less than x% of observed 127 cells (habitats) where at least one ASF + carcass being reported during the simulation period
	for (i= 0; i < n_grid; ++i) {
		if (SS_grid[i][0] == 1) {
			j +=1;
			if (SS_grid[i][1] == 1) {
				k +=1;
			}
		}
	}
	val= (double)k / (double)j;
//	printf("%2.3f\n", val);
	
	if (val < spatial_lim) {
		_log_lik= 0.0;
	}
*/


	// Null the log likelihood if more than 50% of cells in the simulated area was infected
	int j= 0;
	double val;
	for (i= 0; i < n_grid; ++i) {
		if (SS_grid[i][1] == 1) {
			j +=1;
		}
	}
	val= (double)j / (double)n_grid;
//	printf("%2.3f\n", val);
	
	if (val > spatial_lim) {
		_log_lik= 0.0;
	}

	// Null the log likelihood if there was no simulated ASF + carcasses in any of the cells where three ASF + carcasses (i.e. eastmost, westmost, and farthest from the northern borderline) were observed
	if (SS_grid[1032][1] != 1) {
//	if (SS_grid[1032][1] != 1 || SS_grid[357][1] != 1 || SS_grid[1467][1] != 1) {
		_log_lik= 0.0;
	}

	// Free the memory
	for (i= 0; i < n_grid; ++i) {
		free(n_wb[i]);
		free(n_wb_tmp[i]);
		free(carcass[i]);
		free(carcass_tmp[i]);
		free(w_demo[i]);
		free(SS_grid[i]);
	}
	free(n_wb);
	free(n_wb_tmp);
	free(carcass);
	free(carcass_tmp);
	free(w_demo);
	free(SS_grid);
	
	return (_log_lik);
}

// Proposal function of parameters
void _Propose_parameters(double **_parameters, int para_id) {
	double _Random_normal_simple(double mu, double sd);
	
	double val, sd;
	int i;

	for (i= 0; i < n_parameter; ++i) {
		_parameters[New][i]= _parameters[Old][i];
	}
	
	i= para_id;
	// Betas
	if (i <= 1) {
		sd= 2.0;
		val= exp(_Random_normal_simple(log(_parameters[Old][i]), sd));
		while (val > 1.0 || val < 0.0001) {
			val= exp(_Random_normal_simple(log(_parameters[Old][i]), sd));
		}
	// alpha
	} else if (i == 2) {
		sd= 0.5;
//		sd= 2.0;
		val= _Random_normal_simple(_parameters[Old][i], sd);
//		val= exp(_Random_normal(log(_parameters[Old][i]), sd));
		while (val > 2.5 || val < 0.1) {
			val= _Random_normal_simple(_parameters[Old][i], sd);
//			val= exp(_Random_normal(log(_parameters[Old][i]), sd));
		}
	// pi_between, p_road, p_water, p_fence
	} else {
		sd= 0.25;
//		sd= 2.0;
		val= _Random_normal_simple(_parameters[Old][i], sd);
//		val= exp(_Random_normal(log(_parameters[Old][i]), sd));
		while (val > 0.999 || val < 0.001) {
//		while (val > 0.99 || val < 0.01) {
			val= _Random_normal_simple(_parameters[Old][i], sd);
//			val= exp(_Random_normal(log(_parameters[Old][i]), sd));
		}
	}

	_parameters[New][i]= val;
}	

// Selection probability
double _Selection_prob(double **_parameters, double *loglikelihood, int para_id) {
	double _PDF_normal(double x, double mu, double sd);
	double _CDF_normal(double x, double mu, double sd);
	int i;
	double sd, prob;

	i= para_id;	
	if (i <= 1) {
		sd= 2.0;
		prob= exp(loglikelihood[New] - loglikelihood[Old] + log((1.0/_parameters[Old][i])*_PDF_normal(log(_parameters[Old][i]),_parameters[New][i], sd)) - log((1.0/_parameters[New][i])*_PDF_normal(log(_parameters[New][i]),_parameters[Old][i], sd)));
	} else if (i == 2) {
		sd= 0.5;
		prob= exp(loglikelihood[New] - loglikelihood[Old] + log(_CDF_normal(2.5, _parameters[New][i], sd) - _CDF_normal(0.1, _parameters[New][i], sd)) - log(_CDF_normal(2.5, _parameters[Old][i], sd) - _CDF_normal(0.1, _parameters[Old][i], sd)));
	} else {
		sd= 0.25;
		prob= exp(loglikelihood[New] - loglikelihood[Old] + log(_CDF_normal(0.999, _parameters[New][i], sd) - _CDF_normal(0.001, _parameters[New][i], sd)) - log(_CDF_normal(0.999, _parameters[Old][i], sd) - _CDF_normal(0.001, _parameters[Old][i], sd)));
	}
	
	return prob;
}

/*
// Select the best particle index based on the summary statistics
int _Best_particle(double *log_lik) {
	int _Random_binomial(int n, double p);
	
	int i, j;
	double tmp= log_lik[0];

	i= 0;
	for (j= 1; j < n_particle; ++j) {
		if (log_lik[j] >= tmp) {
			if (log_lik[j] == tmp) {
				if (_Random_binomial(1, 0.5) == 1) {
					tmp= log_lik[j];
					i= j;
				}
			} else {
				tmp= log_lik[j];
				i= j;
			} 
		}
	}
	return i;
}
*/
// Read cell id NA status
void _Read_grid_info(char file_name[], int **grid_info) {
	FILE *_file= fopen(file_name, "r");	
	if (_file == NULL) {printf("No file to import - grid info!\n");}

	int i, idx, nsow, nnonsow, fence;

	for (i= 0; i< n_grid; ++i) {
		fscanf(_file, "%i,%i,%i,%i", &idx, &nsow, &nnonsow, &fence);
		grid_info[i][0]= idx;
		grid_info[i][1]= nsow;
		grid_info[i][2]= nnonsow;
		grid_info[i][3]= fence;
	}

	fclose(_file);
}

// Read summary statistics
void _Read_SS(char file_name[], int **SS) {
	FILE *_file= fopen(file_name, "r");	
	if (_file == NULL) {printf("No file to import - summary statistics!\n");}

	int i, j, val;

	for (i= 0; i< n_grid; ++i) {
		for (j= 0; j < w_SS; ++j) {
			if (j < w_SS - 1) {
				fscanf(_file, "%i,", &val);
				SS[i][j]= val;
			} else {
				fscanf(_file, "%i", &val);
				SS[i][j]= val;
			}
		}
	}

	fclose(_file);
}

void _Read_SS_grid(char file_name[], int **SS_grid) {
	FILE *_file= fopen(file_name, "r");	
	if (_file == NULL) {printf("No file to import - summary statistics!\n");}

	int i, obs, sim;

	for (i= 0; i< n_grid; ++i) {
		fscanf(_file, "%i,%i\n", &obs, &sim);
		SS_grid[i][0]= obs;
		SS_grid[i][1]= sim;
	}

	fclose(_file);
}

// Export the number of WBs in each cell
void _Write_init_wb(char file_name[], int **n_wb, int **carcass) {
	FILE *WBs= fopen(file_name, "w");
	if (WBs == NULL) {printf("No file to export!\n");}

	int i;
	
	for (i= 0; i< n_grid; ++i) {
		fprintf(WBs, "%i,%i,%i,%i,%i,%i,%i\n", 
			n_wb[i][g_demo * _s + _fp], n_wb[i][g_demo * _s + _fy], n_wb[i][g_demo * _s + _sow], n_wb[i][g_demo * _s + _mp], n_wb[i][g_demo * _s + _my], n_wb[i][g_demo * _s + _boar], carcass[i][_s]);
	}

	fclose(WBs);
}

// Read the number of WBs in each cell
void _Read_init_wb(char file_name[], int **grid_info, int **n_wb, int **carcass) {
	FILE *_file= fopen(file_name, "r");
	if (_file == NULL) {printf("No file to import - init wb!\n");}

	int i, n_fp, n_fy, n_sow, n_mp, n_my, n_boar, carca;
	
	for (i= 0; i< n_grid; ++i) {
		fscanf(_file, "%i,%i,%i,%i,%i,%i,%i", &n_fp, &n_fy, &n_sow, &n_mp, &n_my, &n_boar, &carca);
		n_wb[i][g_demo * _s + _fp]= n_fp;
		n_wb[i][g_demo * _s + _fy]= n_fy;
		n_wb[i][g_demo * _s + _sow]= n_sow;
		n_wb[i][g_demo * _s + _mp]= n_mp;
		n_wb[i][g_demo * _s + _my]= n_my;
		n_wb[i][g_demo * _s + _boar]= n_boar;
		carcass[i][_s]= carca;
	}

	fclose(_file);
}

// Read the adjacent cids
void _Read_grid_adj(char file_name[], int **adj_cid) {
	FILE *_file= fopen(file_name, "r");
	if (_file == NULL) {printf("No file to import - adj cid!\n");}

	int i, adj0, adj1, adj2, adj3, adj4, adj5, adj6, adj7;
	for (i= 0; i < n_grid; ++i) {
		fscanf(_file, "%i,%i,%i,%i,%i,%i,%i,%i\n", &adj0, &adj1, &adj2, &adj3, &adj4, &adj5, &adj6, &adj7);
		adj_cid[i][0]= adj0;
		adj_cid[i][1]= adj1;
		adj_cid[i][2]= adj2;
		adj_cid[i][3]= adj3;
		adj_cid[i][4]= adj4;
		adj_cid[i][5]= adj5;
		adj_cid[i][6]= adj6;
		adj_cid[i][7]= adj7;
	}

	fclose(_file);
}

// Read the coordinates of the centroids of grid
void _Read_grid_xy(char file_name[], double **grid_xy) {
	FILE *_file= fopen(file_name, "r");
	if (_file == NULL) {printf("No file to import - grid xy!\n");}

	int i;
	double x, y;
	
	for (i= 0; i< n_grid; ++i) {
		fscanf(_file, "%lf,%lf", &x, &y);
		grid_xy[i][0]= x;
		grid_xy[i][1]= y;
	}

	fclose(_file);
}

void _Read_grid_cross(char file_name[], int **grid_cross, int n_row) {
	FILE *_file= fopen(file_name, "r");
	if (_file == NULL) {printf("No file to import - grid cross!\n");}

	int i, org, dst;
	
	for (i= 0; i< n_row; ++i) {
		fscanf(_file, "%i,%i\n", &org, &dst);
		grid_cross[org][dst]= 1;
		grid_cross[dst][org]= 1;
	}

	fclose(_file);
}

// Export the result
void _Write_parameters(char file_name[], double **parameters) {
	FILE *_file= fopen(file_name, "w");
	if (_file == NULL) {printf("No file to export!\n");}

	int i;
	for (i= 0; i < n_chain * n_parameter + 1; ++i) {
		fprintf(_file, "%f,%f,%f,%f,%f,%f,%f\n", parameters[i][0], parameters[i][1], parameters[i][2], parameters[i][3], parameters[i][4], parameters[i][5], parameters[i][6]);
	}
	
	fclose(_file);
}
