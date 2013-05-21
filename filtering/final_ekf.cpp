#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include "./recipes/nr.h"
#include "filters.h"
#include <vector>
#include "filter_utils.h"
#include <ctime>

using namespace std;

enum VolModel {HESTON = 1, GARCH = 2, THREE_TWO = 3, VAR_P = 4};
enum RunMode {SIMPLE = 1, NORMAL_RESIDUALS = 2, UNCORRELATED_RESIDUALS = 3, NORMAL_AND_UNCORRELATED_RESIDUALS = 4};

int n_stock_prices = 0;
double *log_stock_prices, *u, *v, *estimates;
double muS = 0.0;

double p_heston = 0.50;
double p_garch  = 1.00;
double p_3_2    = 1.50;
double p_variable = 1; //This is the starting case situation where we try to optimize for variable p
int max_simulations = 1;
VolModel model;
RunMode runmode;

int call_counter = 0;
int simulation_counter = 1;

ofstream paramter_file;
ofstream residual_file;
ofstream log_file;

DP ftol = 1.00e-8;


void init_log(const VolModel model, const RunMode runmode, ofstream& log_stream);
string get_run_mode(const RunMode runmode);
string get_model_name(const VolModel model);

void parse_args(int argc, char** argv, 
	string& input_file_name,
	string& parameter_file_name,
	string& residual_file_name);

DP minimize_target_ekf(Vec_I_DP & input);

void log_parameters(double& omega,
	double& theta,
	double& xi,
	double& rho,
	double& p,
	double& mle);

int main(int argc, char** argv) {
	
	srand(time(NULL)); //seed the random number generator

	vector<double> prices;
	string input_file_name, parameter_file_name = "", residual_file_name = "";
	
	parse_args(argc, argv,
		input_file_name,
		parameter_file_name,
		residual_file_name);

	//setting cout to show 10 decimal places by default
	std::cout.setf(ios::fixed);
	std::cout<<setprecision(10); 

	read_lines(input_file_name, prices);
	cout<<"Found "<<prices.size()<<" prices"<<endl;
	
	n_stock_prices = prices.size();

	log_stock_prices = new double[n_stock_prices];
	u = new double[n_stock_prices];
	v = new double[n_stock_prices];
	estimates = new double[n_stock_prices + 1];

	for(unsigned int i = 0; i < prices.size(); i++) {
		log_stock_prices[i] = log(prices[i]);
	}

	//param variables are omega, theta, rho, p, mle, chi_squared_statistic 
	vol_params current_params(0.2, 1.00, 0.50, -0.20, 1.00, pow(10.00, 6), pow(10.00, 6));
	vol_params best_params(current_params), params_buffer(current_params);
	
	Vec_IO_DP* start_;
	
	Mat_IO_DP* identity_matrix;

	if(model!=VAR_P) //if we aren't trying to estimate p, we only need a 4x4 matrix
	{
		//Initializing the starting point
		start_ = new Vec_IO_DP(4);
		identity_matrix = new Mat_IO_DP(4, 4);
	} else { //we are trying to estimate p, so we want a 5x5 matrix
		//Initializing the starting point
		start_ = new Vec_IO_DP(5);
		identity_matrix = new Mat_IO_DP(5, 5);
	}

	Vec_IO_DP& start = *start_; //aliasing the pointer for easy reference
	
	paramter_file.open(parameter_file_name.c_str());
	residual_file.open(residual_file_name.c_str());
	//setting the output to 10 decimal places
	paramter_file.setf(ios::fixed);
	paramter_file<<setprecision(10);
	residual_file.setf(ios::fixed);
	residual_file<<setprecision(10);

	while(simulation_counter <= max_simulations)
	{
		int nvars = (model!=VAR_P) ? 4 : 5; 
		for(int i = 0; i < nvars; i++) 
		{
			for(int j = 0; j < nvars; j++)
			{
				(*identity_matrix)[i][j] = (i==j) ? 1 : 0;
			}
		}


		current_params.populate_starting_vector(start);
		for(int i = 0; i < 5; i++)
			cout<<start[i]<<endl;

		int iter;
		DP mle = 0.00;
		
		NR::powell(start, *identity_matrix, ftol, iter, mle, minimize_target_ekf);
		cout<<"Ran in "<<call_counter<<" iterations with return value "<<mle<<endl;
		
		params_buffer.extract_params_from_vector(start);

		cout<<"Parameters are "<<" omega = "<<params_buffer.get_omega()
					<<" theta = "<<params_buffer.get_theta()
					<<" xi = "<<params_buffer.get_xi()
					<<" rho = "<<params_buffer.get_roe();
		if(model==VAR_P)
			cout<<" p = "<<params_buffer.get_p();
		cout<<endl<<endl;

		if(runmode == NORMAL_RESIDUALS) 
		{
			double* residuals = new double[prices.size()];
			for(unsigned int i = 0; i < prices.size(); i++) {
				residuals[i] = prices[i] - exp(estimates[i]);
			}

			double chi2;

			if(!is_normal(residuals, prices.size(), chi2)) 
			{
 
				if(best_params.get_mle() > mle && best_params.get_chi2() > chi2)
				{
					log_file<<"We have found better parameters on simulation # "<<simulation_counter << " out of "<< max_simulations
						<<" new mle = "<< mle <<" old mle = "<<best_params.get_mle()
						<<" new chi2 = "<<chi2<<" old chi2 = "<<best_params.get_chi2()<<endl;
					best_params.set_best_estimate(params_buffer, mle, chi2, u, v, estimates, prices.size());
				} 
				else 
				{
					log_file<<"We have NOT found better parameters on simulation # "<<simulation_counter << " out of "<< max_simulations
						<<" new mle = "<< mle <<" old mle = "<<best_params.get_mle()
						<<" new chi2 = "<<chi2<<" old chi2 = "<<best_params.get_chi2()<<endl
						<<"Resetting the parameters to the previous best estimate, and then preturbing..."<<endl;
					params_buffer.copy_params(best_params);
				}

				cout<<"Residuals are NOT normal... preturbing the starting variables and beginning again"<<endl;

				//Preturbing omega, theta and xi with N(0,2), rho with U(-0.9, 0.9) and p with U(0.5,1.5)
				current_params.set_omega(gaussrand() * sqrt(2.00) + params_buffer.get_omega());
				current_params.set_theta(gaussrand() * sqrt(2.00) + params_buffer.get_theta());
				current_params.set_xi(gaussrand() * sqrt(2.00) + params_buffer.get_xi());
				double rand_max = RAND_MAX;
				current_params.set_roe(((rand()/rand_max) * 2.00) - 1.00);
				if(model==VAR_P) {//p is alse being estimated
					current_params.set_p((((double)rand())/rand_max) + 0.5);
				}

				log_file<<"omega preturbed from "<<params_buffer.get_omega()<<" to "<<current_params.get_omega()<<endl;
				log_file<<"theta preturbed from "<<params_buffer.get_theta()<<" to "<<current_params.get_theta()<<endl;
				log_file<<"xi    preturbed from "<<params_buffer.get_xi()<<" to "<<current_params.get_xi()<<endl;
				log_file<<"rho   preturbed from "<<params_buffer.get_roe()<<" to "<<current_params.get_roe()<<endl;
				if(model==VAR_P)
					log_file<<"p changed from     "<<params_buffer.get_p()<<" to "<<current_params.get_p()<<endl;

				call_counter = 0;
			} else {
				log_file<<"Results are Normal, exiting"<<endl;
				break;
			}
		} 
		
		simulation_counter++;
		paramter_file.flush(); //flush the buffer to disc
	}
	
	residual_file<<"Prices"<<","<<"Estimates"<<","<<"Error"<<endl;
	if(runmode==SIMPLE) {
		log_file<<"Optimum values are "<<endl
			<<"--omega "<<params_buffer.get_omega()<<endl
			<<"--theta "<<params_buffer.get_theta()<<endl
			<<"--roe   "<<params_buffer.get_roe()<<endl
			<<"--xi    "<<params_buffer.get_xi()<<endl;
		if(model==VAR_P)
			log_file<<"--p     "<<params_buffer.get_p()<<endl;
		for(unsigned int i = 0; i < prices.size(); i++) {
			residual_file<<prices[i]<<","<<exp(estimates[i])<<","<<prices[i] - exp(estimates[i])<<endl;
		}
	} else {
		log_file<<"Optimum values are "<<endl
			<<"--omega "<<best_params.get_omega()<<endl
			<<"--theta "<<best_params.get_theta()<<endl
			<<"--roe   "<<best_params.get_roe()<<endl
			<<"--xi    "<<best_params.get_xi()<<endl;
		if(model==VAR_P)
			log_file<<"--p     "<<best_params.get_p()<<endl;
		for(unsigned int i = 0; i < prices.size(); i++) {
			residual_file<<prices[i]<<","<<exp(best_params.get_estimates()[i])<<","<<prices[i] - exp(best_params.get_estimates()[i])<<endl;
		}
	}

	residual_file.close();
	paramter_file.close();

	delete[] log_stock_prices, u, v, estimates;
	delete start_, identity_matrix;
	/*
	cout<<"Press any key to continue"<<endl;
	cin.get();
	*/
	return 0;
}

void parse_args(int argc, char** argv,
	string& input_file_name,
	string& parameter_file_name,
	string& residual_file_name) {

	string usage = "syntax is program_name <input_file> <parameter_output_file> <residual_output_file> "
		"<MODEL 1=HESTON, 2=GARCH, 3=3_2, 4=variable p> "
		"<RUNMODE 1=SIMPLE, 2=NORMAL_RESIDUALS, 3=UNCORRELATED_RESIDUALS, 4 = NORMAL_AND_UNCORRELATED_RESIDUALS>. "
		"<MAXSIMULATIONS: Required if RUNMODE!=1> \n "
		"The input file should be a 1 columned csv price file, the output file(s) will be created";

	//Parse the command line args.
	try {
		if(argc < 6) 
			throw usage;
		input_file_name = string(argv[1]);
		
		parameter_file_name = string(argv[2]);
		residual_file_name = string(argv[3]);
		
		string model_ = string(argv[4]);
		string runmode_ = string(argv[5]);
		if(argc==7) 
		{
			max_simulations = atoi(argv[6]);
			cout<<"max_simulations is set to "<<max_simulations<<endl;
		}

		std::cout<<"input file is "<<input_file_name<<endl;
		std::cout<<"parameter_output_file is "<<parameter_file_name<<endl;

		std::cout<<"residual_output_file is "<<residual_file_name<<endl;

		if(model_!="1" &&
			model_!="2" &&
			model_!="3" &&
			model_!="4")
			throw "model is not one of 1,2,3,4. Usage is \n" + usage;

		if(runmode_!="1" &&
			runmode_!="2" &&
			runmode_!="3" &&
			runmode_!="4")
			throw "runmode is not one of 1,2,3,4. Usage is \n" + usage;

		if(runmode_!="1" && argc!=7)
			throw "RUNMODE!=1 and MAXSIMULATIONS not provided";

		//convert them into int's
		model = (VolModel) atoi(model_.c_str()); 
		runmode = (RunMode) atoi(runmode_.c_str());
		
		ifstream ifile(input_file_name.c_str());
		if(!ifile)
			throw "could not open input file " + input_file_name;

		ifile.close();

		init_log(model, runmode, log_file);
	} catch (string& exception) {
		std::cout<<"CLI Parsing error"<<endl<<exception<<endl;
		std::cout<<"please hit a enter to continue..."<<endl;
		cin.get();
		exit(-1);
	}
}

DP minimize_target_ekf(Vec_I_DP & input) 
{
	double omega = input[0];
	double theta = input[1];
	double xi = input[2];
	double rho = input[3];
	double p;
	switch(model) 
	{
	case HESTON: p = p_heston; break;
	case GARCH: p = p_garch; break;
	case THREE_TWO: p = p_3_2; break;
	default: p = input[4];
	}

	estimate_ekf_parm_1_dim(log_stock_prices, 
		muS, 
		n_stock_prices, 
		omega, 
		theta, 
		xi, 
		rho, 
		p, //p as specified by the model
		u, 
		v, 
		estimates);

	double mle = 0.00;

	for(int i1 = 0; i1 < n_stock_prices; i1++)
		mle+=(log(v[i1])+u[i1]*u[i1]/v[i1]);

	//test for nan
	if(mle!=mle || abs(rho) > 0.95 || abs(p) > 3)
		mle = pow(10.00, 5.00) + rand(); //add a random number to it so it doesn't settle there

	log_parameters(omega, theta, xi, rho, p, mle);

	return mle;
}


void log_parameters(double& omega,
	double& theta,
	double& xi,
	double& rho,
	double& p,
	double& mle) {

	//print the header first time around
	if(call_counter==0) {
		paramter_file
			<<"simumation #"<<","
			<<"iteration"<<","
			<<"omega"<<","
			<<"theta"<<","
			<<"xi"<<","
			<<"rho"<<","
			<<"p"<<","
			<<"likelihood"<<"\n";
	}

	paramter_file
			<<simulation_counter<<","
			<<call_counter<<","
			<<omega<<","
			<<theta<<","
			<<xi<<","
			<<rho<<","
			<<p<<","
			<<mle<<"\n";

	//Printing out a status message to see whats going on
	if(call_counter++ % 100 == 0)
		cout<<" simulation counter = "<<simulation_counter
			<<" call_counter = "<<call_counter 
			<<" omega = "<<omega
			<<" theta = "<<theta
			<<" xi = "<<xi
			<<" rho = "<<rho
			<<" p = "<<p
			<<" sum = "<<mle
			<<endl;

}

void init_log(const VolModel model, const RunMode runmode, ofstream& log_stream)
{
	string logfile = "./log_" + get_model_name(model) + "_" + get_run_mode(runmode);
	cout<<"log file is: "<<logfile<<endl;
	log_stream.open(logfile.c_str());
	log_stream.setf(ios::fixed);
	log_stream<<setprecision(10);
}

string get_run_mode(const RunMode runmode) 
{
	switch(runmode)
	{
	case SIMPLE: return "simple";
	case NORMAL_RESIDUALS: return "normal_residuals";
	case UNCORRELATED_RESIDUALS: return "uncorrelated_normals";
	default: return "normal_and_uncorrelated_residuals";
	}
}

string get_model_name(const VolModel model) 
{
	switch(model)
	{
	case HESTON: return "heston";
	case GARCH: return "garch";
	case THREE_TWO: return "three_two";
	default: return "var_p";
	}
}


/*

//The function to minimize for the heston case
DP minimize_target_ekf_heston(Vec_I_DP & input);
//The function to minimize for the GARCH case
DP minimize_target_ekf_garch(Vec_I_DP & input);
//The function to minimize for the GARCH case
DP minimize_target_ekf_3_2(Vec_I_DP & input);
//The function to minimize for the variable p case
DP minimize_target_ekf_variable_p(Vec_I_DP & input);


//The function to minimize for the heston case
DP minimize_target_ekf_heston(Vec_I_DP & input) {
	double omega = input[0];
	double theta = input[1];
	double xi = input[2];
	double rho = input[3];

	estimate_ekf_parm_1_dim(log_stock_prices, 
		muS, 
		n_stock_prices, 
		omega, 
		theta, 
		xi, 
		rho, 
		p_heston, //p_heston = 0.5
		u, 
		v, 
		estimates);

	double sum = 0;
	for(int i1 = 0; i1 < n_stock_prices; i1++)
		sum+=(log(v[i1])+u[i1]*u[i1]/v[i1]);

	//test for nan
	if(sum!=sum)
		return pow(10.00, 10.00);

	log_parameters(omega, theta, xi, rho, sum);

	return sum;
}

//The function to minimize for the GARCH case
DP minimize_target_ekf_garch(Vec_I_DP & input) {
	double omega = input[0];
	double theta = input[1];
	double xi = input[2];
	double rho = input[3];

	estimate_ekf_parm_1_dim(log_stock_prices, 
		muS, 
		n_stock_prices, 
		omega, 
		theta, 
		xi, 
		rho, 
		p_garch, // p_garch = 1
		u, 
		v, 
		estimates);

	double sum = 0;
	for(int i1 = 0; i1 < n_stock_prices; i1++)
		sum+=(log(v[i1])+u[i1]*u[i1]/v[i1]);

	//test for nan
	if(sum!=sum)
		return pow(10.00, 10.00);

	log_parameters(omega, theta, xi, rho, sum);

	return sum;
}

//The function to minimize for the GARCH case
DP minimize_target_ekf_3_2(Vec_I_DP & input) {
	double omega = input[0];
	double theta = input[1];
	double xi = input[2];
	double rho = input[3];

	estimate_ekf_parm_1_dim(log_stock_prices, 
		muS, 
		n_stock_prices, 
		omega, 
		theta, 
		xi, 
		rho, 
		p_3_2, // p_3_2 = 1.50
		u, 
		v, 
		estimates);

	double sum = 0;
	for(int i1 = 0; i1 < n_stock_prices; i1++)
		sum+=(log(v[i1])+u[i1]*u[i1]/v[i1]);

	//test for nan & rho going wild
	if(sum!=sum)
		sum = pow(10.00, 10.00);

	log_parameters(omega, theta, xi, rho, sum);

	return sum;
}

//The function to minimize for the GARCH case
DP minimize_target_ekf_variable_p(Vec_I_DP & input) {
	double omega = input[0];
	double theta = input[1];
	double xi = input[2];
	double rho = input[3];
	double p = input[4];

	double sum = 0;

	//constrining rho 
	if(abs(rho) < 0.95) {
		estimate_ekf_parm_1_dim(log_stock_prices, 
			muS, 
			n_stock_prices, 
			omega, 
			theta, 
			xi, 
			rho, 
			p, //p is variable
			u, 
			v, 
			estimates);

		for(int i1 = 0; i1 < n_stock_prices; i1++)
			sum+=(log(v[i1])+u[i1]*u[i1]/v[i1]);
	} else {
		sum = pow(10.00, 6.00) * abs(rand()) * 1000.00; //if the optimizer tries unusual parameters, blow it up
	}

	//test for nan
	if(sum!=sum)
		return pow(10.00, 10.00);

	log_parameters(omega, theta, xi, rho, sum);

	return sum;
}
*/


