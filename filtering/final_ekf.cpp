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

//The function to minimize for the heston case
DP minimize_target_ekf_heston(Vec_I_DP & input);
//The function to minimize for the GARCH case
DP minimize_target_ekf_garch(Vec_I_DP & input);
//The function to minimize for the GARCH case
DP minimize_target_ekf_3_2(Vec_I_DP & input);
//The function to minimize for the variable p case
DP minimize_target_ekf_variable_p(Vec_I_DP & input);

int n_stock_prices = 0;
double *log_stock_prices, *u, *v, *estimates;
double muS = 0.0;

double p_heston = 0.50;
double p_garch  = 1.00;
double p_3_2    = 1.50;
double p_variable = 1; //This is the starting case situation where we try to optimize for variable p

int call_counter = 0;

string input_file_name, output_file_name = "", estimate_file_name = "";
ofstream output_file;
ofstream estimate_file;

DP ftol = 1.00e-14;

enum VolModel {HESTON = 1, GARCH = 2, THREE_TWO = 3, VAR_P = 4};
enum RunMode {SIMPLE = 1, NORMAL_RESIDUALS = 2, UNCORRELATED_RESIDUALS = 3, NORMAL_AND_UNCORRELATED_RESIDUALS = 4};

int main(int argc, char** argv) {
	
	srand(time(NULL)); //seed the random number generator

	vector<double> prices;
	cout<<setprecision(15);
	string usage = "syntax is program_name <input_file> <parameter_output_file> <residual_output_file> "
		"<MODEL 1=HESTON, 2=GARCH, 3=3_2, 4=variable p> "
		"<RUNMODE 1=SIMPLE, 2=NORMAL_RESIDUALS, 3=UNCORRELATED_RESIDUALS, 4 = NORMAL_AND_UNCORRELATED_RESIDUALS>. \n "
		"The input file should be a 1 columned csv price file, the output file(s) will be created";

	int model, runmode;

	//Parse the command line args.
	try {
		if(argc!=6) 
			throw usage;
		input_file_name = string(argv[1]);
		
		output_file_name = string(argv[2]);
		estimate_file_name = string(argv[3]);

		string model_ = string(argv[4]);
		string runmode_ = string(argv[5]);

		cout<<"input file is "<<input_file_name<<endl;
		cout<<"parameter_output_file is "<<output_file_name<<endl;

		cout<<"estimate_output_file is "<<estimate_file_name<<endl;

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

		//convert them into int's
		model = atoi(model_.c_str()); 
		runmode = atoi(runmode_.c_str());
		
		ifstream ifile(input_file_name.c_str());
		if(!ifile)
			throw "could not open input file " + input_file_name;

		ifile.close();
	} catch (string & exception) {
		cout<<exception<<endl;
		cout<<"please hit a key to continue..."<<endl;
		cin.get();
		exit(-1);
	}

	//setting cout to show 10 decimal places by default
	cout.setf(ios::fixed);
	cout<<setprecision(10); 

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

	bool keep_running = true;

	double omega = 0.2;
	double theta = 1.00;
	double xi = 0.50;
	double rho = -0.20;
	double p = 1.00;

	
	Vec_IO_DP* start_;
	
	if(model!=VAR_P) {
		start_ = new Vec_IO_DP(4);
	} else {
		start_ = new Vec_IO_DP(5);
	}
	Vec_IO_DP& start = *start_;
	
	start[0] = omega; start[1] = theta; start[2] = xi; start[3] = rho;
	if(model==VAR_P)
		start[4] = p;

	for(int i = 0; i < 5; i++)
		cout<<start[i]<<endl;

	while(keep_running)
	{
		output_file.open(output_file_name.c_str());
		estimate_file.open(estimate_file_name.c_str());

		//setting the output to 10 decimal places
		output_file.setf(ios::fixed);
		output_file<<setprecision(10);
		estimate_file.setf(ios::fixed);
		estimate_file<<setprecision(10);

		if(model!=VAR_P) //if we aren't trying to estimate p
		{
			//Initializing the starting point
			double identity_4_4[16]  = 
			{ 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1, 0,
			  0, 0, 0, 1};
			Mat_IO_DP identity_matrix(identity_4_4, 4, 4);
		
			int iter;
			DP fret;
			switch(model) 
			{
			case HESTON:
				NR::powell(start, identity_matrix, ftol, iter, fret, minimize_target_ekf_heston);
				break;
			case GARCH:
				NR::powell(start, identity_matrix, ftol, iter, fret, minimize_target_ekf_garch);
				break;
			case THREE_TWO:
				NR::powell(start, identity_matrix, ftol, iter, fret, minimize_target_ekf_3_2);
				break;
			}
			cout<<"Ran succesfully in "<<call_counter<<" iterations with return value "<<fret<<endl;
	
			cout<<"Parameters are "<<" omega = "<<omega
					<<" theta = "<<theta
					<<" xi = "<<xi
					<<" rho = "<<rho<<endl<<endl;

			if(!output_file_name.empty())
				output_file.close();
		} else { //we are trying to estimate p
			//Initializing the starting point
			
			double identity_5_5[25]  = 
			{ 1, 0, 0, 0, 0,
			  0, 1, 0, 0, 0,
			  0, 0, 1, 0, 0,
			  0, 0, 0, 1, 0,
			  0, 0, 0, 0, 1};
			Mat_IO_DP identity_matrix(identity_5_5, 5, 5);
	
			int iter;
			DP fret;
			NR::powell(start, identity_matrix, ftol, iter, fret, minimize_target_ekf_variable_p);
			cout<<"Ran succesfully in "<<call_counter<<" iterations with return value "<<fret<<endl;
	
			double p = start[4];

			cout<<"Parameters are "<<" omega = "<<omega
					<<" theta = "<<theta
					<<" xi = "<<xi
					<<" rho = "<<rho
					<<" p = "<<p<<endl<<endl;

			output_file.close();
		}

		omega = start[0];
		theta = start[1];
		xi = start[2];
		rho = start[3];
		if(model==VAR_P)
			p = start[4];

		if(runmode == NORMAL_RESIDUALS) 
		{
			double* residuals = new double[prices.size()];
			for(unsigned int i = 0; i < prices.size(); i++) {
				residuals[i] = prices[i] - exp(estimates[i]);
			}

			if(!is_normal(residuals, prices.size())) 
			{
				cout<<"Results are NOT normal... preturbing the starting variables and beginning again"<<endl;
				//Preturbing omega, theta and xi with N(0,1), rho with U(-0.9, 0.9) and p with U(0.5,1.5)
				start[0] = gaussrand() + omega;
				start[1] = gaussrand() + theta;
				start[2] = gaussrand() + xi;
				double rand_max = RAND_MAX;
				start[3] = (((rand()/rand_max) * 2.00) - 1.00);
				if(model==VAR_P) {//p is alse being estimated
					start[4] = ((((double)rand())/rand_max) + 0.5);
				}

				cout<<"omega "<<omega<<" "<<start[0]<<endl;
				cout<<"theta "<<theta<<" "<<start[1]<<endl;
				cout<<"xi    "<<xi<<" "<<start[2]<<endl;
				cout<<"rho   "<<rho<<" "<<start[3]<<endl;
				if(model==VAR_P)
					cout<<"p      "<<p<<" "<<start[4]<<endl;

				call_counter = 0;
					
			} else {
				cout<<"Results are Normal, exiting"<<endl;
				keep_running = false;
			}
		} else {
			keep_running = false;
		}

	}

	estimate_file<<"Prices"<<","<<"Estimates"<<","<<"Error"<<endl;
	//first line, there is no estimate or error
	estimate_file<<setprecision(10);
	
	for(unsigned int i = 0; i < prices.size(); i++) {
		estimate_file<<prices[i]<<","<<exp(estimates[i])<<","<<prices[i] - exp(estimates[i])<<endl;
	}

	estimate_file.close();

	delete[] log_stock_prices, u, v, estimates;
	delete start_;
	/*
	cout<<"Press any key to continue"<<endl;
	cin.get();
	*/
	return 0;
}

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

	//print the header first time around
	if(call_counter==0 && !output_file_name.empty()) {
		output_file<<"iteration"<<","
			<<"omega"<<","
			<<"theta"<<","
			<<"xi"<<","
			<<"rho"<<","
			<<"likelihood"<<endl;
	}

	output_file<<call_counter<<","
			<<omega<<","
			<<theta<<","
			<<xi<<","
			<<rho<<","
			<<sum<<endl;

	//Printing out a status message to see whats going on
	if(call_counter++ % 100 == 0)
		cout<<" call_counter = "<<call_counter 
			<<" omega = "<<omega
			<<" theta = "<<theta
			<<" xi = "<<xi
			<<" rho = "<<rho
			<<" sum = "<<sum
			<<endl;

	//test for nan
	if(sum!=sum)
		return pow(10.00, 10.00);

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

	//print the header first time around
	if(call_counter==0 && !output_file_name.empty()) {
		output_file<<"iteration"<<","
			<<"omega"<<","
			<<"theta"<<","
			<<"xi"<<","
			<<"rho"<<","
			<<"likelihood"<<endl;
	}

	output_file<<call_counter<<","
			<<omega<<","
			<<theta<<","
			<<xi<<","
			<<rho<<","
			<<sum<<endl;

	//Printing out a status message to see whats going on
	if(call_counter++ % 100 == 0)
		cout<<" call_counter = "<<call_counter 
			<<" omega = "<<omega
			<<" theta = "<<theta
			<<" xi = "<<xi
			<<" rho = "<<rho
			<<" sum = "<<sum
			<<endl;

	//test for nan
	if(sum!=sum)
		return pow(10.00, 10.00);

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

	//print the header first time around
	if(call_counter==0 && !output_file_name.empty()) {
		output_file<<"iteration"<<","
			<<"omega"<<","
			<<"theta"<<","
			<<"xi"<<","
			<<"rho"<<","
			<<"likelihood"<<endl;
	}

	output_file<<call_counter<<","
			<<omega<<","
			<<theta<<","
			<<xi<<","
			<<rho<<","
			<<sum<<endl;

	//test for nan
	if(sum!=sum)
		sum = pow(10.00, 10.00);

	//Printing out a status message to see whats going on
	if(call_counter++ % 100 == 0)
		cout<<" call_counter = "<<call_counter 
			<<" omega = "<<omega
			<<" theta = "<<theta
			<<" xi = "<<xi
			<<" rho = "<<rho
			<<" sum = "<<sum
			<<endl;

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

	//print the header first time around
	if(call_counter==0 && !output_file_name.empty()) {
		output_file<<"iteration"<<","
			<<"omega"<<","
			<<"theta"<<","
			<<"xi"<<","
			<<"rho"<<","
			<<"p"<<","
			<<"likelihood"<<endl;
	}

	output_file<<call_counter<<","
			<<omega<<","
			<<theta<<","
			<<xi<<","
			<<rho<<","
			<<p<<","
			<<sum<<endl;

	//Printing out a status message to see whats going on
	if(call_counter++ % 100 == 0)
		cout<<" call_counter = "<<call_counter 
			<<" omega = "<<omega
			<<" theta = "<<theta
			<<" xi = "<<xi
			<<" rho = "<<rho
			<<" p = "<<p
			<<" sum = "<<sum
			<<endl;

	//test for nan
	if(sum!=sum)
		return pow(10.00, 10.00);

	return sum;
}


/*

for(int i = 0; i < 3; i++)
		cout<<"log("<<prices[i]<<") "<<log_stock_prices[i]<<endl;


	double omega = 0.599217, 
		theta = 6.44888, 
		xi = 0.305174, 
		rho = -0.717586, 
		p = 0.5;

	estimate_ekf_parm_1_dim_heston(log_stock_prices, 
		muS, 
		n_stock_prices,
		omega,
		theta,
		xi,
		rho,
		u, 
		v, 
		estimates);

	double heston_sum = 0;
	for(int i1 = 0; i1 < n_stock_prices; i1++)
		heston_sum+=(log(v[i1])+u[i1]*u[i1]/v[i1]);

	estimate_ekf_parm_1_dim(log_stock_prices, 
		muS, 
		n_stock_prices,
		omega,
		theta,
		xi,
		rho,
		p,
		u, 
		v, 
		estimates
	

	double general_sum = 0;
	for(int i1 = 0; i1 < n_stock_prices; i1++)
		general_sum+=(log(v[i1])+u[i1]*u[i1]/v[i1]);

	cout<<"heston sum "<<heston_sum<<endl;
	cout<<"general sum "<<general_sum<<endl;
	*/
