#include <iostream>
#include <cstdlib>

double T = 10.0;
double agents[25];
double random_number(double nmax);
double minimum(double a, double b);
void initialize_agents(double ag[25], double T);
void show_agents(double ag[25]);
void interaction(double ag[25]);

int main(){
	
	srand(0);
	//std::cout << RAND_MAX << std::endl;
	
	initialize_agents(agents, T);
	interaction(agents);
	//show_agents(agents);
 
	return 0;
}

void initialize_agents(double ag[25], double T){
	
	for(int i=0; i<25; i++){
		ag[i] = T;
	}

}

void show_agents(double ag[25]){
	
	for(int i=0; i<25; i++){
		std::cout<< ag[i] << "\t";
	}

}

double random_number(double nmax){
	
	double rnumber = T*rand()/RAND_MAX;
	return rnumber;
}

double minimum(double a, double b){
	
	if(a<b) return a;
	else return b;
}

void interaction(double ag[25]){
	
	double Delta_m = 0;
	
	for(int i=0; i<25; i++){
		for(int j=0; j<25; j++){
			double min = minimum(agents[i], agents[j]);
			Delta_m = random_number(min);
			ag[i] += Delta_m;
			ag[j] -= Delta_m;
			//std::cout << min << "\t" << Delta_m << std::endl;
		}
	}

}

/*void interaction(double ag[25]){
	
	double Delta_m = 0;
	
	for(int i=0; i<25; i++){
		for(int j=0; j<25; j++){
			if(i=j){
			} else {
				//double min = minimum(agents[i], agents[j]);
				Delta_m = random_number(10.0);
				agents[i] += Delta_m;
				agents[j] -= Delta_m;
			}
		}
	}

}*/

