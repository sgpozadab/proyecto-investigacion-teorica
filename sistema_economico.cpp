#include <iostream>
#include <cstdlib>

double T = 10.0; //average money per agent
double lambda = 0.5; // constant saving propensity
const int N = 25; //number of agents

double agents[N];
double random_number(double nmax);
double random_agent();
double normal_rnum();
double minimum(double a, double b);
void initialize_agents(double ag[N], double T);
void show_agents(double ag[N]);
void interaction(double ag[N], int NSTEPS);

int main(){
	
	srand(45); //rand function seed
	
	/*for(int k=0; k<20; k++){
		std::cout << random_agent() << std::endl;
	}*/
	
	initialize_agents(agents, T);
	interaction(agents, 1000000);
	show_agents(agents);
 
	return 0;
}

void initialize_agents(double ag[N], double T){
	
	for(int i=0; i<N; i++){
		ag[i] = T;
	}
}

void show_agents(double ag[N]){
	
	for(int i=0; i<N; i++){
		std::cout<< ag[i] << "\n";
	}
}

//rand function generates random numbers between 0 and RAND_MAX=32767
//rand()/RAND_MAX generates numbers between 0 and 1
double random_number(double nmax){
	
	double rnumber = 2*nmax*rand()/RAND_MAX-nmax;
	//generates numbers between -nmax and +nmax
	return rnumber;
}

double random_agent(){
	
	int ragent = rand()%N;
	//generates integer numbers between 0 and N-1
	return ragent;
}

double normal_rnum(){
	
	double repsilon = 1.0*rand()/RAND_MAX;
	return repsilon;
}

double minimum(double a, double b){
	
	if(a<b) return a;
	else return b;
}

/*void interaction(double ag[N]){
	
	double Delta_m = 0;
	
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			double min = minimum(agents[i], agents[j]);
			Delta_m = random_number(min);
			ag[i] += Delta_m;
			ag[j] -= Delta_m;
			//std::cout << min << "\t" << Delta_m << std::endl;
		}
	}
}*/

void interaction(double ag[N], int NSTEPS){
	
	for(int k=0; k<NSTEPS; k++){
		double Delta_m = 0;
		double epsilon = normal_rnum();
		int a = random_agent();
		int b = random_agent();
		
		if(a!=b){
			Delta_m = (1.0-lambda)*(epsilon*ag[b]-(1.0-epsilon)*ag[a]);
			ag[a] += Delta_m;
			ag[b] -= Delta_m;
		}
		//std::cout << a << "\t" << b << "\t" << epsilon << "\t" << Delta_m << std::endl;
	}
}

