#include <iostream>
#include <cstdlib>

double lambda = 0.5; //constant saving propensity
const int N = 250; //number of agents
const int W = 13; //number of states
double T = 10.0; //average money per agent
//M = T*N is the total money 

double agents[N]; //money per agent
int states[W] = {}; //occupation number per state
double random_agent();
double normal_rnum(); //epsilon
void initialize_agents(double ag[], double T); // each agent starts with the same quantity T
void initialize_agents2(double ag[], double T); //ag[0] starts with all the money
void show_agents(double ag[]);
void show_states(int st[]);
void interaction(double ag[N], int NSTEPS); //exchange rule between two random agents
void counting(double ag[N], int st[W]); //evaluates the occupation numbers

int main(){
	
	srand(0); //rand function seed
	
	/*for(int k=0; k<20; k++){
		std::cout << random_agent() << std::endl;
	}*/
	
	initialize_agents(agents, T);
	interaction(agents, 1000000);
	counting(agents, states);
	show_states(states);
 
	return 0;
}

void initialize_agents(double ag[], double T){
	
	for(int k=0; k<N; k++){
		ag[k] = T;
	}
}

void initialize_agents2(double ag[], double T){
	ag[0] = T*N;  
	
	for(int k=1; k<N; k++){
		ag[k] = 0;
	}
}

void show_agents(double ag[]){
	
	for(int k=0; k<N; k++){
		std::cout<< ag[k] << "\n";
	}
}

void show_states(int st[]){
	
	for(int k=0; k<W; k++){
		std::cout<< st[k] << "\n";
	}
}

//rand function generates random numbers between 0 and RAND_MAX=32767
double random_agent(){
	
	int ragent = rand()%N;
	//generates integer numbers between 0 and N-1
	return ragent;
}

//rand()/RAND_MAX generates numbers between 0 and 1
double normal_rnum(){
	
	double repsilon = 1.0*rand()/RAND_MAX;
	return repsilon;
}

void interaction(double ag[N], int NSTEPS){
	
	for(int k=0; k<NSTEPS; k++){
		double Delta_m = 0;
		double epsilon = normal_rnum();
		int i = random_agent();
		int j = random_agent();
		
		if(i!=j){
			Delta_m = (1.0-lambda)*(epsilon*ag[j]-(1.0-epsilon)*ag[i]);
			ag[i] += Delta_m;
			ag[j] -= Delta_m;
		}
		//std::cout << a << "\t" << b << "\t" << epsilon << "\t" << Delta_m << std::endl;
	}
}

//ag[k]/T is the normalized money of agent-k
void counting(double ag[N], int st[W]){
	
	for(int k=0; k<N; k++){
		double x = 3.0*ag[k]/T;
		
		if(x<12){
			st[(int)x] += 1;
		}
		else{
			st[12] += 1;
		}
	}
}
