#include <iostream>
#include <cstdlib> //to use rand function

const int N = 8; //number of agents
double T = 10.0; //average money per agent -> M = T*N is the total money
double lambda = 0.1; //constant saving propensity
const int NSIM = 100; //number of simulations
const int W = 125; //number of states    elegí W = beta*5 = 25*5

double agents[N]; //agent income (u vector)
//double credit[N]; //agent credit (v vector)
double distrib[N] = {}; //auxiliar agent vector to get the final distribution
int states[W] = {}; //occupation number per state
double random_agent();
double normal_rnum(); //epsilon
void initialize_agents(double ag[], double T); // each agent starts with the same quantity T
void initialize_agents2(double ag[], double T); //ag[0] starts with all the money
void initialize_credit(double cr[]); //each agent starts with a credit of 0
void show_dvector_row(double vec[]);
void show_agents(double ag[]);
void show_states(int st[]);
void interaction(double ag[N], int NSTEPS); //exchange rule between two random agents
void insertion_sort(double v[N]); //sorts a vector by insertion
void counting(double ag[N], int st[W]); //evaluates the occupation numbers

int main(){
	
	srand(0); //rand function seed
	//show_dvector_row(distrib);
	
	for(int k=0; k<NSIM; k++){
		initialize_agents(agents, T);
		interaction(agents, 10000);
		insertion_sort(agents);
		
		for(int l=0; l<N; l++){
			distrib[l] += agents[l];
		}
	}
	
	for(int i=0; i<N; i++){
		distrib[i] = distrib[i]/NSIM;
	}
	
	show_dvector_row(distrib);
	
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

void initialize_credit(double cr[]){
	
	for(int k=0; k<N; k++){
		cr[k] = 0.0;
	}
}

void show_dvector_row(double vec[]){
	
	for(int k=0; k<N; k++){
		std::cout<< vec[k] << " ";
	}
	std::cout<< "\n";
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
		//std::cout << i << "\t" << j << "\t" << epsilon << "\t" << Delta_m << std::endl;
	}
}

void insertion_sort(double v[N]){
	int i, pos;
	double aux;
	for(i=0; i<N; i++){
		pos = i;
		aux = v[i];
		while((pos > 0) && (v[pos-1] > aux)){
			v[pos] = v[pos-1];
			pos--;
		}
		v[pos] = aux;
	}
}

//ag[k]/T is the normalized money of agent-k
void counting(double ag[N], int st[W]){
	
	for(int k=0; k<N; k++){
		// beta is the number of states you want to define between 0 and 1 in a normalized graph
		double x = 25*ag[k]/T;  // x = beta*ag[k]/T
		
		if(x<124){
			st[(int)x] += 1;
		}
		else{
			st[124] += 1;
		}
	}
}
