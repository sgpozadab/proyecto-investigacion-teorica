#include <iostream>
#include <cstdlib> //to use rand function

const int N = 1000; //number of agents
double T = 10.0; //average money per agent -> M = T*N is the total money
const int NSIM = 1000; //number of simulations 
const int W = 126; //number of states (beta = 25; W = beta*5 + 1)

double agents0[N]; //agent income (u vector)
double agents05[N];
double agents09[N]; 

//double credit[N]; //agent credit (v vector)
double distrib0[N] = {}; //auxiliar agent vector to get the final distribution
double distrib05[N] = {};
double distrib09[N] = {};
int states0[W] = {}; //occupation number per state
int states05[W] = {};
int states09[W] = {};

float xaxis[W]; //x axis

double random_agent();
double normal_rnum(); //epsilon
void initialize_agents(double ag[], double T); // each agent starts with the same quantity T
void initialize_agents2(double ag[], double T); //ag[0] starts with all the money
void initialize_credit(double cr[]); //each agent starts with a credit of 0
void initialize_xaxis(float xa[], int beta); //x axis
void show_dvector_row(double vec[]);
void show_agents(double ag[]);
void show_states(int st[]);
void interaction(double ag[N], float lambda, int NSTEPS); //exchange rule between two random agents
void insertion_sort(double v[N]); //sorts a vector by insertion
void counting(double ag[N], int st[W]); //evaluates the occupation numbers

int main(){
	
	srand(0); //rand function seed
	initialize_xaxis(xaxis, 25);
	
	for(int k=0; k<NSIM; k++){
		initialize_agents(agents0, T);
		interaction(agents0, 0.0, 100000);
		insertion_sort(agents0);
		
		for(int l=0; l<N; l++){
			distrib0[l] += agents0[l];
		}
	}
	for(int i=0; i<N; i++){
		distrib0[i] = distrib0[i]/NSIM;
	}
	counting(distrib0, states0);

	for(int k=0; k<NSIM; k++){
		initialize_agents(agents05, T);
		interaction(agents05, 0.5, 100000);
		insertion_sort(agents05);
		
		for(int l=0; l<N; l++){
			distrib05[l] += agents05[l];
		}
	}
	for(int i=0; i<N; i++){
		distrib05[i] = distrib05[i]/NSIM;
	}
	counting(distrib05, states05);

	for(int k=0; k<NSIM; k++){
		initialize_agents(agents09, T);
		interaction(agents09, 0.9, 100000);
		insertion_sort(agents09);
		
		for(int l=0; l<N; l++){
			distrib09[l] += agents09[l];
		}
	}
	for(int i=0; i<N; i++){
		distrib09[i] = distrib09[i]/NSIM;
	}
	counting(distrib09, states09);
	
	for(int k=0; k<W-1; k++){
	  std::cout<< xaxis[k] << "\t" << states0[k] << "\t" << states05[k] << "\t" << states09[k] <<std::endl;
	}
	
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

void initialize_xaxis(float xa[], int beta){
	for(int k=0; k<W-1; k++){
	  xa[k] = 1.0*(k+1)/beta;
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
	
	for(int k=0; k<W-1; k++){
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

void interaction(double ag[N], float lambda, int NSTEPS){
	
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
		
		if(x<125){
			st[(int)x] += 1;
		}
		else{
		        st[125] += 1;
		}
	}
}
