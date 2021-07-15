#include <iostream>
#include <cstdlib> //to use rand()
#include <cmath> //to use abs()

//---------------------------------------------CONSTANTS AND VECTORS
const int N = 9; //number of agents
double T = 10.0; //average money per agent -> M = T*N is the total money
double d = -T; //maximum debt allowed
const int NSIM = 100; //number of simulations 
const int W = 126; //number of states (beta = 25; W = beta*5 + 1)

double agents0[N]; //agent income (u vector)
double cr[N]; //credit vector
double agents05[N];
double agents09[N]; 

double distrib0[N] = {}; //auxiliar agent vector to get the final distribution
double distrib05[N] = {};
double distrib09[N] = {};
int states0[W] = {}; //occupation number per state
int states05[W] = {};
int states09[W] = {};

float xaxis[W]; //x axis
double Y[N] = {}; //Y variable of the Gini coefficient

//----------------------------------------------FUNCTION DECLARATION
double random_agent();
double normal_rnum(); //epsilon
void initialize_agents(double ag[], double T); // each agent starts with the same quantity T
void initialize_agents2(double ag[], double T); //ag[0] starts with all the money
void initialize_agents3(double ag[], double T); //warning: use this function only if N=9
void initialize_credit(double cr[]); //each agent starts with a credit of 0
void initialize_xaxis(float xa[], int beta); //x axis
void show_dvector_row(double vec[]);
void show_agents(double ag[]);
void show_states(int st[]);
void interaction(double ag[N], float lambda, int NSTEPS); //exchange rule between two random agents
void interaction_cr(double ag[N], double cr[N], float lambda, int NSTEPS);
void insertion_sort(double v[N]); //sorts a vector by insertion
void counting(double ag[N], int st[W]); //evaluates the occupation numbers
double Gini_coef(double vec[]);

//--------------------------------------------------MAIN FUNCTION
int main(){
	
	srand(0); //rand function seed
	//initialize_xaxis(xaxis, 25);
	
	initialize_agents(agents0, T);
	//initialize_credit(cr);
	show_dvector_row(agents0);
	interaction(agents0, 0.9, 100000);
	//interaction_cr(agents0, cr, 0.5, 20);
	show_dvector_row(agents0);
	//show_dvector_row(cr);
	
	double a = 0;
	for(int i=0; i<N; i++){
		a += agents0[i];
	}
	
	std::cout<<a<<"\n";
	
/*	double n1=3;
	double n2=-1;
	if(-1<n1<1 && n2<0){
		std::cout<<"raios"<<"\n";
	}
	if(-1<n1<4 && n2<0){
		std::cout<<"genial"<<"\n";
	}*/
	
	/*for(int k=0; k<NSIM; k++){
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
	//counting(distrib0, states0);*/

/*	for(int k=0; k<NSIM; k++){
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
	//counting(distrib05, states05);

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
	//counting(distrib09, states09);
	
	double Gini0 = Gini_coef(distrib0);      //0,4994
	double Gini05 = Gini_coef(distrib05);    //0,2738
	double Gini09 = Gini_coef(distrib09);    //0,1060
	
	std::cout<< Gini0 << "\n"
	         << Gini05 << "\n"
			 << Gini09 << "\n"; */
	
/*	for(int k=0; k<W-1; k++){
	  std::cout<< xaxis[k] << "\t" << states0[k] << "\t" << states05[k] << "\t" << states09[k] <<std::endl;
	}*/
	
	return 0;
}

//----------------------------------------------------FUNCTIONS
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

void initialize_agents3(double ag[], double T){
	for(int k=0; k<5; k++){
		ag[k] = 0;
	}
	for(int k=5; k<N; k++){
		ag[k] = N*T/4;
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
	
	double repsilon = 1.0 - (2.0*rand()/RAND_MAX);
	return repsilon;
}

void interaction(double ag[N], float lambda, int NSTEPS){
	
	for(int k=0; k<NSTEPS; k++){
		double Delta_m = 0;
		double epsilon = normal_rnum()/2.0; //with lambda=0 the money grows fast
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

void interaction_cr(double ag[N], double cr[N], float lambda, int NSTEPS){
	
	for(int k=0; k<NSTEPS; k++){
		double Delta_m = 0, fi = 0, fj = 0, absD_m = 0.0;
		double epsilon = normal_rnum()/1.0; //with lambda=0 the money grows fast
		int i = random_agent();
		int j = random_agent();
		
		if(i!=j){
			Delta_m = (1.0-lambda)*(epsilon*ag[j]-(1.0-epsilon)*ag[i]);
			fi = ag[i] + Delta_m;
			fj = ag[j] - Delta_m;
			absD_m = fabs(Delta_m);
			
			if(fi>0 && fj>0){
				ag[i] = fi;
				ag[j] = fj;
				//std::cout<<"Holi"<<"\n";
			}
			
			if(d<fi<0 && fj>fabs(fi)){
				cr[i] -= absD_m;
				cr[j] += absD_m;
				std::cout<<"Holi"<<"\n";
			}
			
			if(d<fj<0 && fi>fabs(fj)){
				cr[i] += absD_m;
				cr[j] -= absD_m;
				std::cout<<"Holi"<<"\n";
			}
		}
		std::cout << i <<"\t"<< j <<"\t"<< absD_m <<"\t"<< fi <<"\t"<< fj <<std::endl;
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

double Gini_coef(double vec[N]){
	double y = 0;
	
	for(int k=0; k<N; k++){
		y += vec[k];
		Y[k] = y;
	}
	
	double A = Y[0];
	
	for(int k=2; k<=N; k++){
		A += (Y[k-1] + Y[k-2]);
	}
	
	double G = 1 - A/(N*Y[N-1]);
	return G;
}

/* coeficiente de Gini
¿quiénes son X_i?... pues X_0=0, X_N=1 ---> X_i=i/N para i=0,1,...,N
para Y_i... Y_0=0, Y_N=1=m_max/m_max ---> Y_i=Y[i-1]/Y[N-1] para i=1,2,...,N; Y_0=0

luego (X_i-X_{i-1}) = i/N-(i-1)/N = (i-i+1)/N = 1/N
	  (Y_i+Y_{i-1}) = (Y[i-1]+Y[i-2])/Y[N-1] para i=2,3,...,N. (Y_i-Y_{i-1}) = Y[0]/Y[N-1] para i=1
	  (X_i-X_{i-1})(Y_i+Y_{i-1}) = (1/N)(Y[i-1]+Y[i-2])/Y[N-1] para i=2,3,...,N
	  
G = 1 - \sum_{1}^{N} (X_i-X_{i-1})(Y_i+Y_{i-1})
  = 1 - (Y[0] + \sum_{2}^{N} (Y[i-1]+Y[i-2]))/(N*Y[N-1])
  = 1 - A/(N*Y[N-1])
  
A = Y[0] + \sum_{2}^{N} (Y[i-1]+Y[i-2])
*/
