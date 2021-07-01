#include<iostream>
#include<cmath>
#include<vector>
#include<random>
#include<cstdlib>
/*
//Funciones que definen el modelo de CC, si lambda0=0 obtengo el modelo de DY
double Xi(double lambda0,double xi0,double xj0,double epsilon0){
    return lambda0*xi0+epsilon0*(1.-lambda0)*(xi0+xj0);
}

double Xj(double lambda0,double xi0,double xj0,double epsilon0){
    return lambda0*xj0+(1.-epsilon0)*(1.-lambda0)*(xi0+xj0);
}
*/


//Funcion que define el modelo de BM
double Xi(double J0,double xi0,double xj0,double eta0,double etamean0, double sigma0, double Wnormalized0){
    return xi0+((eta0-etamean0-sigma0*sigma0)*Wnormalized0+J0*(1-Wnormalized0));
}

void init_x(std::vector<double> & y){
for(int i=0;i<y.size();i++){
y[i]=1;
}

}
void print(std::vector<double> & y){
for(int i=0;i<y.size();i++){
std::cout<<y[i]<<std::endl;
}

}

int main(int argc, char * argv[]){
int N=atoi(argv[1]);
double epsilon=1;
double eta=1;
double etamean=0;
double sigma=0.5;
double Wmean=1;
double Winormalized=1;
int Ni=1,Nj=2;
double tempi=0.1,tempj=0.2;
std::vector<double> x;
x.resize(N);
init_x(x);
int it=atoi(argv[2]);
std::vector<double> lambda1 {0.95,0.9};
int runs=atoi(argv[3]);
//double landa=atof(argv[4]);
double J=atof(argv[4]);
srand((int) time(0));

//Clase que me genera numeros aleatorios de una distribucion normal para la variable estocastica eta
class Generator {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution;
    double min;
    double max;
public:
    Generator(double mean, double stddev, double min, double max):
        distribution(mean, stddev), min(min), max(max)
    {}

    double operator ()() {
        while (true) {
            double number = this->distribution(generator);
            if (number >= this->min && number <= this->max)
                return number;
        }
    }
};

//Proceso en donde alterno entre el modelo de CC y de BM comentando y descomentando algunas lineas
 
//for(int k=0;k<lambda1.size();k++){ 
	for(int i=0;i<runs;i++){
		init_x(x);
                Generator g(etamean,sigma,-3,3);
                eta=g();
		for(int n=0;n<it;n++){
			Ni=(rand() % N);     
			Nj=(rand() % N);
                        //Wmean=accumulate( x.begin(), x.end(), 0.0)/x.size();     
			//epsilon=(double)rand() / (double)RAND_MAX;
			//std::cout<<epsilon<<std::endl;
			tempi=x[Ni];
			//tempj=x[Nj];
                        Winormalized=tempi/Wmean;
			//x[Ni]=Xi(landa,tempi,tempj,epsilon);
			//x[Nj]=Xj(landa,tempi,tempj,epsilon);
			x[Ni]=Xi(J,tempi,tempj,eta,etamean,sigma,Winormalized);
                        while(x[Ni]<=0){
                        x[Ni]=Xi(J,tempi,tempj,eta,etamean,sigma,Winormalized);
			//x[Nj]=Xj(J,tempi,tempj,eta,etamean,sigma,Wjnormalized);
                        }
                        Wmean=Wmean+(tempi-x[Ni]/N);
		}
		print(x);
	}
//}
return 0;
}
