#include "ReplicatorModel.h"
#include "Tests.h"

int main(int argc, char *argv[]) {
	TestCaseBratus test;
	test.RunTest();
	return 0;

	double dt = 0.01;
	double maxTime = 20.0;
	const double xMin = 0.0;
	const double xMax = 1.0;	
	const int nX = 1;
	const int n = 2;
	const int m = 2;


	//Specify game matrices
	double E = 0.1;	//
	double G = 1.0; //
	double C = 1.1; //
	//Check setting
	//0 < E < G < C < 2*(G-E)
	bool isGoodSetting;
	if ((0 < E) && (E<G) && (G<C) && (C < 2*(G-E))) {
		isGoodSetting = true;
	} else {
		isGoodSetting = false;
	};
	if (!isGoodSetting) {
		std::cout<<"Setting payoffs is wrong";
		return 0;
	};
	//Payoff matrix for male species
	DenseMatrix A(n, m);
	A[0][0] = 0.0;
	A[0][1] = G;
	A[1][0] = G - C/2 - E;
	A[1][1] = G - C/2;
	//Payoff matrix for female species
	DenseMatrix B(n, m);
	B[0][0] = 0.0;
	B[0][1] = G - C/2 - E;
	B[1][0] = G - C;
	B[1][1] = G - C/2;

	//Fill initial conditions	
	std::vector<double> u;
	u.push_back(0.4);
	u.push_back(1.0 - u[0]);
	std::vector<double> v;
	v.push_back(0.6);
	v.push_back(1.0 - v[0]);	

	//For replicator model class
	std::vector<int> nStrategies;
	nStrategies.push_back(2);
	nStrategies.push_back(2);
	ReplicatorModel kernel(2, nStrategies, 100);

	kernel.SetPayoffMatrix(0, 1, A);
	kernel.SetPayoffMatrix(1, 0, B);

	std::vector<std::vector<double> > initValues;
	initValues.push_back(u);
	initValues.push_back(v);
	kernel.SetInitialConditions(initValues);
	kernel.SetMethod(ExplicitEuler);
	//kernel.SetDiffusionCoefficient(1.0);

	kernel.RunCalculation(100000, 0.0, 0.1, maxTime, 0.1, 0.1);
	return 0;
};
