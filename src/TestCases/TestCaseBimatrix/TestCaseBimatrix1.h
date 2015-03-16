#ifndef ReplicatorNKT_TestCases_TestCase2MatrixDiffusion_TestCaseBimatrix1
#define ReplicatorNKT_TestCases_TestCase2MatrixDiffusion_TestCaseBimatrix1

#include "TestCase.h"
#include "ReplicatorModel.h"

class TestCaseBimatrix1 : public TestCase {
public:
	void RunTest() {		
		const double maxTime = 20.0; //Maximal simulation time
		const double snapshotTime = 0.1; //Save solution every snapshot time
		const double xMin = 0.0; //Left border coordinate
		const double xMax = 1.0; //Right border coordinate	
		const double CFL = 0.1; //CFL number
		const double maxTimeStep = 0.01; //
		const int n = 2; //Number of strategies for first species
		const int m = 2; //Number of strategies for second species
		const int nCells = 100; //Number of cells

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
			return;
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

		//For replicator model class
		std::vector<int> nStrategies; 
		nStrategies.push_back(n);
		nStrategies.push_back(m);
		ReplicatorModel kernel(2, nStrategies, nCells);

		kernel.SetPayoffMatrix(0, 1, A);	
		kernel.SetPayoffMatrix(1, 0, B);

		//Integration method
		kernel.SetMethod(ExplicitEuler);	

		//Output phase portrait without diffusion
		kernel.WritePhasePortrait("game.dat", 100);

		//Diffusion
		const double D = 1; //On
		kernel.SetDiffusionCoefficient(0, 0, 0.02 * D); //0.03
		kernel.SetDiffusionCoefficient(0, 1, 0.02 * D); //0.02
		kernel.SetDiffusionCoefficient(1, 0, 0.02 * D); //0.03
		kernel.SetDiffusionCoefficient(1, 1, 0.02 * D); //0.02	

		//Fill initial conditions					 
		double u0avg = 0.50;
		double u1avg = 0.50;
		double v0avg = 0.65;
		double v1avg = 0.35;

		auto normalDistribution = [&](double x, double my, double sd) {
			const double pi = 3.14159265359;
			double pd = (1.0 / sd) * (1.0 / std::sqrt(2.0 * pi)) * std::exp( -0.5 * std::pow(x - my, 2) / std::pow(sd, 2) );
			return pd;
		};

		auto dist = std::bind(normalDistribution, std::placeholders::_1, 0.4, 0.5);
		//auto dist = [](double x) { return 1.0; };

		std::vector<double> avg;
		avg.push_back(u0avg);
		avg.push_back(u1avg);
		avg.push_back(v0avg);
		avg.push_back(v1avg);

		kernel.SetInitialConditions(avg, dist);

		//Run simulation
		kernel.RunCalculation(10000000, 0.0, snapshotTime, maxTime, CFL, maxTimeStep);
	};
};

#endif