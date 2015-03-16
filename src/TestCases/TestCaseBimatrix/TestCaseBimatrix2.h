#ifndef ReplicatorNKT_TestCases_TestCase2MatrixDiffusion_TestCaseBimatrix2
#define ReplicatorNKT_TestCases_TestCase2MatrixDiffusion_TestCaseBimatrix2

#include "TestCase.h"
#include "ReplicatorModel.h"

class TestCaseBimatrix2 : public TestCase {
public:
	void RunTest() {		
		const double maxTime = 100.0; //Maximal simulation time
		const double snapshotTime = 0.1; //Save solution every snapshot time
		const double xMin = 0.0; //Left border coordinate
		const double xMax = 1.0; //Right border coordinate	
		const double CFL = 0.1; //CFL number
		const double maxTimeStep = 0.01; //
		const int n = 2; //Number of strategies for first species
		const int m = 2; //Number of strategies for second species
		const int nCells = 100; //Number of cells

		//Specify game matrices
		double a = 3;

		double b = a + 2;
		double d = a + 1;

		double c = a + 2;
		double e = a + 4;

		//Payoff matrix A Case 2
		DenseMatrix A(n, m);
		A[0][0] = a;
		A[0][1] = b;
		A[1][0] = c;
		A[1][1] = 0;

		//Payoff matrix B Case 2
		DenseMatrix B(n, m);
		B[0][0] = a;
		B[0][1] = d;
		B[1][0] = e;
		B[1][1] = 0;

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
		const double D = 1;
		kernel.SetDiffusionCoefficient(0, 0, 0.02 * D); //0.03
		kernel.SetDiffusionCoefficient(0, 1, 0.02 * D); //0.02
		kernel.SetDiffusionCoefficient(1, 0, 0.02 * D); //0.03
		kernel.SetDiffusionCoefficient(1, 1, 0.02 * D); //0.02	

		//Output phase portrait with diffusion
		//kernel.WritePhasePortrait("gameDiffusion.dat", 10);

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

		auto dist = std::bind(normalDistribution, std::placeholders::_1, 0.4, 0.1);
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