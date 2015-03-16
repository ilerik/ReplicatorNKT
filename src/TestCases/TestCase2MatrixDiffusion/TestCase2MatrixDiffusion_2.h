#ifndef ReplicatorNKT_TestCases_TestCase2MatrixDiffusion_TestCase2MatrixDiffusion_2
#define ReplicatorNKT_TestCases_TestCase2MatrixDiffusion_TestCase2MatrixDiffusion_2

#include "TestCase.h"
#include "ReplicatorModel.h"

class TestCase2MatrixDiffusion_2 : public TestCase {
public:
	void RunTest() {		
		const double maxTime = 2.0; //Maximal simulation time
		const double snapshotTime = 0.01; //Save solution every snapshot time
		const double xMin = 0.0; //Left border coordinate
		const double xMax = 1.0; //Right border coordinate	
		const double CFL = 0.1; //CFL number
		const double maxTimeStep = 0.01; //
		const int n = 2; //Number of strategies for first species
		const int m = 2; //Number of strategies for second species
		const int nCells = 1; //Number of cells

		//Case 2
		//Payoff matrix 
		DenseMatrix A(n, m);
		A[0][0] = 5;
		A[0][1] = 1;
		A[1][0] = 4;
		A[1][1] = 3;


		//Payoff matrix 
		DenseMatrix B(n, m);
		B[0][0] = 4;
		B[0][1] = 6;
		B[1][0] = 0;
		B[1][1] = -2;


		//For replicator model class
		std::vector<int> nStrategies; 
		nStrategies.push_back(n);
		nStrategies.push_back(m);
		ReplicatorModel kernel(2, nStrategies, nCells);

		kernel.SetPayoffMatrix(0, 1, A);	
		kernel.SetPayoffMatrix(1, 0, B);

		//Diffusion
		const double D = 0.3;
		kernel.SetDiffusionCoefficient(0, 0, 0.02 * D); //0.03
		kernel.SetDiffusionCoefficient(0, 1, 0.02 * D); //0.02
		kernel.SetDiffusionCoefficient(1, 0, 0.02 * D); //0.03
		kernel.SetDiffusionCoefficient(1, 1, 0.02 * D); //0.02

		//Initial values		
		//Fill initial conditions					
		double u0avg = 0.55;
		double u1avg = 0.45;
		double v0avg = 0.6;
		double v1avg = 0.4;
		std::vector<std::vector<double> > initValues;
		for (int i = 0; i< nCells; i++) {
			std::vector<double> u;
			const double du0 = 0.5;
			double u0 = u0avg + 2 * du0 * (kernel.getVertex(i).x - 0.5);
			double u1 = u1avg;
			double v0 = v0avg;
			double v1 = v1avg;
			if (( i == nCells / 2) || (i == (nCells-1.0) / 2)) {
				u0 = 0.5 * nCells * u0avg;
			};
			u.push_back(u0);
			u.push_back(u1);
			u.push_back(v0);
			u.push_back(v1);			
			initValues.push_back(u);	
		};

		kernel.SetInitialConditions(initValues);

		//Integration method
		kernel.SetMethod(ExplicitEuler);		

		//Run simulation
		kernel.RunCalculation(10000000, 0.0, snapshotTime, maxTime, CFL, maxTimeStep);
	};
};

#endif