#ifndef ReplicatorNKT_TestCases_TestCase2MatrixDiffusion_TestCase2MatrixDiffusion
#define ReplicatorNKT_TestCases_TestCase2MatrixDiffusion_TestCase2MatrixDiffusion

#include "TestCase.h"
#include "ReplicatorModel.h"

class TestCase2MatrixDiffusion : public TestCase {
public:
	void RunTest() {		
		const double maxTime = 10.0; //Maximal simulation time
		const double snapshotTime = 0.1; //Save solution every snapshot time
		const double xMin = 0.0; //Left border coordinate
		const double xMax = 1.0; //Right border coordinate		
		const int n = 2; //Number of strategies for first species
		const int m = 2; //Number of strategies for second species
		const int nCells = 100; //Number of cells

		//Specify game matrices

		//Payoff matrix 
		DenseMatrix A(n, m);
		A[0][0] = 0.8;
		A[0][1] = 1.1;
		A[1][0] = 1.2;
		A[1][1] = 0.9;
		
		DenseMatrix B(m, n);
		A[0][0] = 0.8;
		A[0][1] = 1.1;
		A[1][0] = 1.2;
		A[1][1] = 0.9;


		//For replicator model class
		std::vector<int> nStrategies; 
		nStrategies.push_back(n);
		nStrategies.push_back(m);
		ReplicatorModel kernel(2, nStrategies, nCells);

		kernel.SetPayoffMatrix(0, 1, A);	
		kernel.SetPayoffMatrix(1, 0, B);

		//Diffusion
		kernel.SetDiffusionCoefficient(0, 0, 0.03); //0.03
		kernel.SetDiffusionCoefficient(0, 1, 0.02); //0.02
		kernel.SetDiffusionCoefficient(1, 0, 0.03); //0.03
		kernel.SetDiffusionCoefficient(1, 1, 0.02); //0.02
		//kernel.SetDiffusionCoefficient(0, 0, 0); //0.03
		//kernel.SetDiffusionCoefficient(0, 1, 0); //0.02

		//Initial values		
		//Fill initial conditions					
		double u0avg = 0.55;
		double u1avg = 0.45;
		std::vector<std::vector<double> > initValues;
		for (int i = 0; i< nCells; i++) {
			std::vector<double> u;
			const double du0 = -0.5;
			double u0 = u0avg + 2 * du0 * (kernel.getVertex(i).x - 0.5);
			double u1 = u1avg;
			//Delta
			/*if (i == 0) {
				u0 = u0avg / nCells;
				u1 = u1avg / nCells; 
			} else {
				u0 = 0;
				u1 = 0;
			};*/
			u.push_back(u0);
			u.push_back(u1);			
			initValues.push_back(u);	
		};

		//Temporary for 2 cell test
	/*	initValues.clear();
		std::vector<double> us;
		us.push_back(u0avg * 2.0 / 3.0);
		us.push_back(u1avg * 2.0 / 3.0);
		initValues.push_back(us);
		us.clear();
		us.push_back(u0avg * 1.0 / 3.0);
		us.push_back(u1avg * 1.0 / 3.0);
		initValues.push_back(us);*/

		kernel.SetInitialConditions(initValues);

		//Integration method
		kernel.SetMethod(ExplicitEuler);		

		//Run simulation
		kernel.RunCalculation(100000, 0.0, snapshotTime, maxTime);
	};
};

#endif