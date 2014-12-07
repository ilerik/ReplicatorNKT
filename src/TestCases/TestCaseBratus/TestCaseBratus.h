#ifndef ReplicatorNKT_TestCases_TestCaseBratus_TestCaseBratus
#define ReplicatorNKT_TestCases_TestCaseBratus_TestCaseBratus

#include "TestCase.h"
#include "ReplicatorModel.h"

class TestCaseBratus : public TestCase {
public:
	void PrintTestInfo() {
		std::string info = "Example. 1 from paper by Bratus et al.";
		std::cout<<info;
		std::cout<<"\n";
	};

	void RunTest() {		
		const double maxTime = 10.0; //Maximal simulation time
		const double snapshotTime = 0.1; //Save solution every snapshot time
		const double xMin = 0.0; //Left border coordinate
		const double xMax = 1.0; //Right border coordinate		
		const double CFL = 0.1; //CFL number
		const double maxTimeStep = 0.01; //
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


		//For replicator model class
		std::vector<int> nStrategies; 
		nStrategies.push_back(n);
		ReplicatorModel kernel(1, nStrategies, nCells);

		kernel.SetPayoffMatrix(0, 0, A);	

		//Diffusion
		kernel.SetDiffusionCoefficient(0, 0, 0.03); //0.03
		kernel.SetDiffusionCoefficient(0, 1, 0.02); //0.02

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

		kernel.SetInitialConditions(initValues);

		//Integration method
		kernel.SetMethod(ExplicitEuler);		

		//Run simulation
		kernel.RunCalculation(100000, 0.0, snapshotTime, maxTime, CFL, maxTimeStep);
	};
};

#endif