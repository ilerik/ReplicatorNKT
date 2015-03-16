#ifndef ReplicatorNKT_ReplicatorModel_ReplicatorModel
#define ReplicatorNKT_ReplicatorModel_ReplicatorModel

#include <fstream>
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <cmath>
#include <sstream>
#include <cassert>
#include <functional>
#include "bicgstab.h"

double computeFitness() {
	return 0;
};

double computeAverageFitness() {
	return 0;
};

class DenseMatrix {
	typedef std::vector<double> DenseMatrixRow;
	std::vector<DenseMatrixRow > _a;
	int _n;
	int _m;
public:
	//Accessor
	DenseMatrixRow& operator[](int i) {
		return _a[i];
	};

	DenseMatrix() {		
	};

	//Constructor
	DenseMatrix(int n, int m) {
		_a.resize(n);
		_n = n;
		_m = m;
		for (int i = 0; i<n; i++) _a[i].resize(m, 0);
	};

	//Multiply by vector
	std::vector<double> operator*(const std::vector<double>& x) {
		std::vector<double> res(_n);
		for (int i = 0; i<_n; i++) {
			res[i] = 0;
			for (int j = 0; j<_m; j++) res[i] += _a[i][j] * x[j];
		};
		return res;
	};

	//Multiply by vector
	std::vector<double> operator*(const double* x) {
		std::vector<double> res(_n);
		for (int i = 0; i<_n; i++) {
			res[i] = 0;
			for (int j = 0; j<_m; j++) res[i] += _a[i][j] * x[j];
		};
		return res;
	};
};

double operator*(const std::vector<double>&a, const std::vector<double>& b) {
	int size = a.size();
	if (size != b.size()) throw 0;
	double res = 0;
	for (int i = 0; i<size; i++) {
		res += a[i] * b[i];		
	};
	return res;
};

double mult(const double *a, const double* b, int size) {	
	double res = 0;
	for (int i = 0; i<size; i++) {
		res += a[i] * b[i];		
	};
	return res;
};

struct Vertex {
public:
	double x;	//coordinate
	double h;	//size
	std::vector<int> neighbours; //List of adjacent vertices indexes
};

enum IntegrationMethod {
	ExplicitEuler,
	ExplicitExponentialEuler
};

class ReplicatorModel {	
	//model parameters
	int _nSpecies;
	std::vector<int> _nStrategies;
	int _nVariables;
	std::vector<int> _shift;
	int _nCells;

	//internal storage
	std::vector<double> _values; //values
	std::vector<double> _residual; //residuals

	//interaction data (payoff matrices)
	std::map<int, std::map<int, DenseMatrix> > _games;

	//diffusion coefficients 
	std::map<int, std::map<int, double>> _d;

	//grid topology
	//for now simple
	std::vector<Vertex> _vertices; //list of vertices		

	//indexing accessors
	inline int _getIndex(int vertexInd, int specieInd, int strategyInd) {
		int ind = vertexInd * _nVariables + _shift[specieInd] + strategyInd;
		return ind;
	};	

	//Calculation settings
	double _CFL;
	double _maxTime;
	int _maxIter;
	IntegrationMethod _method;

	//Calculation state
	double _time;
	double _nextSnapshotTime;
	int _iter;
	std::vector<std::map<int, double> > _avgFitness;

public:	
	//Public accessors
	inline Vertex getVertex(int index) {
		return _vertices[index];
	};

	//Constructor
	ReplicatorModel(int nSp, std::vector<int> nSt, int nCells) {
		_nSpecies = nSp;
		_nStrategies.resize(_nSpecies);
		_shift.resize(_nSpecies);		
		if (nSt.size() != _nSpecies) throw 1;

		//Initilize internal data
		for (int i = 0; i<_nSpecies; i++) _nStrategies[i] = nSt[i];

		//Calculate variables shift in _values vector
		_nVariables = 0;
		for (int i = 0; i<_nSpecies; i++) {
			_shift[i] = _nVariables;
			_nVariables += _nStrategies[i];
		};

		//Initialize values vector
		_nCells = nCells;
		_values.resize(_nCells * _nVariables);

		//Create simple topology in R1
		double xMin = 0;
		double xMax = 1.0;
		double h = (xMax - xMin) / nCells;
		_vertices.resize(nCells);
		for (int i = 0; i<nCells; i++) {
			Vertex newVertex;
			newVertex.x = i * h + 0.5*h;
			newVertex.h = h;
			newVertex.neighbours.clear();
			if (i!=0) newVertex.neighbours.push_back(i-1);
			if (i!=nCells-1) newVertex.neighbours.push_back(i+1);
			_vertices[i] = newVertex;
		};
	};

	//Set payoff matices for each interaction
	//TO DO for now 2 species and assymetric game assumed
	void SetPayoffMatrix(int specieOne, int specieTwo, DenseMatrix payoffMatrix) {
		_games[specieOne][specieTwo] = payoffMatrix;
	};

	//Compute residual
	void ComputeResidual(std::vector<double>& R, const std::vector<double> U) {		
		//"conservative part"
		if (_method == IntegrationMethod::ExplicitEuler) {

			//We compute average fittnes for each specie and each game			
			_avgFitness.clear();			
			for (int i = 0; i<_nSpecies; i++) { //Species
				std::map<int, double> avgF;								
				avgF.clear();
				for(auto game : _games[i]) { //Games							
					double avgGame = 0;
					double sumS = 0;
					int sp2Ind = game.first;
					DenseMatrix& A = game.second;					
					for (int j = 0; j<_nCells; j++) { //Vertices
						double S = _vertices[j].h;						
						int player1 = _getIndex(j, i, 0);
						int player2 = _getIndex(j, sp2Ind, 0);	
						
						std::vector<double> u1(&U[player1], &U[player1] + _nStrategies[i]);
						std::vector<double> u2(&U[player2], &U[player2] + _nStrategies[sp2Ind]);
						double u1norm = 0;
						for (double& value : u1) u1norm += value;
						double u2norm = 0;
						for (double& value : u2) u2norm += value;
						//Assert
						assert(u1norm > 0.0);
						assert(u2norm > 0.0);
						double localFitness = mult(&U[player1], &(A * (&U[player2]))[0], _nStrategies[i]);
						//localFitness /= u1norm * u2norm;
						avgGame += S * localFitness;					
						sumS += S;
					};
					avgGame /= sumS;					
					avgF[sp2Ind] = avgGame;
				};
				_avgFitness.push_back(avgF);
			};

			//Compute residual for each cell
			double sumR = 0;
			R.resize(_nCells * _nVariables);
			//"convective" part of residual
			for (int cellInd = 0; cellInd < _nCells; cellInd++) {
				Vertex& cell = _vertices[cellInd];
				int startInd = _getIndex(cellInd, 0, 0);
				for (int i = 0; i<_nVariables; i++) R[startInd + i] = 0; //Initialize
				for (int iSp = 0; iSp < _nSpecies; iSp++) {
					for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {
						int index = _getIndex(cellInd, iSp, iSt);
						for (auto game : _games[iSp]) { //Games													
							int iSp2 = game.first;
							int player2 = _getIndex(cellInd, iSp2, 0);		
							DenseMatrix& A = game.second;	

							std::vector<double> u2(&U[player2], &U[player2] + _nStrategies[iSp2]);
							double u2norm = 0;
							for (double& value : u2) u2norm += value;
							//Assert
							assert(u2norm > 0.0);
							double localFittnes = (A * (&U[player2]))[iSt];
							//localFittnes /= u2norm;
							double u = U[index];						
							R[index] =  u * (localFittnes - _avgFitness[iSp][iSp2]);
							sumR += R[index];
						};
					};
				};
			};

			//"diffusion term"
			//between each pair of cells
			for (int cellInd = 0; cellInd < _nCells - 1; cellInd++) {				
				double dx = _vertices[cellInd].x - _vertices[cellInd+1].x;
				for (int iSp = 0; iSp < _nSpecies; iSp++) {
					for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {
						int leftIndex = _getIndex(cellInd, iSp, iSt);
						int rightIndex = _getIndex(cellInd+1, iSp, iSt);
						
						//Compute gradient
						double dU = U[leftIndex] - U[rightIndex];
						double dUdx = dU / dx;
						
						//Compute diffusion
						double d = _d[iSp][iSt];
						double S = 1.0;
						double dFlux = S * d * dUdx;

						//Distribute
						R[leftIndex] += dFlux / _vertices[cellInd].h ;
						R[rightIndex] -= dFlux / _vertices[cellInd+1].h;
					};
				};
			};

		}; //Explicit euler		

		return;
	};

	//Multiply implicit operator by vector containing cell values differences
	void MultiplyJacobianByVector(const double *x, double *r) {					
	};		

	//Explicit time step
	void ExplicitStep(double CFL, double maxTimeStep) {
		//Compute residuals
		ComputeResidual(_residual, _values);

		//Determine time step
		double dt = std::numeric_limits<double>::max();
		dt = maxTimeStep;

		//Maximum time step
		for (int i = 0; i<_nCells * _nVariables; i++) {
			double R = _residual[i];			
			double a = _values[i];			
			double dtMax = 0;
		/*	if (R > 0) {
				dtMax = (1.0 - a) / R;
				if (dtMax < dt) dt = dtMax;				
			};*/
			if (R < 0) {
				dtMax = -a / R;
				if (dtMax < dt) dt = dtMax;				
			};
		};

		//CFL coefficient
		dt *= CFL;

		//Consider next snapshot time
		dt = std::min(_nextSnapshotTime - _time, dt);

		//Update time		
		_time += dt; 

		//Distribute residual
		double norm = 0;
		for (int i = 0; i<_nCells; i++) {
			double h = _vertices[i].h;
			for (int j = 0; j<_nVariables; j++) {			
				double R = _residual[i * _nVariables + j];
				double newValue = _values[i * _nVariables + j] + R * dt;
				_values[i * _nVariables + j] = newValue;
				norm += newValue * h;
			};
		};	

		//Normilize distribution
		for (int iSp = 0; iSp < _nSpecies; iSp++) {
			double sum = 0;
			for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {
				for (int i = 0; i< _nCells; i++) {
					int startIndex = _getIndex(i, iSp, iSt);
					sum += _values[startIndex];
				};		
			};

			double A = 1.0 / sum;
			for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {
				for (int i = 0; i< _nCells; i++) {
					int startIndex = _getIndex(i, iSp, iSt);
					_values[startIndex] *= A;
				};	
			};
		};
	};

	//Initial conditions
	void SetInitialConditions(std::vector<std::vector<double> > initValues) {
		//Prepare initial values
		/*std::vector<double> inV;
		for (int i = 0; i<initValues.size(); i++) {
			for (int j = 0; j<initValues[i].size(); j++) {
				inV.push_back(initValues[i][j]);
			};
		};	*/	

		//Write to each vertex
 		for (int i = 0; i< _nCells; i++) {
			double x = _vertices[i].x;
			int startIndex = _getIndex(i, 0, 0);
			for (int j = 0; j<_nVariables; j++) {
				_values[startIndex + j] = initValues[i][j];//inV[j];
			};
		};
	};

	void SetInitialConditions(std::vector<double> avg, std::function<double(double)> dist) {
		//Write to each vertex
		double sum = 0;
 		for (int i = 0; i< _nCells; i++) {
			double x = _vertices[i].x;
			double cDistribution = dist(_vertices[i].x);
			int startIndex = _getIndex(i, 0, 0);
			for (int j = 0; j<_nVariables; j++) {
				_values[startIndex + j] = cDistribution;
				sum += cDistribution;
			};
		};

		for (int iSp = 0; iSp < _nSpecies; iSp++) {
			double sum = 0;
			double avgSum = 0;
			for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {
				int index = _getIndex(0, iSp, iSt);
				avgSum += avg[index];
				for (int i = 0; i< _nCells; i++) {
					int startIndex = _getIndex(i, iSp, iSt);
					sum += _values[startIndex];
				};		
			};

			double A = avgSum / sum;
			for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {
				for (int i = 0; i< _nCells; i++) {
					int index = _getIndex(0, iSp, iSt);
					int startIndex = _getIndex(i, iSp, iSt);
					_values[startIndex] *= avg[index] * A;
				};	
			};
		};		

		//Normilize
		/*double norm = 1.0 * _nSpecies;
		double A = norm / sum;
		for (int i = 0; i< _nCells; i++) {
			int startIndex = _getIndex(i, 0, 0);
			for (int j = 0; j<_nVariables; j++) {
				_values[startIndex + j] *= avg[j] * A;
			};
		};*/
	};

	//Set diffusion coefficients
	//iSp - Index of specie
	//iSt - Index of strategy
	//value - Diffusion coeff value
	void SetDiffusionCoefficient(int iSp, int iSt, double value) {
		_d[iSp][iSt] = value;
	};

	//Set method parameters
	void SetMethod(IntegrationMethod method) {
		_method = method;
	};
	
	//Set 
	void SetParameter(std::string parname, std::string value) {		
	};

	//Run calculation
	void RunCalculation(int maxIter, double startTime, double snapshotTime, double maxTime, double CFL, double maxTimeStep, bool verbose = true, bool output = true, std::string fname = "out.dat") {
		_maxIter = maxIter;
		_maxTime = maxTime;
		_time = startTime;	
		_nextSnapshotTime = startTime;
		std::ofstream ofs;
		std::ofstream ofsTime("history.dat");
		std::ofstream ofsConvergence("historyConvergence.dat");
		std::ofstream ofsX("x.txt");
		std::ofstream ofsValues("u.txt");
		std::stringstream uname;
		if (output) {			
			ofs.open(fname, std::ios_base::out);
			//Header
			ofs<<"VARIABLES= \"t\" ";
			ofs<<"\"x\" ";			
			for (int iSp = 0; iSp < _nSpecies; iSp++) {
				for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {					
					uname.str(std::string());
					uname<<"u_"<<iSp<<"_"<<iSt;
					ofs<<"\""<<uname.str()<<"\""<<" ";			
				};
			};						
			ofs<<"\n";

			//Convergence history file header
			ofsConvergence<<"VARIABLES= \"t\" ";
			ofsConvergence<<"\"total_residual\" ";
			ofsConvergence<<"\"total_sd\" ";			
			for (int iSp = 0; iSp < _nSpecies; iSp++) {
				for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {					
					uname.str(std::string());
					uname<<"avg_u_"<<iSp<<"_"<<iSt;
					ofsConvergence<<"\""<<uname.str()<<"\""<<" ";			
				};
			};				
			for (int iSp = 0; iSp < _nSpecies; iSp++) {
				for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {					
					uname.str(std::string());
					uname<<"sd_u_"<<iSp<<"_"<<iSt;
					ofsConvergence<<"\""<<uname.str()<<"\""<<" ";			
				};
			};			
			ofsConvergence<<std::endl;

			ofsTime<<"VARIABLES= \"t\" ";
			for (int iSp = 0; iSp < _nSpecies; iSp++) {
				for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {					
					uname.str(std::string());
					uname<<"u_"<<iSp<<"_"<<iSt;
					ofsTime<<"\""<<uname.str()<<"\""<<" ";			
				};
			};		

			//Average fitness for each game
			for (auto p1 : _games) {
				int iSp1 = p1.first;
				for (auto p2 : p1.second) {
					int iSp2 = p2.first;
					uname.str(std::string());
					uname<<"f_"<<iSp1<<"_"<<iSp2;
					ofsTime<<"\""<<uname.str()<<"\""<<" ";
				};
			};

			ofsTime<<"\n";

			//Coords to file
			for (int i = 0; i<_nCells; i++) {
				ofsX<<_vertices[i].x<<" ";
			};			
		};

		//Main computational cycle
		for (_iter = 0; _iter<_maxIter; _iter++) {
			double prevTime = _time;			
			if (_method == ExplicitEuler) {
				ExplicitStep(CFL, maxTimeStep);
			};

			//Compute L2 norm of residual
			double L2norm = 0;
			for (int i = 0; i< _nCells; i++) {
				L2norm += pow(_residual[i],2);
			};
			L2norm = sqrt(L2norm);

			//Compute average values and standart deviations
			std::vector<double> sd(_nVariables);
			std::vector<double> avg(_nVariables);
			for (int iSp = 0; iSp < _nSpecies; iSp++) {
				for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {
					int index = _getIndex(0, iSp, iSt);
					
					avg[index] = 0;
					double totalS = 0;
					for (int cellInd = 0; cellInd < _nCells; cellInd++) {		
						int startIndex = _getIndex(cellInd, iSp, iSt);
						Vertex& v = getVertex(cellInd);
						double S = v.h;
						totalS += S;
						avg[index] += S * _values[startIndex];
					};
					avg[index] /= totalS;
					
					sd[index] = 0;
					totalS = 0;
					for (int cellInd = 0; cellInd < _nCells; cellInd++) {		
						int startIndex = _getIndex(cellInd, iSp, iSt);
						Vertex& v = getVertex(cellInd);
						double S = v.h;
						totalS += S;
						sd[index] += std::pow(S * (_values[startIndex] - avg[index]), 2); 	
					};
					sd[index] = std::sqrt(sd[index]) / totalS;
				};
			};


			//Output convergence history
			double total_sd = 0;
			for (double& sdval : sd) total_sd += sdval;
			ofsConvergence<<_time<<" ";
			ofsConvergence<<L2norm<<" ";
			ofsConvergence<<total_sd<<" ";	
			
			for (int iSp = 0; iSp < _nSpecies; iSp++) {
				for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {
					int index = _getIndex(0, iSp, iSt);
					ofsConvergence<<avg[index]<<" ";			
				};
			};				
			for (int iSp = 0; iSp < _nSpecies; iSp++) {
				for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {
					int index = _getIndex(0, iSp, iSt);
					ofsConvergence<<sd[index]<<" ";			
				};
			};		
			ofsConvergence<<std::endl;

			//Output data snapshot
			if ((output) && (_time == _nextSnapshotTime)) {				
				//Save values to file
				//ofs<<"ZONE T=\"Time = " << _time << " s\", I=" << _nCells << ", J=1, DATAPACKING=POINT\n";
				for (int i = 0; i<_nCells; i++) {				
					ofs<<_time<<" ";
					ofs<<_vertices[i].x<<" ";
					for (int j = 0; j<_nSpecies; j++) {
						for (int k = 0; k<_nStrategies[j]; k++) {
							int index = _getIndex(i, j, k);
							ofs<<_values[index]<<" ";
						};
					};
					ofs<<"\n";
				};

				//Save time to file
				ofsTime<<_time<<" ";

				//Save average values to file				
				for (int j = 0; j<_nSpecies; j++) {
					for (int k = 0; k<_nStrategies[j]; k++) {
						double sum = 0;
						for (int i = 0; i<_nCells; i++) {														
							int index = _getIndex(i, j, k);
							sum+=_values[index];
						};
						sum /= _nCells;
						ofsTime<<sum<<" ";
					};					
				};

				//Save average fitness for each game
				for (auto p1 : _games) {
					int iSp1 = p1.first;
					for (auto p2 : p1.second) {
						int iSp2 = p2.first;
						ofsTime<<_avgFitness[iSp1][iSp2]<<" ";
					};
				};
				
				ofsTime<<"\n";
				
				//Save values to file				
				for (int j = 0; j<_nSpecies; j++) {
					for (int k = 0; k<_nStrategies[j]; k++) {
						for (int i = 0; i<_nCells; i++) {														
							int index = _getIndex(i, j, k);
							ofsValues<<_values[index]<<" ";
						};
						ofsValues<<"\n";
					};					
				};							
			};

			//Update next snapshot time
			if (_time == _nextSnapshotTime) _nextSnapshotTime += snapshotTime;

			double avgFitness = _avgFitness[0][1];
			if (verbose) std::cout<<"Iteration "<<_iter<<"; ";
			if (verbose) std::cout<<"Time : " <<_time << "; "<<" dt = "<< _time - prevTime <<" L2-residual = "<<L2norm<<" Average fitness = "<<avgFitness<<"\n";
			if (L2norm < 1e-20) {
				if (verbose) std::cout<<"Iterations stoped : Steady state reached, judging by L2-residual norm.\n";
				break;
			};
			if (_time > _maxTime) {
				if (verbose) std::cout<<"Iterations stoped : Maximum time reached\n";
				break;
			};
		};
		if ((verbose) && (_iter == _maxIter)) std::cout<<"Iterations stoped : Maximum number of iterations reached\n";

		if (output) {
			ofs.close();
			ofsX.close();
			ofsTime.close();
			ofsValues.close();
			ofsConvergence.close();
		};
	};

	//Export data to tecplot
	void WriteToTecplot() {
	};

	void WritePhasePortrait(std::string fname, int nRes) {
		std::ofstream ofs(fname);

		//Output part header
		std::stringstream uname;
		ofs<<"VARIABLES= ";
		for (int iSp = 0; iSp < _nSpecies; iSp++) {
			for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {					
				uname.str(std::string());
				uname<<"u_"<<iSp<<"_"<<iSt;
				ofs<<"\""<<uname.str()<<"\""<<" ";			
			};
		};		
		for (int iSp = 0; iSp < _nSpecies; iSp++) {
			for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {					
				uname.str(std::string());
				uname<<"du_"<<iSp<<"_"<<iSt;
				ofs<<"\""<<uname.str()<<"_dt\""<<" ";			
			};
		};		
		ofs<<std::endl;

		//Compute residual in every point averaged over 
		std::vector<double> u(0);

		std::function<void(int, int)> outputPoint = [&](int iSp, int iSt) {
			//Iterative part
			if (iSp < this->_nSpecies) {
				double uvalue = 0;
				//If it's the last strategy use normalization criterium to compute value
				if (iSt == _nStrategies[iSp] - 1) {
					uvalue = 1.0;
					for (int i = 1; i < _nStrategies[iSp]; i++) uvalue -= *(std::end(u) - i);

					//If normilized value is negative cut that execution branch
					if (uvalue < 0) return;
					
					//If it's the last strategy move to the next specie
					u.push_back(uvalue);
					outputPoint(iSp+1, 0);
					u.pop_back();
					return;
				};

				for (int i = 0; i < nRes; i++) {
					//Compute variable variants
					uvalue = 1.0 * i / nRes;

					//If possible move to the next strategy
					if (iSt < _nStrategies[iSp] - 1) {
						u.push_back(uvalue);
						outputPoint(iSp, iSt+1);
						u.pop_back();
					};
				};

				return;
			};

			//Calculation part
			std::vector<double> U(_nCells * _nVariables);
			std::vector<double> R(_nCells * _nVariables);
			//Write values
			for (int cellInd = 0; cellInd < _nCells; cellInd++) {
				//Get cell info
				Vertex& cell = _vertices[cellInd];
				int startInd = _getIndex(cellInd, 0, 0);
				for (int i = 0; i < _nVariables; i++) U[startInd + i] = u[i];
			};

			ComputeResidual(R, U); //Compute residual

			std::vector<double> Ravg(_nVariables);
			for (int i = 0; i<_nVariables; i++) Ravg[i] = 0; //Initialize
			for (int cellInd = 0; cellInd < _nCells; cellInd++) {
				//Get cell info
				Vertex& cell = _vertices[cellInd];
				int startInd = _getIndex(cellInd, 0, 0);
				for (int iSp = 0; iSp < _nSpecies; iSp++) {
					for (int iSt = 0; iSt < _nStrategies[iSp] - 1; iSt++) {
						int index = _getIndex(cellInd, iSp, iSt);
						//
						for (int i = 0; i<_nVariables; i++) Ravg[i] += R[startInd + i]; //Sum
					};
				};
			};
			for (int i = 0; i<_nVariables; i++) Ravg[i] /= 1.0 * _nCells; //Average

			//Output state variables
			for (int i = 0; i<_nVariables; i++) {
				ofs<<u[i]<<" ";
			};

			//Output averaged over cells residual
			for (int i = 0; i<_nVariables; i++) {
				ofs<<-Ravg[i]<<" ";
			};
			ofs<<std::endl;

			//Leave
			return;
		}; //outputPoint()
	
		//Make call
		outputPoint(0,0);
		
		//CLose output stream
		ofs.close();
	};
};

void mult( ReplicatorModel &A, const double *v, double *w ) {
	A.MultiplyJacobianByVector(v, w);
};

#endif