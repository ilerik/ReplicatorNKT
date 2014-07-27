#include <fstream>
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <math.h>

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
	int _iter;

public:	
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
			std::vector<std::map<int, double> > avgFittnes(_nSpecies); //f1 and f2
			avgFittnes.clear();
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
						//double *v = 
						avgGame += S * mult(&U[player1], &(A * (&U[player2]))[0], _nStrategies[i]);					
						sumS += S;
					};
					avgGame /= sumS;					
					avgF[sp2Ind] = avgGame;
				};
				avgFittnes.push_back(avgF);
			};

			//Compute residual for each cell
			R.resize(_nCells * _nVariables);
			//"convective" part of residual
			for (int cellInd = 0; cellInd < _nCells; cellInd++) {
				int startInd = _getIndex(cellInd, 0, 0);
				for (int i = 0; i<_nVariables; i++) R[startInd + i] = 0; //Initialize
				for (int iSp = 0; iSp < _nSpecies; iSp++) {
					for (int iSt = 0; iSt < _nStrategies[iSp]; iSt++) {
						int index = _getIndex(cellInd, iSp, iSt);
						for (auto game : _games[iSp]) { //Games													
							int iSp2 = game.first;
							int player2 = _getIndex(cellInd, iSp2, 0);		
							DenseMatrix& A = game.second;	
							double localFittnes = (A * (&U[player2]))[iSt];
							R[index] = U[index] * (localFittnes - avgFittnes[iSp][iSp2]);							
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
						const double d = 0.1;
						double dFlux = d * dUdx;

						//Distribute
						R[leftIndex] += dFlux;
						R[rightIndex] -= dFlux;
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
	void ExplicitStep() {
		double CFL = 0.0005;	

		//Compute residuals
		ComputeResidual(_residual, _values);

		//Determine time step
		double dt = std::numeric_limits<double>::max();

		//Maximum time step
		for (int i = 0; i<_nCells * _nVariables; i++) {
			double R = _residual[i];			
			double a = _values[i];			
			double dtMax = 0;
			if (R > 0) {
				dtMax = (1.0 - a) / R;
				if (dtMax < dt) dt = dtMax;				
			};
			if (R < 0) {
				dtMax = -a / R;
				if (dtMax < dt) dt = dtMax;				
			};
		};

		//Update time
		dt *= CFL;
		_time += dt; 

		//Distribute residual
		for (int i = 0; i<_nCells * _nVariables; i++) {
			double R = _residual[i];
			_values[i] += R * dt;
		};		
	};

	//Initial conditions
	void SetInitialConditions(std::vector<std::vector<double> > initValues) {
		//Prepare initial values
		std::vector<double> inV;
		for (int i = 0; i<initValues.size(); i++) {
			for (int j = 0; j<initValues[i].size(); j++) {
				inV.push_back(initValues[i][j]);
			};
		};		

		//Write to each vertex
		for (int i = 0; i< _nCells; i++) {
			double x = _vertices[i].x;
			int startIndex = _getIndex(i, 0, 0);
			for (int j = 0; j<_nVariables; j++) {
				_values[startIndex + j] = inV[j];
			};
		};
	};

	//Set method parameters
	void SetMethod(IntegrationMethod method) {
		_method = method;
	};
	
	//Set 
	void SetParameter(std::string parname, std::string value) {		
	};

	//Run calculation
	void RunCalculation(int maxIter, double startTime, double maxTime, bool verbose = true, bool output = true, std::string fname = "out.dat") {
		_maxIter = maxIter;
		_maxTime = maxTime;
		_time = startTime;
		std::ofstream ofs;
		if (output) {
			ofs.open(fname, std::ios_base::out);
			//Header
			ofs<<"VARIABLES= \"t\"\n";
			ofs<<"\"x\"\n";
			ofs<<"\"u0\"\n";
			ofs<<"\"u1\"\n";
			ofs<<"\"v0\"\n";
			ofs<<"\"v1\"\n";			
		};
		//Main computational cycle
		for (_iter = 0; _iter<_maxIter; _iter++) {
			double prevTime = _time;			
			if (_method == ExplicitEuler) {
				ExplicitStep();
			};

			if (output) {
				ofs<<"ZONE T=\"Time = " << _time << " s\", I=" << _nCells << ", J=1, DATAPACKING=POINT\n";
				for (int i = 0; i<_nCells; i++) {
					ofs<<_time<<"\n";
					ofs<<_vertices[i].x<<"\n";
					for (int j = 0; j<_nSpecies; j++) {
						for (int k = 0; k<_nStrategies[j]; k++) {
							int index = _getIndex(i, j, k);
							ofs<<_values[index]<<"\n";
						};
					};
				};
			};

			if (verbose) std::cout<<"Iteration "<<_iter<<"; ";
			if (verbose) std::cout<<"Time : " <<_time << "; "<<" dt = "<< _time - prevTime <<"\n";
			if (_time > _maxTime) {
				if (verbose) std::cout<<"Iterations stoped : Maximum time reached\n";
				break;
			};
		};
		if ((verbose) && (_iter == _maxIter)) std::cout<<"Iterations stoped : Maximum number of iterations reached\n";

		if (output) ofs.close();
	};

	//Export data to tecplot
	void WriteToTecplot() {
	};
};

void mult( ReplicatorModel &A, const double *v, double *w ) {
	A.MultiplyJacobianByVector(v, w);
};

int main(int argc, char *argv[]) {
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

	kernel.RunCalculation(100000, 0.0, maxTime);
	return 0;
};
