#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <String>
#include <ctime>
#include <Math.h>
#include <fstream>
#include <list>
#include <vector>
#include <iomanip>

using namespace std;

class Matrix
{
public:
	double** elements;
	int col;
	int row;

	Matrix()
	{
		elements = nullptr;
		col = 0;
		row = 0;
	}

	Matrix(int col, int row)
	{
		elements = new double*[row];
		for (int i = 0; i < row; i++)
		{
			elements[i] = new double[col];
		}

		SetToNull();
		this->col = col;
		this->row = row;
	}

	Matrix(const Matrix& _matrix)
	{
		col = _matrix.col;
		row = _matrix.row;

		elements = new double*[row];
		for (int j = 0; j < row; j++)
		{
			elements[j] = new double[col];
		}

		for (int j = 0; j < row; j++)
		{
			for (int i = 0; i < col; i++)
			{
				elements[j][i] = _matrix.elements[j][i];
			}
		}
	}

	~Matrix()
	{
		for (int i = 0; i < row; i++)
		{
			delete[] elements[i];
		}

		delete[] elements;
	}
	
	void SetToNull()
	{
		for (int j = 0; j < row; j++)
		{
			for (int i = 0; i < col; i++)
			{
				elements[j][i] = 0;
			}
		}
	}

	//??????? ???????????????? ???????
	void Transpose()
	{
		double** tElements;
		tElements = new double*[col];
		for (int i = 0; i < col; i++)
		{
			tElements[i] = new double[row];
			for (int j = 0; j < row; j++)
			{
				tElements[i][j] = elements[j][i];
			}
		}

		for (int i = 0; i < row; i++)
		{
			delete[] elements[i];
		}

		delete[] elements;

		int c = 0;
		c = col;
		col = row;
		row = c;

		elements = tElements;
	}
};

void PrintMatrix(Matrix matrix)
{
	for (int j = 0; j < matrix.row; j++)
	{
		for (int i = 0; i < matrix.col; i++)
		{
			cout << setw(5) << matrix.elements[j][i] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

double AvgTrustedInterval(double& avg, double*& times, int& cnt)
{
	double sd = 0, newAVg = 0;
	int newCnt = 0;
	for (int i = 0; i < cnt; i++)
	{
		sd += (times[i] - avg) * (times[i] - avg);
	}
	sd /= (cnt - 1.0);
	sd = sqrt(sd);
	for (int i = 0; i < cnt; i++)
	{
		if (avg - sd <= times[i] && times[i] <= avg + sd)
		{
			newAVg += times[i];
			newCnt++;
		}
	}
	if (newCnt == 0) newCnt = 1;
	return newAVg / newCnt;
}

typedef double(*PiFunctTempl)(long, double&);

double CalculatePiConsistenty(long num_steps, double& time)
{
	double step, pi, x, sum = 0.0;
	int i;
	step = 1.0 / (double)num_steps;
	double startTime, endTime;
	startTime = omp_get_wtime();
	for (i = 0; i < num_steps; i++)
	{
		x = (i + 0.5) * step;
		sum = sum + 4.0 / (1.0 + x * x);
	}
	pi = step * sum;
	endTime = omp_get_wtime();

	time = endTime - startTime;
	return pi;
}

double CalculatePiParallelForStatic(long num_steps, double& time)
{
	double step, pi, x, sum = 0.0;
	step = 1.0 / (double)num_steps;
	double startTime, endTime;
	startTime = omp_get_wtime();
//#pragma omp parallel
#pragma omp parallel for schedule(static, num_steps / 12) reduction(+:sum) private(x)
	for (int i = 0; i < num_steps; i++)
	{
		x = (i + 0.5) * step;
		sum = sum + 4.0 / (1.0 + x * x);
	}
		
	pi = step * sum;
	endTime = omp_get_wtime();

	time = endTime - startTime;
	return pi;
}

double CalculatePiParallelForDynamic(long num_steps, double& time)
{
	double step, pi, x, sum = 0.0;
	step = 1.0 / (double)num_steps;
	double startTime, endTime;
	startTime = omp_get_wtime();

#pragma omp parallel for schedule(dynamic, num_steps / 6) reduction(+:sum) private(x)
	for (int i = 0; i < num_steps; i++)
	{
		x = (i + 0.5) * step;
		sum = sum + 4.0 / (1.0 + x * x);
	}

	pi = step * sum;
	endTime = omp_get_wtime();

	time = endTime - startTime;
	return pi;
}

double CalculatePiParallelForGuided(long num_steps, double& time)
{
	double step, pi, x, sum = 0.0;
	step = 1.0 / (double)num_steps;
	double startTime, endTime;
	startTime = omp_get_wtime();
#pragma omp parallel for schedule(guided, 10000) reduction(+:sum) private(x)
	for (int i = 0; i < num_steps; i++)
	{
		x = (i + 0.5) * step;
		sum = sum + 4.0 / (1.0 + x * x);
	}
	pi = step * sum;
	endTime = omp_get_wtime();

	time = endTime - startTime;
	return pi;
}

double CalculatePiParallelSections(long num_steps, double& time)
{
	double step, pi, sum = 0.0;
	step = 1.0 / (double)num_steps;
	int threadsNum;
		
#pragma omp parallel
	{
		threadsNum = omp_get_num_threads();
	}
	double startTime, endTime;
	int s1 = num_steps / threadsNum, s2 = num_steps * 2 / threadsNum, s3 = num_steps * 3 / threadsNum;
	
	if (threadsNum == 3)
		s3 = num_steps;
	
	double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
	startTime = omp_get_wtime();
#pragma omp sections
	{
#pragma omp section
		{
			double x;
			for (int i = 0; i < s1; i++)
			{
				x = (i + 0.5) * step;
				sum1 = sum1 + 4.0 / (1.0 + x * x);
			}
		}
#pragma omp section
		{
			double x;
			for (int i = s1; i < s2; i++)
			{
				x = (i + 0.5) * step;
				sum2 = sum2 + 4.0 / (1.0 + x * x);
			}
		}
#pragma omp section
		{
			double x;
			if (threadsNum > 2)
				for (int i = s2; i < s3; i++)
				{
					x = (i + 0.5) * step;
					sum3 = sum3 + 4.0 / (1.0 + x * x);
				}
		}
#pragma omp section
		{
			double x;
			if (threadsNum > 3)
				for (int i = s3; i < num_steps; i++)
				{
					x = (i + 0.5) * step;
					sum4 = sum4 + 4.0 / (1.0 + x * x);
				}
		}
	}
	sum = sum1 + sum2 + sum3 + sum4;
	pi = step * sum;

	endTime = omp_get_wtime();
	time = endTime - startTime;

	return pi;
}

double CalculatePiMyAlgorithm(long num_steps, double& time)
{
	double step, pi, x, sum = 0.0;
	step = 1.0 / (double)num_steps;
	int part1 = num_steps / 2;
	double startTime, endTime;
	startTime = omp_get_wtime();
	//#pragma omp parallel


	int threadsNum;

#pragma omp parallel
	{
		threadsNum = omp_get_num_threads();
	}
	int s1 = part1 / threadsNum, s2 = part1 * 2 / threadsNum, s3 = part1 * 3 / threadsNum;

	if (threadsNum == 3)
		s3 = part1;

	double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;

#pragma omp sections
	{
#pragma omp section
		{
			double x;
			for (int i = 0; i < s1; i++)
			{
				x = (i + 0.5) * step;
				sum1 = sum1 + 4.0 / (1.0 + x * x);
			}
		}
#pragma omp section
		{
			double x;
			for (int i = s1; i < s2; i++)
			{
				x = (i + 0.5) * step;
				sum2 = sum2 + 4.0 / (1.0 + x * x);
			}
		}
#pragma omp section
		{
			double x;
			if (threadsNum > 2)
				for (int i = s2; i < s3; i++)
				{
					x = (i + 0.5) * step;
					sum3 = sum3 + 4.0 / (1.0 + x * x);
				}
		}
#pragma omp section
		{
			double x;
			if (threadsNum > 3)
				for (int i = s3; i < part1; i++)
				{
					x = (i + 0.5) * step;
					sum4 = sum4 + 4.0 / (1.0 + x * x);
				}
		}
	}

#pragma omp parallel for schedule(static, num_steps / 12) reduction(+:sum) private(x)
	for (int i = part1; i < num_steps; i++)
	{
		x = (i + 0.5) * step;
		sum = sum + 4.0 / (1.0 + x * x);
	}

	sum += sum1 + sum2 + sum3 + sum4;
	pi = step * sum;
	endTime = omp_get_wtime();

	time = endTime - startTime;
	return pi;
}

double CalculatePiFuncTrustedTime(void* func, long num_steps, double& time, int iterations)
{
	double curtime = 0, avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	double* Times = new double[iterations];
	double pi = 0;
	for (int i = 0; i < iterations; i++)
	{
		pi = (*(PiFunctTempl)func)(num_steps, curtime);
		Times[i] = curtime;
		avgTime += curtime;
	}
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	time = avgTimeT;
	return pi;
}

void TaskCalculatePi()
{
	string* funcsNames = new string[6]{ "???????????????? ??????????", "???????????? ?????????? FOR (static)", "???????????? ?????????? FOR (dynamic)", "???????????? ?????????? FOR (guided)", "???????????? ?????????? Sections", "???? ????????" };
	void** funcs = new void*[6]{ CalculatePiConsistenty, CalculatePiParallelForStatic, CalculatePiParallelForDynamic, CalculatePiParallelForGuided, CalculatePiParallelSections, CalculatePiMyAlgorithm };

	static long* nums_steps = new long[4]{ 1000000, 5000000, 10000000, 15000000 };

	std::ofstream resultsFile;
	resultsFile.open("Task1Results.csv", std::ios_base::app);
	resultsFile << "???????;??????;??1 - " << nums_steps[0] << ";??2 - " << nums_steps[1] << ";??3 - " << nums_steps[2] << ";??4 - " << nums_steps[3] << ";\n";

	double time, pi;
	for (int i = 0; i < 6; i++)
	{
		resultsFile << funcsNames[i] << ";";
		for (int t = 2; t < 5; t++)
		{
			if (t != 2)
			{
				resultsFile << ";";
			}

			if (i != 0)
			{
				resultsFile << t << ";";
			}
			else
			{
				resultsFile << ";";
			}

			omp_set_num_threads(t);
			cout << "????? ???????:" << t << endl;
			for (int d = 0; d < 4; d++)
			{
				cout << "??:" << nums_steps[d] << endl;
				time = 0;
				pi = CalculatePiFuncTrustedTime(funcs[i], nums_steps[d], time, 50);
				cout << funcsNames[i] << ". ????? ?? = " << pi << ". ??????? ?????????: " << time * 1000 << endl;
				resultsFile << time * 1000 << "^;";
			}

			resultsFile << "\n";
			if (i == 0)
				break;

			cout << endl;
		}
	}

	resultsFile.close();
}

typedef double(*FillFunctTempl)(Matrix, double&);

void FillLineConsistently(Matrix matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	for (int j = 0; j < matrix.row; j++)
	{
		for (int i = 0; i < matrix.col; i++)
		{
			matrix.elements[j][i] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
		}
	}
	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void FillLineParallelForStatic(Matrix matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
#pragma omp parallel for schedule(static, matrix.row / 12)
	for (int j = 0; j < matrix.row; j++)
	{
		for (int i = 0; i < matrix.col; i++)
		{
			matrix.elements[j][i] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
		}
	}
	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void FillLineParallelForDynamic(Matrix matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, matrix.row / 6)
	for (int j = 0; j < matrix.row; j++)
	{
		for (int i = 0; i < matrix.col; i++)
		{
			matrix.elements[j][i] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
		}
	}
	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void FillLineParallelSections(Matrix matrix, double& time)
{
	int threadsNum;
#pragma omp parallel
	{
		threadsNum = omp_get_num_threads();
	}
	double startTime, endTime;
	int s1 = matrix.row / threadsNum, s2 = matrix.row * 2 / threadsNum, s3 = matrix.row * 3 / threadsNum;

	if (threadsNum == 3)
		s3 = matrix.row;

	startTime = omp_get_wtime();
#pragma omp sections
	{
#pragma omp section
		{
			for (int j = 0; j < s1; j++)
			{
				for (int i = 0; i < matrix.col; i++)
				{
					matrix.elements[j][i] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
				}
			}
		}
#pragma omp section
		{
			for (int j = s1; j < s2; j++)
			{
				for (int i = 0; i < matrix.col; i++)
				{
					matrix.elements[j][i] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
				}
			}
		}
#pragma omp section
		{
			if (threadsNum > 2)
			{
				for (int j = s2; j < s3; j++)
				{
					for (int i = 0; i < matrix.col; i++)
					{
						matrix.elements[j][i] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
					}
				}
			}
		}
#pragma omp section
		{
			if (threadsNum > 3)
			{
				for (int j = s3; j < matrix.row; j++)
				{
					for (int i = 0; i < matrix.col; i++)
					{
						matrix.elements[j][i] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
					}
				}
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void GetBlockStartPosition(int i, int m, int sizeX, int& shiftX, int& shiftY)
{
	shiftX = (i * sizeX) % m;
	shiftY = trunc((i * sizeX) / m) * shiftY;
}

void FillBlockConsistently(Matrix matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	int blockSizeX = 5;
	int blockSizeY = 5;
	int blocksAmount = matrix.row * matrix.col / (blockSizeX * blockSizeY);
	int shiftX = 0, shiftY = 0;
	int m = matrix.col;

	for (int i = 0; i < blocksAmount; i++)
	{
		GetBlockStartPosition(i, m, blockSizeX, shiftX, shiftY);
		for (int t = 0; t < blockSizeY; t++)
		{
			int tIndex = t + shiftY;
			for (int d = 0; d < blockSizeX; d++)
			{
				int dIndex = d + shiftX;
				matrix.elements[tIndex][dIndex] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1);
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void FillBlockParallelForStatic(Matrix matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	int blockSizeX = 5;
	int blockSizeY = 5;
	int blocksAmount = matrix.row * matrix.col / (blockSizeX * blockSizeY);
	int shiftX = 0, shiftY = 0;
	int m = matrix.col;

#pragma omp parallel for schedule(static, blocksAmount / 12)
	for (int i = 0; i < blocksAmount; i++)
	{
		GetBlockStartPosition(i, m, blockSizeX, shiftX, shiftY);
		for (int t = 0; t < blockSizeY; t++)
		{
			int tIndex = t + shiftY;
			for (int d = 0; d < blockSizeX; d++)
			{
				int dIndex = d + shiftX;
				matrix.elements[tIndex][dIndex] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1);
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void FillBlockParallelForDynamic(Matrix matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	int blockSizeX = 5;
	int blockSizeY = 5;
	int blocksAmount = matrix.row * matrix.col / (blockSizeX * blockSizeY);
	int shiftX = 0, shiftY = 0;
	int m = matrix.col;

#pragma omp parallel for schedule(dynamic, blocksAmount / 6)
	for (int i = 0; i < blocksAmount; i++)
	{
		GetBlockStartPosition(i, m, blockSizeX, shiftX, shiftY);
		for (int t = 0; t < blockSizeY; t++)
		{
			int tIndex = t + shiftY;
			for (int d = 0; d < blockSizeX; d++)
			{
				int dIndex = d + shiftX;
				matrix.elements[tIndex][dIndex] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1);
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void CalculateFillFuncTrustedTime(void* func, int m, int n, double& time, int iterations)
{
	Matrix matrixA(m, n);
	double curtime = 0, avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	double* Times = new double[iterations];

	cout << "[";
	for (int i = 0; i < iterations; i++)
	{
		(*(FillFunctTempl)func)(matrixA, curtime);
		Times[i] = curtime;
		avgTime += curtime;
		cout << "*";
	}
	cout << "]\n";
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	time = avgTimeT;
}

typedef double(*MultFunctTempl)(Matrix, Matrix, Matrix, double&);

void MultiplicationConsistently(Matrix matrixA, Matrix matrixB, Matrix matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	matrixB.Transpose();

	//? ????? ?????????? ??????? ?? ??????? ??????? C
	for (int j = 0; j < matrixB.row; j++)
	{
		//? ????? ?????????? ??????? ?? ???????? ??????? C
		for (int i = 0; i < matrixA.row; i++)
		{
			matrixC.elements[j][i] = 0;
			//? ????? ?????????? ???????????? ????????? ?????? A ? B ?? ??????????????? ??????
			for (int t = 0; t < matrixA.col; t++)
			{
				matrixC.elements[j][i] += matrixA.elements[i][t] * matrixB.elements[j][t];
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void MultiplicationParallelForStatic(Matrix matrixA, Matrix matrixB, Matrix matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	matrixB.Transpose();

#pragma omp parallel for schedule(static, matrixA.row / 12)
	//? ????? ?????????? ??????? ?? ??????? ??????? C
	for (int j = 0; j < matrixB.row; j++)
	{
		//? ????? ?????????? ??????? ?? ???????? ??????? C
		for (int i = 0; i < matrixA.row; i++)
		{
			matrixC.elements[j][i] = 0;
			//? ????? ?????????? ???????????? ????????? ?????? A ? B ?? ??????????????? ??????
			for (int t = 0; t < matrixA.col; t++)
			{
				matrixC.elements[j][i] += matrixA.elements[i][t] * matrixB.elements[j][t];
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void MultiplicationParallelForDynamic(Matrix matrixA, Matrix matrixB, Matrix matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	matrixB.Transpose();

#pragma omp parallel for schedule(dynamic, matrixA.row / 6)
	//? ????? ?????????? ??????? ?? ??????? ??????? C
	for (int j = 0; j < matrixB.row; j++)
	{
		//? ????? ?????????? ??????? ?? ???????? ??????? C
		for (int i = 0; i < matrixA.row; i++)
		{
			matrixC.elements[j][i] = 0;
			//? ????? ?????????? ???????????? ????????? ?????? A ? B ?? ??????????????? ??????
			for (int t = 0; t < matrixA.col; t++)
			{
				matrixC.elements[j][i] += matrixA.elements[i][t] * matrixB.elements[j][t];
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void MultiplicationParallelSections(Matrix matrixA, Matrix matrixB, Matrix matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	matrixB.Transpose();

	int numThreads = 0;
	int s1, s2, s3;

#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
	}
	
	s1 = matrixB.row / numThreads;
	s2 = matrixB.row * 2 / numThreads;
	s3 = matrixB.row * 3 / numThreads;

#pragma omp sections
	{
#pragma omp section
		{
			//? ????? ?????????? ??????? ?? ??????? ??????? C
			for (int j = 0; j < s1; j++)
			{
				//? ????? ?????????? ??????? ?? ???????? ??????? C
				for (int i = 0; i < matrixA.row; i++)
				{
					matrixC.elements[j][i] = 0;
					//? ????? ?????????? ???????????? ????????? ?????? A ? B ?? ??????????????? ??????
					for (int t = 0; t < matrixA.col; t++)
					{
						matrixC.elements[j][i] += matrixA.elements[i][t] * matrixB.elements[j][t];
					}
				}
			}
		}
#pragma omp section
		{
			//? ????? ?????????? ??????? ?? ??????? ??????? C
			for (int j = s1; j < s2; j++)
			{
				//? ????? ?????????? ??????? ?? ???????? ??????? C
				for (int i = 0; i < matrixA.row; i++)
				{
					matrixC.elements[j][i] = 0;
					//? ????? ?????????? ???????????? ????????? ?????? A ? B ?? ??????????????? ??????
					for (int t = 0; t < matrixA.col; t++)
					{
						matrixC.elements[j][i] += matrixA.elements[i][t] * matrixB.elements[j][t];
					}
				}
			}
		}
#pragma omp section
		{
			if (numThreads > 2)
				//? ????? ?????????? ??????? ?? ??????? ??????? C
				for (int j = s2; j < s3; j++)
				{
					//? ????? ?????????? ??????? ?? ???????? ??????? C
					for (int i = 0; i < matrixA.row; i++)
					{
						matrixC.elements[j][i] = 0;
						//? ????? ?????????? ???????????? ????????? ?????? A ? B ?? ??????????????? ??????
						for (int t = 0; t < matrixA.col; t++)
						{
							matrixC.elements[j][i] += matrixA.elements[i][t] * matrixB.elements[j][t];
						}
					}
				}
		}
#pragma omp section
		{
			if (numThreads > 3)
				//? ????? ?????????? ??????? ?? ??????? ??????? C
				for (int j = s3; j < matrixB.row; j++)
				{
					//? ????? ?????????? ??????? ?? ???????? ??????? C
					for (int i = 0; i < matrixA.row; i++)
					{
						matrixC.elements[j][i] = 0;
						//? ????? ?????????? ???????????? ????????? ?????? A ? B ?? ??????????????? ??????
						for (int t = 0; t < matrixA.col; t++)
						{
							matrixC.elements[j][i] += matrixA.elements[i][t] * matrixB.elements[j][t];
						}
					}
				}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

int GetDimentionDegree2(int m, int n, int k)
{
	int maxNum = std::fmax(m, std::fmax(n, k));
	int i = 1;
	int log2Elements = 0;
	while (maxNum > log2Elements)
	{
		log2Elements = pow(2, i);
		i++;
	}

	return log2Elements;
}

Matrix*& MatrixExtension(Matrix matrix, int size) 
{
	Matrix* result = new Matrix(size, size);

	for (int j = 0; j < size; j++) {
		for (int i = 0; i < size; i++) {
			if (i < matrix.col && j < matrix.row)
				//?????? ??? ?????????? ??????? MatrixC ?? 4 ??????? 1 ????????
				result->elements[j][i] = matrix.elements[j][i];
			else
				result->elements[j][i] = 0;
		}
	}

	return result;
}

struct Position
{
	int x;
	int y;

	Position()
	{
		x = 0;
		y = 0;
	}

	Position(int _x, int _y)
	{
		x = _x;
		y = _y;
	}
};

Matrix SumMatrix(Matrix matrixA, Position posA, Matrix matrixB, Position posB, int size)
{
	Matrix matrixC(size, size);
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrixC.elements[j][i] = matrixA.elements[j + posA.y][i + posA.x] + matrixB.elements[j + posB.y][i + posB.x];
		}
	}

	return matrixC;
}

Matrix SubMatrix(Matrix matrixA, Position posA, Matrix matrixB, Position posB, int size)
{
	Matrix matrixC(size, size);
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrixC.elements[j][i] = matrixA.elements[j + posA.y][i + posA.x] - matrixB.elements[j + posB.y][i + posB.x];
		}
	}

	return matrixC;
}

Matrix SumTMatrix(Matrix matrixA, Matrix matrixB)
{
	int size = matrixA.col;
	Matrix matrixC(size, size);
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrixC.elements[j][i] = matrixA.elements[j][i] + matrixB.elements[j][i];
		}
	}

	return matrixC;
}

Matrix SubTMatrix(Matrix matrixA, Matrix matrixB)
{
	int size = matrixA.col;
	Matrix matrixC(size, size);
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrixC.elements[j][i] = matrixA.elements[j][i] - matrixB.elements[j][i];
		}
	}

	return matrixC;
}

Matrix GetPart(Matrix matrixA, Position posA, int size)
{
	Matrix matrixC(size, size);
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrixC.elements[j][i] = matrixA.elements[j + posA.y][i + posA.x];
		}
	}

	return matrixC;
}

Matrix UniteMatrix(Matrix matrix11, Matrix matrix12, Matrix matrix21, Matrix matrix22)
{
	int size = matrix11.col;
	Matrix matrixC(size * 2, size * 2);
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrixC.elements[j][i] = matrix11.elements[j][i];
			matrixC.elements[j + size][i] = matrix21.elements[j][i];
			matrixC.elements[j][i + size] = matrix12.elements[j][i];
			matrixC.elements[j + size][i + size] = matrix12.elements[j][i];
		}
	}

	return matrixC;
}

Matrix MultMatrix(Matrix matrixA, Matrix matrixB)
{
	int size = matrixA.row;
	Matrix matrixC(size, size);
	matrixB.Transpose();

	//? ????? ?????????? ??????? ?? ??????? ??????? C
	for (int j = 0; j < matrixB.row; j++)
	{
		//? ????? ?????????? ??????? ?? ???????? ??????? C
		for (int i = 0; i < matrixA.row; i++)
		{
			matrixC.elements[i][j] = 0;
			//? ????? ?????????? ???????????? ????????? ?????? A ? B ?? ??????????????? ??????
			for (int t = 0; t < matrixA.col; t++)
			{
				matrixC.elements[j][i] += matrixA.elements[i][t] * matrixB.elements[j][t];
			}
		}
	}

	matrixB.Transpose();

	return matrixC;
}

Matrix MultStrassen(Matrix matrixA, Matrix matrixB, int size)
{
	if (size <= 64) {
		return MultMatrix(matrixA, matrixB);
	}

	size = size / 2;

	//x, y
	Position pos11(0, 0);
	//x, y + size
	Position pos12(0, size);
	//x + size, y
	Position pos21(size, 0);
	//x + size, y + size
	Position pos22(size, size);

	Matrix p1 = MultStrassen(SumMatrix(matrixA, pos11, matrixB, pos22, size), SumMatrix(matrixA, pos11, matrixB, pos22, size), size);
	Matrix p2 = MultStrassen(SumMatrix(matrixA, pos21, matrixB, pos22, size), GetPart(matrixB, pos11, size), size);
	Matrix p3 = MultStrassen(GetPart(matrixA, pos11, size), SubMatrix(matrixB, pos12, matrixB, pos22, size), size);
	Matrix p4 = MultStrassen(GetPart(matrixA, pos22, size), SubMatrix(matrixB, pos21, matrixB, pos11, size), size);
	Matrix p5 = MultStrassen(SumMatrix(matrixA, pos11, matrixA, pos12, size), GetPart(matrixB, pos22, size), size);
	Matrix p6 = MultStrassen(SubMatrix(matrixA, pos21, matrixA, pos11, size), SumMatrix(matrixB, pos11, matrixB, pos12, size), size);
	Matrix p7 = MultStrassen(SubMatrix(matrixA, pos12, matrixA, pos22, size), SumMatrix(matrixB, pos21, matrixB, pos22, size), size);

	Matrix c11 = SumTMatrix(SumTMatrix(p1, p4), SubTMatrix(p7, p5));
	Matrix c12 = SumTMatrix(p3, p5);
	Matrix c21 = SumTMatrix(p2, p4);
	Matrix c22 = SumTMatrix(SubTMatrix(p1, p2), SumTMatrix(p3, p6));

	Matrix matrixC = UniteMatrix(c11, c12, c21, c22);

	return matrixC;
}

Matrix SumMatrixParallel(Matrix matrixA, Position posA, Matrix matrixB, Position posB, int size)
{
	Matrix matrixC(size, size);

#pragma omp parallel for
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrixC.elements[j][i] = matrixA.elements[j + posA.y][i + posA.x] + matrixB.elements[j + posB.y][i + posB.x];
		}
	}

	return matrixC;
}

Matrix SubMatrixParallel(Matrix matrixA, Position posA, Matrix matrixB, Position posB, int size)
{
	Matrix matrixC(size, size);

#pragma omp parallel for
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrixC.elements[j][i] = matrixA.elements[j + posA.y][i + posA.x] - matrixB.elements[j + posB.y][i + posB.x];
		}
	}

	return matrixC;
}

Matrix SumTMatrixParallel(Matrix matrixA, Matrix matrixB)
{
	int size = matrixA.col;
	Matrix matrixC(size, size);

#pragma omp parallel for
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrixC.elements[i][j] = matrixA.elements[j][i] + matrixB.elements[j][i];
		}
	}

	return matrixC;
}

Matrix SubTMatrixParallel(Matrix matrixA, Matrix matrixB)
{
	int size = matrixA.col;
	Matrix matrixC(size, size);
#pragma omp parallel for
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrixC.elements[j][i] = matrixA.elements[j][i] - matrixB.elements[j][i];
		}
	}

	return matrixC;
}

Matrix GetPartParallel(Matrix matrixA, Position posA, int size)
{
	Matrix matrixC(size, size);
#pragma omp parallel for
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrixC.elements[j][i] = matrixA.elements[j + posA.y][i + posA.x];
		}
	}

	return matrixC;
}

Matrix UniteMatrixParallel(Matrix matrix11, Matrix matrix12, Matrix matrix21, Matrix matrix22)
{
	int size = matrix11.col;
	Matrix matrixC(size * 2, size * 2);

#pragma omp parallel for
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrixC.elements[j][i] = matrix11.elements[j][i];
			matrixC.elements[j + size][i] = matrix21.elements[j][i];
			matrixC.elements[j][i + size] = matrix12.elements[j][i];
			matrixC.elements[j + size][i + size] = matrix12.elements[j][i];
		}
	}

	return matrixC;
}

Matrix MultMatrixParallel(Matrix matrixA, Matrix matrixB)
{
	int size = matrixA.row;
	Matrix matrixC(size, size);
	matrixB.Transpose();

	//? ????? ?????????? ??????? ?? ??????? ??????? C
#pragma omp parallel for
	for (int j = 0; j < matrixB.row; j++)
	{
		//? ????? ?????????? ??????? ?? ???????? ??????? C
		for (int i = 0; i < matrixA.row; i++)
		{
			matrixC.elements[i][j] = 0;
			//? ????? ?????????? ???????????? ????????? ?????? A ? B ?? ??????????????? ??????
			for (int t = 0; t < matrixA.col; t++)
			{
				matrixC.elements[j][i] += matrixA.elements[i][t] * matrixB.elements[j][t];
			}
		}
	}

	matrixB.Transpose();

	return matrixC;
}

Matrix MultStrassenParallel(Matrix matrixA, Matrix matrixB, int size)
{
	if (size <= 64) {
		return MultMatrixParallel(matrixA, matrixB);
	}

	size = size / 2;

	//x, y
	Position pos11(0, 0);
	//x, y + size
	Position pos12(0, size);
	//x + size, y
	Position pos21(size, 0);
	//x + size, y + size
	Position pos22(size, size);

	Matrix p1 = MultStrassenParallel(SumMatrixParallel(matrixA, pos11, matrixB, pos22, size), SumMatrixParallel(matrixA, pos11, matrixB, pos22, size), size);
	Matrix p2 = MultStrassenParallel(SumMatrixParallel(matrixA, pos21, matrixB, pos22, size), GetPartParallel(matrixB, pos11, size), size);
	Matrix p3 = MultStrassenParallel(GetPartParallel(matrixA, pos11, size), SubMatrixParallel(matrixB, pos12, matrixB, pos22, size), size);
	Matrix p4 = MultStrassenParallel(GetPartParallel(matrixA, pos22, size), SubMatrixParallel(matrixB, pos21, matrixB, pos11, size), size);
	Matrix p5 = MultStrassenParallel(SumMatrixParallel(matrixA, pos11, matrixA, pos12, size), GetPartParallel(matrixB, pos22, size), size);
	Matrix p6 = MultStrassenParallel(SubMatrixParallel(matrixA, pos21, matrixA, pos11, size), SumMatrixParallel(matrixB, pos11, matrixB, pos12, size), size);
	Matrix p7 = MultStrassenParallel(SubMatrixParallel(matrixA, pos12, matrixA, pos22, size), SumMatrixParallel(matrixB, pos21, matrixB, pos22, size), size);

	Matrix c11 = SumTMatrixParallel(SumTMatrixParallel(p1, p4), SubTMatrixParallel(p7, p5));
	Matrix c12 = SumTMatrixParallel(p3, p5);
	Matrix c21 = SumTMatrixParallel(p2, p4);
	Matrix c22 = SumTMatrixParallel(SubTMatrixParallel(p1, p2), SumTMatrixParallel(p3, p6));

	Matrix matrixC = UniteMatrixParallel(c11, c12, c21, c22);

	return matrixC;
}

void MultiplicationStrassenConsistently(Matrix matrixA, Matrix matrixB, Matrix matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();

	int m = matrixA.col;
	int n = matrixA.row;
	int k = matrixB.col;
	int size = GetDimentionDegree2(matrixA.col, matrixA.row, matrixB.col);
	Matrix* copyMatrixA = MatrixExtension(matrixA, size);
	Matrix* copyMatrixB = MatrixExtension(matrixB, size);
	Matrix* copyMatrixC = MatrixExtension(matrixC, size);

	copyMatrixC = &MultStrassen(*copyMatrixA, *copyMatrixB, size);

	endTime = omp_get_wtime();
	delete copyMatrixA, copyMatrixB, copyMatrixC;
	time = endTime - startTime;
}

void MultiplicationStrassenParallel(Matrix matrixA, Matrix matrixB, Matrix matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();

	int m = matrixA.col;
	int n = matrixA.row;
	int k = matrixB.col;
	int size = GetDimentionDegree2(matrixA.col, matrixA.row, matrixB.col);
	//Matrix copyMatrixC;

	//copyMatrixC = MultStrassenParallel(MatrixExtension(matrixA, size), MatrixExtension(matrixB, size), size);
	//copyMatrixC.SetToNull();
	Matrix* copyMatrixA = MatrixExtension(matrixA, size);
	Matrix* copyMatrixB = MatrixExtension(matrixB, size);
	Matrix* copyMatrixC = MatrixExtension(matrixC, size);
	
	copyMatrixC = &MultStrassenParallel(*copyMatrixA, *copyMatrixB, size);

	endTime = omp_get_wtime();
	delete copyMatrixA, copyMatrixB, copyMatrixC;
	time = endTime - startTime;
}

void CalculateMultFuncTrustedTime(void* func, int m, int n, int k, double& time, int iterations)
{
	Matrix matrixA(n, m);
	Matrix matrixB(k, n);

	//?????????? ??????
	FillLineConsistently(matrixA, time);
	FillLineConsistently(matrixB, time);

	Matrix matrixC(m, k);

	double curtime = 0, avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	double* Times = new double[iterations];

	cout << "[";
	for (int i = 0; i < iterations; i++)
	{
		(*(MultFunctTempl)func)(matrixA, matrixB, matrixC, curtime);
		Times[i] = curtime;
		avgTime += curtime;
		cout << "*";
	}
	cout << "]\n";
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	time = avgTimeT;
}

void CalculateAllFillFuncs(int m, int n)
{
	double time;
	string* funcsNames = new string[6]{ "????????? ???????????????? ??????????", "????????? ???????????? ?????????? ? for static", "????????? ???????????? ?????????? ? for dynamic", "????????? ???????????? ? sections", "??????? ????????????????", "??????? ???????????? ?????????? ? for static" };
	void** funcs = new void*[6]{ FillLineConsistently, FillLineParallelForStatic, FillLineParallelForDynamic, FillLineParallelSections, FillBlockConsistently, FillBlockParallelForStatic };

	cout << "??:" << n << "x" << m << endl;
	for (int t = 2; t < 5; t++)
	{
		omp_set_num_threads(t);
		cout << "????? ???????:" << t << endl;
		for (int i = 0; i < 6; i++)
		{
			time = 0;
			CalculateFillFuncTrustedTime(funcs[i], m, n, time, 20);
			cout << funcsNames[i] << "??????? ?????????: " << time * 1000 << endl;
		}
		cout << endl;
	}
}

void CalculateAllMultFuncs(int m, int n, int k)
{
	double time;
	string* funcsNames = new string[6]{ "???????????????? ?????????", "???????????? ????????? ? for static", "???????????? ????????? ? for dynamic", "???????????? ????????? ? sections", "???????????????? ??????? ????????? (????????)", "???????????? ??????? ????????? (????????)"  };
	void** funcs = new void*[6]{ MultiplicationConsistently, MultiplicationParallelForStatic, MultiplicationParallelForDynamic, MultiplicationParallelSections, MultiplicationStrassenConsistently, MultiplicationStrassenParallel };

	cout << "??:" << n << "x" << m << " ? " << k << "x" << n << endl;
	for (int t = 2; t < 5; t++)
	{
		omp_set_num_threads(t);
		cout << "????? ???????:" << t << endl;
		for (int i = 0; i < 6; i++)
		{
			time = 0;
			CalculateMultFuncTrustedTime(funcs[i], m, n, k, time, 20);
			cout << funcsNames[i] << "??????? ?????????: " << time * 1000 << endl;
		}
		cout << endl;
	}
}

void TaskFillAndMultiplyMatrix()
{
	int* n = new int[4]{ 850, 1250, 1500, 2000 };
	int* m = new int[4]{ 700, 1050, 1500, 2000 };
	int* k = new int[4]{ 650, 1000, 1700, 1900 };

	double time;
	string* funcsNamesFill = new string[6]{ "????????? ???????????????? ??????????", "????????? ???????????? ?????????? ? for static", "????????? ???????????? ?????????? ? for dynamic", "????????? ???????????? ? sections", "??????? ????????????????", "??????? ???????????? ?????????? ? for static" };
	void** funcsFill = new void*[6]{ FillLineConsistently, FillLineParallelForStatic, FillLineParallelForDynamic, FillLineParallelSections, FillBlockConsistently, FillBlockParallelForStatic };
	string* funcsNamesMult = new string[6]{ "???????????????? ?????????", "???????????? ????????? ? for static", "???????????? ????????? ? for dynamic", "???????????? ????????? ? sections", "???????????????? ??????? ????????? (????????)", "???????????? ??????? ????????? (????????)" };
	void** funcsMult = new void*[6]{ MultiplicationConsistently, MultiplicationParallelForStatic, MultiplicationParallelForDynamic, MultiplicationParallelSections, MultiplicationStrassenConsistently, MultiplicationStrassenParallel };

	std::ofstream resultsFile;
	resultsFile.open("Task3_Results.csv", std::ios_base::app);
	resultsFile.clear();
	for (int d = 0; d < 4; d++)
	{
		resultsFile << "?? A = " << n[d] << "x" << m[d] << "B = " << k[d] << "x" << n[d] << endl;
		resultsFile << "???????;??????;?????;Sp(n);Ep(n);\n";
		cout << "?? A = " << n[d] << "x" << m[d] << "B = " << k[d] << "x" << n[d] << endl;
		//??????? ?????????? ???????
		for (int i = 0; i < 6; i++)
		{
			resultsFile << funcsNamesFill[i] << ";";
			cout << funcsNamesFill[i] << ":" << endl;
			for (int t = 2; t < 5; t++)
			{
				if (t != 2)
					resultsFile << ";";

				omp_set_num_threads(t);
				cout << "????? ???????:" << t << endl;
				time = 0;
				CalculateFillFuncTrustedTime(funcsFill[i], m[d], n[d], time, 20);
				
				cout << "??????? ?????????: " << time * 1000 << endl;
				if (i == 0 || i == 4)
				{
					resultsFile << ";" << time * 1000 << "^;\n";
					break;
				}
				else
				{
					resultsFile << t << ";" << time * 1000 << "^;\n";
				}
			}
			cout << endl;
		}

		//??????? ???????????? ??????
		for (int i = 0; i < 6; i++)
		{
			resultsFile << funcsNamesMult[i] << ";";
			cout << funcsNamesMult[i] << ":" << endl;
			for (int t = 2; t < 5; t++)
			{
				if (t != 2)
					resultsFile << ";";

				omp_set_num_threads(t);
				cout << "????? ???????:" << t << endl;
				time = 0;
				CalculateMultFuncTrustedTime(funcsMult[i], m[d], n[d], k[d], time, 20);

				cout << "??????? ?????????: " << time * 1000 << endl;
				if (i == 0 || i == 4)
				{
					resultsFile << ";" << time * 1000 << "^;\n";
					break;
				}
				else
				{
					resultsFile << t << ";" << time * 1000 << "^;\n";
				}
			}
			cout << endl;
		}

		resultsFile << endl;
	}

	resultsFile.close();
}

void main()
{
	srand(std::time(0));
	setlocale(LC_ALL, "Russian");
	int choice;

	//double time;
	//CalculateMultFuncTrustedTime(MultiplicationStrassenParallel, 1000, 1000, 1000, time, 10);
	//cout << ". ??????? ?????????: " << time * 1000 << endl;

	cout << "[1]: ?????????? ????? ??\n"
		<< "[2]: ???????????. ?????????? ? ???????????? ??????\n";
	cin >> choice;

	if (choice == 1)
	{
		TaskCalculatePi();
	}
	else
	{
		TaskFillAndMultiplyMatrix();
	}

	cout << endl << "?????? ????????? ?????????" << endl;
	string wait;
	cin >> wait;
}