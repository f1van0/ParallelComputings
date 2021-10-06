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
		elements = new double*[col];
		for (int i = 0; i < col; i++)
		{
			elements[i] = new double[row];
		}

		SetToNull();
		this->col = col;
		this->row = row;
	}

	Matrix(const Matrix& _matrix)
	{
		col = _matrix.col;
		row = _matrix.row;
		elements = new double*[col];
		for (int j = 0; j < col; j++)
		{
			elements[j] = new double[row];
		}

		for (int j = 0; j < row; j++)
		{
			for (int i = 0; i < col; i++)
			{
				elements[i][j] = _matrix.elements[i][j];
			}
		}
	}

	~Matrix()
	{
		for (int i = 0; i < col; i++)
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
				elements[i][j] = 0;
			}
		}
	}

	//Функция транспонирования матрицы
	void Transpose()
	{
		double** tElements;
		tElements = new double*[row];
		for (int i = 0; i < row; i++)
		{
			tElements[i] = new double[col];
			for (int j = 0; j < col; j++)
			{
				tElements[i][j] = elements[j][i];
			}
		}

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
			cout << setw(5) << matrix.elements[i][j] << " ";
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
	int i;
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
#pragma omp sections private(i)
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
	time = 0;
	return 3.14;
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

void Task1()
{
	string* funcsNames = new string[6]{ "Последовательная реализация", "Параллельная реализация FOR (static)", "Параллельная реализация FOR (dynamic)", "Параллельная реализация FOR (guided)", "Параллельная реализация Sections", "Свой алгоритм" };
	void** funcs = new void*[6]{ CalculatePiConsistenty, CalculatePiParallelForStatic, CalculatePiParallelForDynamic, CalculatePiParallelForGuided, CalculatePiParallelSections, CalculatePiMyAlgorithm };

	static long* nums_steps = new long[4]{ 1000000, 5000000, 10000000, 15000000 };
	double time, pi;
	for (int d = 0; d < 4; d++)
	{
		cout << "НД:" << nums_steps[d] << endl;
		for (int t = 2; t < 5; t++)
		{
			omp_set_num_threads(t);
			cout << "Число потоков:" << t << endl;
			for (int i = 0; i < 6; i++)
			{
				time = 0;
				pi = CalculatePiFuncTrustedTime(funcs[i], nums_steps[d], time, 50);
				cout << funcsNames[i] << ". Число Пи = " << pi << ". Времени затрачено: " << time * 1000 << endl;
			}
			cout << endl;
		}
	}
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
			matrix.elements[i][j] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
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
			matrix.elements[i][j] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
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
			matrix.elements[i][j] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
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
					matrix.elements[i][j] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
				}
			}
		}
#pragma omp section
		{
			for (int j = s1; j < s2; j++)
			{
				for (int i = 0; i < matrix.col; i++)
				{
					matrix.elements[i][j] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
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
						matrix.elements[i][j] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
					}
				}
			}
		}
#pragma omp section
		{
			if (threadsNum > 2)
			{
				for (int j = s3; j < matrix.row; j++)
				{
					for (int i = 0; i < matrix.col; i++)
					{
						matrix.elements[i][j] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
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
			for (int d = 0; d < blockSizeX; d++)
			{
				int dIndex = d + shiftX;
				int tIndex = t + shiftY;
				matrix.elements[dIndex][tIndex] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1);
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
			for (int d = 0; d < blockSizeX; d++)
			{
				int dIndex = d + shiftX;
				int tIndex = t + shiftY;
				matrix.elements[dIndex][tIndex] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1);
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
			for (int d = 0; d < blockSizeX; d++)
			{
				int dIndex = d + shiftX;
				int tIndex = t + shiftY;
				matrix.elements[dIndex][tIndex] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1);
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

typedef double(*MultFunctTempl)(Matrix, Matrix, Matrix&, double&);

void MultiplicationConsistently(Matrix matrixA, Matrix matrixB, Matrix& matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	matrixB.Transpose();

	//В цикле происходит переход по строкам матрицы C
	for (int j = 0; j < matrixB.row; j++)
	{
		//В цикле происходит переход по столбцам матрицы C
		for (int i = 0; i < matrixA.row; i++)
		{
			matrixC.elements[i][j] = 0;
			//В цикле происходит перемножение элементов матриц A и B на соответствующей строке
			for (int t = 0; t < matrixA.col; t++)
			{
				matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void MultiplicationParallelForStatic(Matrix matrixA, Matrix matrixB, Matrix& matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	matrixB.Transpose();

#pragma omp parallel for schedule(static, matrixA.row / 12)
	//В цикле происходит переход по строкам матрицы C
	for (int j = 0; j < matrixB.row; j++)
	{
		//В цикле происходит переход по столбцам матрицы C
		for (int i = 0; i < matrixA.row; i++)
		{
			matrixC.elements[i][j] = 0;
			//В цикле происходит перемножение элементов матриц A и B на соответствующей строке
			for (int t = 0; t < matrixA.col; t++)
			{
				matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void MultiplicationParallelForDynamic(Matrix matrixA, Matrix matrixB, Matrix& matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	matrixB.Transpose();

#pragma omp parallel for schedule(dynamic, matrixA.row / 6)
	//В цикле происходит переход по строкам матрицы C
	for (int j = 0; j < matrixB.row; j++)
	{
		//В цикле происходит переход по столбцам матрицы C
		for (int i = 0; i < matrixA.row; i++)
		{
			matrixC.elements[i][j] = 0;
			//В цикле происходит перемножение элементов матриц A и B на соответствующей строке
			for (int t = 0; t < matrixA.col; t++)
			{
				matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void MultiplicationParallelSections(Matrix matrixA, Matrix matrixB, Matrix& matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	matrixB.Transpose();

	int numThreads = 0;
	int s1, s2, s3;

#pragma omp parallel
	{
		numThreads = omp_get_max_threads();
	}
	
	s1 = matrixB.row / numThreads;
	s2 = matrixB.row * 2 / numThreads;
	s3 = matrixB.row * 3 / numThreads;

#pragma omp sections
	{
#pragma omp section
		{
			for (int j = 0; j < s1; j++)
			{
				//В цикле происходит переход по столбцам матрицы C
				for (int i = 0; i < matrixA.row; i++)
				{
					matrixC.elements[i][j] = 0;
					//В цикле происходит перемножение элементов матриц A и B на соответствующей строке
					for (int t = 0; t < matrixA.col; t++)
					{
						matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
					}
				}
			}
		}
#pragma omp section
		{
			for (int j = s1; j < s2; j++)
			{
				//В цикле происходит переход по столбцам матрицы C
				for (int i = 0; i < matrixA.row; i++)
				{
					matrixC.elements[i][j] = 0;
					//В цикле происходит перемножение элементов матриц A и B на соответствующей строке
					for (int t = 0; t < matrixA.col; t++)
					{
						matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
					}
				}
			}
		}
#pragma omp section
		{
			if (numThreads > 2)
				for (int j = s2; j < s3; j++)
				{
					//В цикле происходит переход по столбцам матрицы C
					for (int i = 0; i < matrixA.row; i++)
					{
						matrixC.elements[i][j] = 0;
						//В цикле происходит перемножение элементов матриц A и B на соответствующей строке
						for (int t = 0; t < matrixA.col; t++)
						{
							matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
						}
					}
				}
		}
#pragma omp section
		{
			if (numThreads > 3)
				for (int j = s3; j < matrixB.row; j++)
				{
					//В цикле происходит переход по столбцам матрицы C
					for (int i = 0; i < matrixA.row; i++)
					{
						matrixC.elements[i][j] = 0;
						//В цикле происходит перемножение элементов матриц A и B на соответствующей строке
						for (int t = 0; t < matrixA.col; t++)
						{
							matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
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

Matrix& MatrixExtension(Matrix matrix, int size) 
{
	Matrix* result = new Matrix(size, size);

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i < matrix.col && j < matrix.row)
				//Ошибка при увеличении размера MatrixC на 4 потоках 1 итерация
				result->elements[i][j] = matrix.elements[i][j];
			else
				result->elements[i][j] = 0;
		}
	}

	return *result;
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
			matrixC.elements[i][j] = matrixA.elements[i + posA.x][j + posA.y] + matrixB.elements[i + posB.x][j + posB.y];
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
			matrixC.elements[i][j] = matrixA.elements[i + posA.x][j + posA.y] - matrixB.elements[i + posB.x][j + posB.y];
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
			matrixC.elements[i][j] = matrixA.elements[i][j] + matrixB.elements[i][j];
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
			matrixC.elements[i][j] = matrixA.elements[i][j] - matrixB.elements[i][j];
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
			matrixC.elements[i][j] = matrixA.elements[i + posA.x][j + posA.y];
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
			matrixC.elements[i][j] = matrix11.elements[i][j];
			matrixC.elements[i + size][j] = matrix21.elements[i][j];
			matrixC.elements[i][j + size] = matrix12.elements[i][j];
			matrixC.elements[i + size][j + size] = matrix12.elements[i][j];
		}
	}

	return matrixC;
}

Matrix MultMatrix(Matrix matrixA, Matrix matrixB)
{
	int size = matrixA.row;
	Matrix matrixC(size, size);
	matrixB.Transpose();

	//В цикле происходит переход по строкам матрицы C
	for (int j = 0; j < matrixB.row; j++)
	{
		//В цикле происходит переход по столбцам матрицы C
		for (int i = 0; i < matrixA.row; i++)
		{
			matrixC.elements[i][j] = 0;
			//В цикле происходит перемножение элементов матриц A и B на соответствующей строке
			for (int t = 0; t < matrixA.col; t++)
			{
				matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
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
			matrixC.elements[i][j] = matrixA.elements[i + posA.x][j + posA.y] + matrixB.elements[i + posB.x][j + posB.y];
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
			matrixC.elements[i][j] = matrixA.elements[i + posA.x][j + posA.y] - matrixB.elements[i + posB.x][j + posB.y];
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
			matrixC.elements[i][j] = matrixA.elements[i][j] + matrixB.elements[i][j];
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
			matrixC.elements[i][j] = matrixA.elements[i][j] - matrixB.elements[i][j];
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
			matrixC.elements[i][j] = matrixA.elements[i + posA.x][j + posA.y];
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
			matrixC.elements[i][j] = matrix11.elements[i][j];
			matrixC.elements[i + size][j] = matrix21.elements[i][j];
			matrixC.elements[i][j + size] = matrix12.elements[i][j];
			matrixC.elements[i + size][j + size] = matrix12.elements[i][j];
		}
	}

	return matrixC;
}

Matrix MultMatrixParallel(Matrix matrixA, Matrix matrixB)
{
	int size = matrixA.row;
	Matrix matrixC(size, size);
	matrixB.Transpose();

	//В цикле происходит переход по строкам матрицы C
#pragma omp parallel for
	for (int j = 0; j < matrixB.row; j++)
	{
		//В цикле происходит переход по столбцам матрицы C
		for (int i = 0; i < matrixA.row; i++)
		{
			matrixC.elements[i][j] = 0;
			//В цикле происходит перемножение элементов матриц A и B на соответствующей строке
			for (int t = 0; t < matrixA.col; t++)
			{
				matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
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
/*
Matrix MultStrassenParallel(Matrix matrixA, Matrix matrixB, int _size)
{
	if (_size <= 64) {
		return MultMatrix(matrixA, matrixB);
	}

	_size = _size / 2;

	//x, y
	Position pos11(0, 0);
	//x, y + size
	Position pos12(0, _size);
	//x + size, y
	Position pos21(_size, 0);
	//x + size, y + size
	Position pos22(_size, _size);

	Matrix p1, p2, p3, p4, p5, p6, p7;

#pragma omp sections private(matrixA, matrixB)
	{
#pragma omp section
		{
			p1 = MultStrassen(SumMatrix(matrixA, pos11, matrixB, pos22, _size), SumMatrix(matrixA, pos11, matrixB, pos22, _size), _size);
		}
#pragma omp section
		{
			p2 = MultStrassen(SumMatrix(matrixA, pos21, matrixB, pos22, _size), GetPart(matrixB, pos11, _size), _size);
		}
#pragma omp section
		{
			p3 = MultStrassen(GetPart(matrixA, pos11, _size), SubMatrix(matrixB, pos12, matrixB, pos22, _size), _size);
		}
#pragma omp section
		{
			p4 = MultStrassen(GetPart(matrixA, pos22, _size), SubMatrix(matrixB, pos21, matrixB, pos11, _size), _size);
		}
#pragma omp section
		{
			p5 = MultStrassen(SumMatrix(matrixA, pos11, matrixA, pos12, _size), GetPart(matrixB, pos22, _size), _size);
		}
#pragma omp section
		{
			p6 = MultStrassen(SubMatrix(matrixA, pos21, matrixA, pos11, _size), SumMatrix(matrixB, pos11, matrixB, pos12, _size), _size);
		}
#pragma omp section
		{
			p7 = MultStrassen(SubMatrix(matrixA, pos12, matrixA, pos22, _size), SumMatrix(matrixB, pos21, matrixB, pos22, _size), _size);
		}
	}

	Matrix c11 = SumTMatrix(SumTMatrix(p1, p4), SubTMatrix(p7, p5));
	Matrix c12 = SumTMatrix(p3, p5);
	Matrix c21 = SumTMatrix(p2, p4);
	Matrix c22 = SumTMatrix(SubTMatrix(p1, p2), SumTMatrix(p3, p6));

	Matrix matrixC = UniteMatrix(c11, c12, c21, c22);

	return matrixC;
}
*/

void MultiplicationStrassenConsistently(Matrix matrixA, Matrix matrixB, Matrix& matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();

	int m = matrixA.col;
	int n = matrixA.row;
	int k = matrixB.col;
	int size = GetDimentionDegree2(matrixA.col, matrixA.row, matrixB.col);
	Matrix* copyMatrixA = &MatrixExtension(matrixA, size);
	Matrix* copyMatrixB = &MatrixExtension(matrixB, size);
	Matrix* copyMatrixC = &MatrixExtension(matrixC, size);

	copyMatrixC = &MultStrassen(*copyMatrixA, *copyMatrixB, size);

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void MultiplicationStrassenParallel(Matrix matrixA, Matrix matrixB, Matrix& matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();

	int m = matrixA.col;
	int n = matrixA.row;
	int k = matrixB.col;
	int size = GetDimentionDegree2(matrixA.col, matrixA.row, matrixB.col);
	Matrix* copyMatrixA = &MatrixExtension(matrixA, size);
	Matrix* copyMatrixB = &MatrixExtension(matrixB, size);
	Matrix* copyMatrixC = &MatrixExtension(matrixC, size);

	copyMatrixC = &MultStrassenParallel(*copyMatrixA, *copyMatrixB, size);

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void CalculateMultFuncTrustedTime(void* func, int m, int n, int k, double& time, int iterations)
{
	Matrix matrixA(n, m);
	Matrix matrixB(k, n);

	//Заполнение матриц
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
	string* funcsNames = new string[6]{ "Ленточное последовательное заполнение", "Ленточное параллельное заполнение с for static", "Ленточное параллельное заполнение с for dynamic", "Ленточное параллельный с sections", "Блочное последовательное", "Блочное параллельное заполнение с for static" };
	void** funcs = new void*[6]{ FillLineConsistently, FillLineParallelForStatic, FillLineParallelForDynamic, FillLineParallelSections, FillBlockConsistently, FillBlockParallelForStatic };

	cout << "НД:" << n << "x" << m << endl;
	for (int t = 2; t < 5; t++)
	{
		omp_set_num_threads(t);
		cout << "Число потоков:" << t << endl;
		for (int i = 0; i < 6; i++)
		{
			time = 0;
			CalculateFillFuncTrustedTime(funcs[i], m, n, time, 20);
			cout << funcsNames[i] << ". Времени затрачено: " << time * 1000 << endl;
		}
		cout << endl;
	}
}

void CalculateAllMultFuncs(int m, int n, int k)
{
	double time;
	string* funcsNames = new string[6]{ "Последовательное умножение", "Параллельное умножение с for static", "Параллельное умножение с for dynamic", "Параллельное умножение с sections", "Последовательное быстрое умножение (Штрассен)", "Параллельное быстрое умножение (Штрассен)"  };
	void** funcs = new void*[6]{ MultiplicationConsistently, MultiplicationParallelForStatic, MultiplicationParallelForDynamic, MultiplicationParallelSections, MultiplicationStrassenConsistently, MultiplicationStrassenParallel };

	cout << "НД:" << n << "x" << m << " и " << k << "x" << n << endl;
	for (int t = 2; t < 5; t++)
	{
		omp_set_num_threads(t);
		cout << "Число потоков:" << t << endl;
		for (int i = 0; i < 6; i++)
		{
			time = 0;
			CalculateMultFuncTrustedTime(funcs[i], m, n, k, time, 30);
			cout << funcsNames[i] << ". Времени затрачено: " << time * 1000 << endl;
		}
		cout << endl;
	}
}

void Task2()
{
	int n = 1200;
	int m = 1300;
	int k = 1500;
	Matrix* matrixA = new Matrix(n, m);
	Matrix* matrixB = new Matrix(k, n);
	Matrix* matrixC = new Matrix(k, m);
	double time = 0;

	FillLineConsistently(*matrixA, time);
	FillLineConsistently(*matrixB, time);

	string* funcsNames = new string[1]{ "Последовательное перемножение" };
	void** funcs = new void*[1]{ MultiplicationConsistently };

	for (int d = 0; d < 4; d++)
	{
		for (int t = 2; t < 5; t++)
		{
			omp_set_num_threads(t);
			cout << "Число потоков:" << t << endl;
			for (int i = 0; i < 5; i++)
			{
				time = 0;
				CalculateMultFuncTrustedTime(funcs[i], n, m, k, time, 30);
				cout << funcsNames[i] << ". Времени затрачено: " << time * 1000 << endl;
			}
			cout << endl;
		}
	}
}

void Task3()
{
	int* n = new int[4]{ 85, 1750, 1500, 2000 };
	int* m = new int[4]{ 70, 1050, 1500, 2000 };
	int* k = new int[4]{ 65, 1000, 2000, 1900 };
	
	for (int d = 0; d < 4; d++)
	{
		//CalculateAllFillFuncs(m[d], n[d]);
		CalculateAllMultFuncs(m[d], n[d], k[d]);
	}
}

/*
void FillTest(Matrix& matrix)
{
	for (int j = 0; j < matrix.n; j++)
	{
		for (int i = 0; i < matrix.m; i++)
		{
			matrix.elements[i][j] = rand() % 10;
		}
	}
}
*/

void main()
{
	srand(std::time(0));
	setlocale(LC_ALL, "Russian");
	int choice;

	cout << "[1]: Задание 1\n"
		<< "[2]: Задание 2\n"
		<< "[3]: Задание 3\n";
	cin >> choice;

	if (choice == 1)
	{
		Task1();
	}
	else if (choice == 2)
	{
		Task2();
	}
	else
	{
		Task3();
	}

	cout << endl << "Работа программы завершена" << endl;
	string wait;
}