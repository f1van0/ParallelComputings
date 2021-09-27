#include <omp.h>
#include <stdio.h>
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
	int m;
	int n;

	Matrix()
	{
		elements = nullptr;
		m = 0;
		n = 0;
	}

	Matrix(int m, int n)
	{
		elements = new double*[m];
		for (int i = 0; i < m; i++)
		{
			elements[i] = new double[n];
		}

		this->m = m;
		this->n = n;
	}

	~Matrix()
	{
		for (int i = 0; i < m; i++)
		{
			delete[] elements[i];
		}

		delete[] elements;
	}

	void SetToNull()
	{
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < m; i++)
			{
				elements[i][j] = 0;
			}
		}
	}

	//������� ���������������� �������
	void Transpose()
	{
		double** tElements;
		tElements = new double*[n];
		for (int i = 0; i < n; i++)
		{
			tElements[i] = new double[m];
			for (int j = 0; j < m; j++)
			{
				tElements[i][j] = elements[j][i];
			}
		}

		int c = 0;
		c = m;
		m = n;
		n = c;

		elements = tElements;
	}
};

void PrintMatrix(Matrix matrix)
{
	for (int j = 0; j < matrix.n; j++)
	{
		for (int i = 0; i < matrix.m; i++)
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
	string* funcsNames = new string[6]{ "���������������� ����������", "������������ ���������� FOR (static)", "������������ ���������� FOR (dynamic)", "������������ ���������� FOR (guided)", "������������ ���������� Sections", "���� ��������" };
	void** funcs = new void*[6]{ CalculatePiConsistenty, CalculatePiParallelForStatic, CalculatePiParallelForDynamic, CalculatePiParallelForGuided, CalculatePiParallelSections, CalculatePiMyAlgorithm };

	static long* nums_steps = new long[4]{ 1000000, 5000000, 10000000, 15000000 };
	double time, pi;
	for (int d = 0; d < 4; d++)
	{
		cout << "��:" << nums_steps[d] << endl;
		for (int t = 2; t < 5; t++)
		{
			omp_set_num_threads(t);
			cout << "����� �������:" << t << endl;
			for (int i = 0; i < 6; i++)
			{
				time = 0;
				pi = CalculatePiFuncTrustedTime(funcs[i], nums_steps[d], time, 50);
				cout << funcsNames[i] << ". ����� �� = " << pi << ". ������� ���������: " << time * 1000 << endl;
			}
			cout << endl;
		}
	}
}

typedef double(*FillFunctTempl)(Matrix&, double&);

void FillLineConsistently(Matrix& matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	for (int j = 0; j < matrix.n; j++)
	{
		for (int i = 0; i < matrix.m; i++)
		{
			matrix.elements[i][j] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
		}
	}
	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void FillLineParallelForStatic(Matrix& matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
#pragma omp parallel for schedule(static, matrix.n / 12)
	for (int j = 0; j < matrix.n; j++)
	{
		for (int i = 0; i < matrix.m; i++)
		{
			matrix.elements[i][j] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
		}
	}
	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void FillLineParallelForDynamic(Matrix& matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, matrix.n / 6)
	for (int j = 0; j < matrix.n; j++)
	{
		for (int i = 0; i < matrix.m; i++)
		{
			matrix.elements[i][j] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
		}
	}
	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void FillLineParallelSections(Matrix& matrix, double& time)
{
	int threadsNum;
#pragma omp parallel
	{
		threadsNum = omp_get_num_threads();
	}
	double startTime, endTime;
	int s1 = matrix.n / threadsNum, s2 = matrix.n * 2 / threadsNum, s3 = matrix.n * 3 / threadsNum;

	if (threadsNum == 3)
		s3 = matrix.n;

	startTime = omp_get_wtime();
#pragma omp sections
	{
#pragma omp section
		{
			for (int j = 0; j < s1; j++)
			{
				for (int i = 0; i < matrix.m; i++)
				{
					matrix.elements[i][j] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
				}
			}
		}
#pragma omp section
		{
			for (int j = s1; j < s2; j++)
			{
				for (int i = 0; i < matrix.m; i++)
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
					for (int i = 0; i < matrix.m; i++)
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
				for (int j = s3; j < matrix.n; j++)
				{
					for (int i = 0; i < matrix.m; i++)
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

void FillBlockParallelForStatic(Matrix& matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	int blockSizeX = 5;
	int blockSizeY = 5;
	int blocksAmount = matrix.n * matrix.m / (blockSizeX * blockSizeY);
	int shiftX = 0, shiftY = 0;
	int m = matrix.m;

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

void FillBlockParallelForDynamic(Matrix& matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	int blockSizeX = 5;
	int blockSizeY = 5;
	int blocksAmount = matrix.n * matrix.m / (blockSizeX * blockSizeY);
	int shiftX = 0, shiftY = 0;
	int m = matrix.m;

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
	Matrix* matrixA = new Matrix(m, n);
	double curtime = 0, avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	double* Times = new double[iterations];

	for (int i = 0; i < iterations; i++)
	{
		(*(FillFunctTempl)func)(*matrixA, curtime);
		Times[i] = curtime;
		avgTime += curtime;
	}
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	time = avgTimeT;
}

typedef double(*MultFunctTempl)(Matrix&, Matrix&, Matrix&, double&);

void MultiplicationConsistently(Matrix& matrixA, Matrix& matrixB, Matrix& matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();

	//� ����� ���������� ������� �� ������� ������� C
	for (int j = 0; j < matrixB.n; j++)
	{
		//� ����� ���������� ������� �� �������� ������� C
		for (int i = 0; i < matrixA.n; i++)
		{
			matrixC.elements[i][j] = 0;
			//� ����� ���������� ������������ ��������� ������ A � B �� ��������������� ������
			for (int t = 0; t < matrixA.m; t++)
			{
				matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void MultiplicationParallelForStatic(Matrix& matrixA, Matrix& matrixB, Matrix& matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();

#pragma omp parallel for schedule(static, matrixA.n / 12)
	//� ����� ���������� ������� �� ������� ������� C
	for (int j = 0; j < matrixB.n; j++)
	{
		//� ����� ���������� ������� �� �������� ������� C
		for (int i = 0; i < matrixA.n; i++)
		{
			matrixC.elements[i][j] = 0;
			//� ����� ���������� ������������ ��������� ������ A � B �� ��������������� ������
			for (int t = 0; t < matrixA.m; t++)
			{
				matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void MultiplicationParallelForDynamic(Matrix& matrixA, Matrix& matrixB, Matrix& matrixC, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();

#pragma omp parallel for schedule(dynamic, matrixA.n / 6)
	//� ����� ���������� ������� �� ������� ������� C
	for (int j = 0; j < matrixB.n; j++)
	{
		//� ����� ���������� ������� �� �������� ������� C
		for (int i = 0; i < matrixA.n; i++)
		{
			matrixC.elements[i][j] = 0;
			//� ����� ���������� ������������ ��������� ������ A � B �� ��������������� ������
			for (int t = 0; t < matrixA.m; t++)
			{
				matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void MultiplicationParallelSections(Matrix& matrixA, Matrix& matrixB, Matrix& matrixC, double& time)
{
	double startTime, endTime;
	int numThreads = 0;
	int s1, s2, s3;

#pragma omp parallel
	{
		numThreads = omp_get_max_threads();
	}
	
	s1 = matrixB.n / numThreads;
	s2 = matrixB.n * 2 / numThreads;
	s3 = matrixB.n * 3 / numThreads;

	startTime = omp_get_wtime();
#pragma omp sections
	{
#pragma omp section
		{
			for (int j = 0; j < s1; j++)
			{
				//� ����� ���������� ������� �� �������� ������� C
				for (int i = 0; i < matrixA.n; i++)
				{
					matrixC.elements[i][j] = 0;
					//� ����� ���������� ������������ ��������� ������ A � B �� ��������������� ������
					for (int t = 0; t < matrixA.m; t++)
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
				//� ����� ���������� ������� �� �������� ������� C
				for (int i = 0; i < matrixA.n; i++)
				{
					matrixC.elements[i][j] = 0;
					//� ����� ���������� ������������ ��������� ������ A � B �� ��������������� ������
					for (int t = 0; t < matrixA.m; t++)
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
					//� ����� ���������� ������� �� �������� ������� C
					for (int i = 0; i < matrixA.n; i++)
					{
						matrixC.elements[i][j] = 0;
						//� ����� ���������� ������������ ��������� ������ A � B �� ��������������� ������
						for (int t = 0; t < matrixA.m; t++)
						{
							matrixC.elements[i][j] += matrixA.elements[t][i] * matrixB.elements[t][j];
						}
					}
				}
		}
#pragma omp section
		{
			if (numThreads > 3)
				for (int j = s3; j < matrixB.n; j++)
				{
					//� ����� ���������� ������� �� �������� ������� C
					for (int i = 0; i < matrixA.n; i++)
					{
						matrixC.elements[i][j] = 0;
						//� ����� ���������� ������������ ��������� ������ A � B �� ��������������� ������
						for (int t = 0; t < matrixA.m; t++)
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

void MultiplicationStrassenConsistently()
{

}

void CalculateMultFuncTrustedTime(void* func, int m, int n, int k, double& time, int iterations)
{
	Matrix* matrixA = new Matrix(n, m);
	Matrix* matrixB = new Matrix(k, n);

	//���������� ������
	FillLineConsistently(*matrixA, time);
	FillLineConsistently(*matrixB, time);

	//������� A (n x m) ���������� �� ������� B (k x n)
	//������� ��������������� ������� B
	matrixB->Transpose();
	//������ ������� B ��������������� (n x k)

	Matrix* matrixC = new Matrix(m, k);

	double curtime = 0, avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	double* Times = new double[iterations];

	for (int i = 0; i < iterations; i++)
	{
		(*(MultFunctTempl)func)(*matrixA, *matrixB, *matrixC, curtime);
		Times[i] = curtime;
		avgTime += curtime;
	}
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	time = avgTimeT;
}

void CalculateAllFillFuncs(int m, int n)
{
	double time;
	string* funcsNames = new string[5]{ "��������� ���������������� ����������", "��������� ������������ ���������� � for static", "��������� ������������ ���������� � for dynamic", "��������� ������������ � sections", "������� ������������ ���������� � for static" };
	void** funcs = new void*[5]{ FillLineConsistently, FillLineParallelForStatic, FillLineParallelForDynamic, FillLineParallelSections, FillBlockParallelForStatic };

	cout << "��:" << n << "x" << m << endl;
	for (int t = 2; t < 5; t++)
	{
		omp_set_num_threads(t);
		cout << "����� �������:" << t << endl;
		for (int i = 0; i < 5; i++)
		{
			time = 0;
			CalculateFillFuncTrustedTime(funcs[i], m, n, time, 30);
			cout << funcsNames[i] << ". ������� ���������: " << time * 1000 << endl;
		}
		cout << endl;
	}
}

void CalculateAllMultFuncs(int m, int n, int k)
{
	double time;
	string* funcsNames = new string[4]{ "���������������� ���������", "������������ ��������� � for static", "������������ ��������� � for dynamic", "������������ ��������� � sections" };
	void** funcs = new void*[4]{ MultiplicationConsistently, MultiplicationParallelForStatic, MultiplicationParallelForDynamic, MultiplicationParallelSections };

	cout << "��:" << n << "x" << m << endl;
	for (int t = 2; t < 5; t++)
	{
		omp_set_num_threads(t);
		cout << "����� �������:" << t << endl;
		for (int i = 0; i < 4; i++)
		{
			time = 0;
			CalculateMultFuncTrustedTime(funcs[i], m, n, k, time, 30);
			cout << funcsNames[i] << ". ������� ���������: " << time * 1000 << endl;
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

	string* funcsNames = new string[1]{ "���������������� ������������" };
	void** funcs = new void*[1]{ MultiplicationConsistently };

	for (int d = 0; d < 4; d++)
	{
		for (int t = 2; t < 5; t++)
		{
			omp_set_num_threads(t);
			cout << "����� �������:" << t << endl;
			for (int i = 0; i < 5; i++)
			{
				time = 0;
				CalculateMultFuncTrustedTime(funcs[i], n, m, k, time, 30);
				cout << funcsNames[i] << ". ������� ���������: " << time * 1000 << endl;
			}
			cout << endl;
		}
	}
}

void Task3()
{
	int* n = new int[4]{ 500, 750, 900, 1200 };
	int* m = new int[4]{ 700, 750, 1000, 1300 };
	int* k = new int[4]{ 600, 900, 1000, 1500 };
	
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
	/*
	double time;
	Matrix* matrixA = new Matrix(2, 3);
	Matrix* matrixB = new Matrix(4, 2);
	Matrix* matrixC = new Matrix(4, 3);
	FillTest(*matrixA);
	FillTest(*matrixB);
	PrintMatrix(*matrixA);
	PrintMatrix(*matrixB);
	matrixB->Transpose();
	MultiplicationConsistently(*matrixA, *matrixB, *matrixC, time);
	PrintMatrix(*matrixC);
	*/

	cout << "[1]: ������� 1\n"
		<< "[2]: ������� 2\n"
		<< "[3]: ������� 3\n";
	cin >> choice;

	if (choice == 1)
	{
		Task1();
		/*
		int numThreads;
		cout << "������� ���������� �������. ������������ ��������� ���������� ������� - " << omp_get_max_threads() << endl;
		cin >> numThreads;

		cout << "����� � ������� [1]: printf [2]: cout\n";
		cin >> choice;

			//Task1PrintfHelloWorld(numThreads);
		//else
			//Task1CoutHelloWorld(numThreads);
		if (choice == 1)
		*/
	}
	else if (choice == 2)
	{
		//Task2();
	}
	else
	{
		Task3();
	}

	cout << endl << "������ ��������� ���������" << endl;
	string wait;
	//cin >> wait;
}