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

	//Функция транспонирования матрицы
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
				cout << funcsNames[i] << " .Число Пи = " << pi << " .Времени затрачено: " << time * 1000 << endl;
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
			matrix.elements[i][j] = pow(i, 3 / 4)* cos(i) / atan(i) + pow(i, 3 / 4)* cos(j) / atan(j);
		}
	}
	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void FillLineParallelForStatic(Matrix& matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
#pragma omp parallel for schedule(static, num_steps / 12)
	for (int j = 0; j < matrix.n; j++)
	{
		for (int i = 0; i < matrix.m; i++)
		{
			matrix.elements[i][j] = pow(i, 3 / 4)* cos(i) / atan(i) + pow(i, 3 / 4)* cos(j) / atan(j);
		}
	}
	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void FillLineParallelForDynamic(Matrix& matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, num_steps / 6)
	for (int j = 0; j < matrix.n; j++)
	{
		for (int i = 0; i < matrix.m; i++)
		{
			matrix.elements[i][j] = pow(i, 3 / 4)* cos(i) / atan(i) + pow(i, 3 / 4)* cos(j) / atan(j);
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
					matrix.elements[i][j] = pow(i, 3 / 4)* cos(i) / atan(i) + pow(j, 3 / 4)* cos(j) / atan(j);
				}
			}
		}
#pragma omp section
		{
			for (int j = s1; j < s2; j++)
			{
				for (int i = 0; i < matrix.m; i++)
				{
					matrix.elements[i][j] = pow(i, 3 / 4)* cos(i) / atan(i) + pow(j, 3 / 4)* cos(j) / atan(j);
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
						matrix.elements[i][j] = pow(i, 3 / 4)* cos(i) / atan(i) + pow(j, 3 / 4)* cos(j) / atan(j);
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
						matrix.elements[i][j] = pow(i, 3 / 4)* cos(i) / atan(i) + pow(j, 3 / 4)* cos(j) / atan(j);
					}
				}
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void GetBlockStartPosition(int i, int n, int sizeX, int& shiftX, int& shiftY)
{
	shiftX = (i * sizeX) % n;
	shiftY = trunc((i * sizeX) / n);
}

void FillBlockParallelForStatic(Matrix& matrix, double& time)
{
	double startTime, endTime;
	startTime = omp_get_wtime();
	int blockSizeX = 5;
	int blockSizeY = 5;
	int blocksAmount = matrix.n * matrix.m / (blockSizeX * blockSizeY);
	int shiftX = 0, shiftY = 0;
	int n = matrix.n;

#pragma omp parallel for schedule(static, num_steps / 12)
	for (int i = 0; i < blocksAmount; i++)
	{
		GetBlockStartPosition(i, n, blockSizeX, shiftX, shiftY);
		for (int t = 0; t < blockSizeY; t++)
		{
			for (int d = 0; d < blockSizeX; d++)
			{
				int dIndex = d + shiftX;
				int tIndex = t + shiftY;
				matrix.elements[dIndex][tIndex] = pow(dIndex, 3 / 4)* cos(dIndex) / atan(dIndex) + pow(tIndex, 3 / 4)* cos(tIndex) / atan(tIndex);;
			}
		}
	}

	endTime = omp_get_wtime();
	time = endTime - startTime;
}

void PrintMatrix(Matrix matrix)
{
	for (int j = 0; j < matrix.n; j++)
	{
		for (int i = 0; i < matrix.m; i++)
		{
			cout << setw(5) << (int)matrix.elements[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

double CalculateFillFuncTrustedTime(void* func, int n, int m, double& time, int iterations)
{
	Matrix* matrixA = new Matrix(n, m);
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

void CalculateMultFuncTrustedTime(void* func, int n, int m, int k, double& time, int iterations)
{
	Matrix* matrixA = new Matrix(n, m);
	Matrix* matrixB = new Matrix(m, k);

	//Заполнение матриц
	FillLineConsistently(*matrixA, time);
	FillLineConsistently(*matrixB, time);

	Matrix* matrixC = new Matrix(n, k);

	static long* nums_steps = new long[4]{ 1000000, 5000000, 10000000, 15000000 };
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
				cout << funcsNames[i] << " .Число Пи = " << pi << " .Времени затрачено: " << time * 1000 << endl;
			}
			cout << endl;
		}
	}
}

void Task2()
{
	int* n = new int[4]{ 1000, 2000, 3000, 4000 };
	int* m = new int[4]{ 1500, 2000, 4500, 6000 };
	int* k = new int[4]{ 2000, 3500, 4500, 3000 };
	
	double time;
	//Способы заполнения матриц TODO:
	string* funcNames = new string[5]{ "Ленточное последовательное заполнение", "Ленточное параллельное заполнение с for static", "Ленточное параллельное заполнение с for dynamic", "Ленточное параллельный с sections", "Блочное параллельное заполнение с for static" };
	void** funcs = new void*[5]{ FillLineConsistently, FillLineParallelForStatic, FillLineParallelForDynamic, FillLineParallelSections, FillBlockParallelForStatic };
}

void main()
{
	srand(time(0));
	setlocale(LC_ALL, "Russian");
	int choice;

	cout << "[1]: Задание 1\n"
		<< "[2]: Задание 2\n"
		<< "[3]: Задание 3\n";
	cin >> choice;

	if (choice == 1)
	{
		Task1();
		/*
		int numThreads;
		cout << "Введите количество потоков. Максимальное доступное количество потоков - " << omp_get_max_threads() << endl;
		cin >> numThreads;

		cout << "Вывод с помощью [1]: printf [2]: cout\n";
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
		//Task3();
	}

	cout << endl << "Работа программы завершена" << endl;
	string wait;
	//cin >> wait;
}