#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <String>
#include <ctime>
#include <Math.h>
#include <fstream>
#include <list>

using namespace std;

string GetTime(double time)
{
	char buffer[_MAX_U64TOSTR_BASE2_COUNT];
	string result = _itoa(trunc(time), buffer, 10);
	result += '.';
	for (int i = 0; i < 4; i++)
	{
		time = time * 10;
		result += _itoa((int)trunc(time) % 10, buffer, 10);
	}

	return result;
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
	int i;
	step = 1.0 / (double)num_steps;
	double startTime, endTime;
	startTime = omp_get_wtime();
#pragma omp schedule(static) for
	{
		for (i = 0; i < num_steps; i++)
		{
			x = (i + 0.5) * step;
			sum = sum + 4.0 / (1.0 + x * x);
		}
	}
	pi = step * sum;
	endTime = omp_get_wtime();

	time = endTime - startTime;
	return pi;
}

double CalculatePiParallelForDynamic(long num_steps, double& time)
{
	double step, pi, x, sum = 0.0;
	int i;
	step = 1.0 / (double)num_steps;
	double startTime, endTime;
	startTime = omp_get_wtime();
#pragma omp schedule(dynamic) for
	{
		for (i = 0; i < num_steps; i++)
		{
			x = (i + 0.5) * step;
			sum = sum + 4.0 / (1.0 + x * x);
		}
	}
	pi = step * sum;
	endTime = omp_get_wtime();

	time = endTime - startTime;
	return pi;
}

double CalculatePiParallelForGuided(long num_steps, double& time)
{
	double step, pi, x, sum = 0.0;
	int i;
	step = 1.0 / (double)num_steps;
	double startTime, endTime;
	startTime = omp_get_wtime();
#pragma omp schedule(guided) for
	{
		for (i = 0; i < num_steps; i++)
		{
			x = (i + 0.5) * step;
			sum = sum + 4.0 / (1.0 + x * x);
		}
	}
	pi = step * sum;
	endTime = omp_get_wtime();

	time = endTime - startTime;
	return pi;
}

double CalculatePiParallelSections(long num_steps, double& time)
{
	double step, pi, x, sum = 0.0;
	int i;
	step = 1.0 / (double)num_steps;
	int threadsNum = omp_get_num_threads();
	double startTime, endTime;
	int s1 = 0, s2 = 0, s3 = 0;
	if (threadsNum == 1)
	{
		s1 = num_steps;
	}
	else if (threadsNum >= 2)
	{
		s1 = num_steps / threadsNum;
		s2 = num_steps;
		s3 = num_steps;

		if (threadsNum >= 3)
		{
			s2 = num_steps / threadsNum * 2;

			if (threadsNum >= 4)
			{
				s3 = num_steps / threadsNum * 3;
			}
		}
	}

	startTime = omp_get_wtime();
#pragma omp sections
	{
#pragma omp section
		{
			for (i = 0; i < s1; i++)
			{
				x = (i + 0.5) * step;
				sum = sum + 4.0 / (1.0 + x * x);
			}
		}
#pragma omp section
		{
			for (i = s1; i < s2; i++)
			{
				x = (i + 0.5) * step;
				sum = sum + 4.0 / (1.0 + x * x);
			}
		}
#pragma omp section
		{
			for (i = s2; i < s3; i++)
			{
				x = (i + 0.5) * step;
				sum = sum + 4.0 / (1.0 + x * x);
			}
		}
#pragma omp section
		{
			for (i = s3; i < num_steps; i++)
			{
				x = (i + 0.5) * step;
				sum = sum + 4.0 / (1.0 + x * x);
			}
		}
	}
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
	string* funcsNames = new string[6]{"Последовательная реализация", "Параллельная реализация FOR (static)", "Параллельная реализация FOR (dynamic)", "Параллельная реализация FOR (guided)", "Параллельная реализация Sections", "Свой алгоритм"};
	void** funcs = new void* [6]{ CalculatePiConsistenty, CalculatePiParallelForStatic, CalculatePiParallelForDynamic, CalculatePiParallelForGuided, CalculatePiParallelSections, CalculatePiMyAlgorithm};
	
	static long* nums_steps = new long[4]{1000000, 5000000, 10000000, 15000000};
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
				cout << funcsNames[i] << " .Число Пи = " << pi << " .Времени затрачено: " << GetTime(time) << endl;
			}
			cout << endl;
		}
	}
}

int main()
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
	cin >> wait;
}