#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <String>
#include <ctime>
#include <Math.h>
#include <fstream>
#include <list>

using namespace std;

double CalculatePiConsistenty(long num_steps, double& time)
{
	double step, pi, x, sum = 0.0;
	int i;
	step = 1.0 / (double)num_steps;
	for (i = 0; i < num_steps; i++)
	{
		x = (i + 0.5) * step;
		sum = sum + 4.0 / (1.0 + x * x);
	}
	pi = step * sum;

	return pi;
}

double CalculatePiParallelForStatic(long num_steps, double& time)
{
	double step, pi, x, sum = 0.0;
	int i;
	step = 1.0 / (double)num_steps;
#pragma omp for static
	{
		for (i = 0; i < num_steps; i++)
		{
			x = (i + 0.5) * step;
			sum = sum + 4.0 / (1.0 + x * x);
		}
	}
	pi = step * sum;

	return pi;
}

double CalculatePiParallelSections(long num_steps, double& time)
{
	double step, pi, x, sum = 0.0;
	int i;
	step = 1.0 / (double)num_steps;
	int = 
#pragma omp sections
	{
#pragma omp section
		{

		}
		for (i = 0; i < num_steps; i++)
		{
			x = (i + 0.5) * step;
			sum = sum + 4.0 / (1.0 + x * x);
		}
	}
	pi = step * sum;

	return pi;
}

void Task1()
{
	void** fincs = new void* [1]{ CalculatePiConsistenty };
	static long num_steps = 500000;
	printf("Pi = %f\n", pi);
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
		int numThreads;
		cout << "Введите количество потоков. Максимальное доступное количество потоков - " << omp_get_max_threads() << endl;
		cin >> numThreads;

		cout << "Вывод с помощью [1]: printf [2]: cout\n";
		cin >> choice;

		if (choice == 1)
			Task1PrintfHelloWorld(numThreads);
		else
			Task1CoutHelloWorld(numThreads);
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
	cin >> wait;
}