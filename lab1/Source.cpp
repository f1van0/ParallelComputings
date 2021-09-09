#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <String>
#include <ctime>
#include <Math.h>

using namespace std;

void TestPrintfHelloWorld()
{
	printf("Hello, World!\n");
#pragma omp parallel num_threads (4)
	{
		int i, n;
		i = omp_get_thread_num();
		n = omp_get_num_threads();
		printf("I’m thread %d of %d\n", i, n);
	}

	int wait;
	scanf("%d", &wait);
}

void Task1PrintfHelloWorld(int numThreads)
{
	//вывод с помощью printf и parallel for
#pragma omp parallel num_threads(numThreads)
	{
		printf("\nHello World! from thread %d of %d", omp_get_thread_num(), omp_get_max_threads());
	}
}

void Task1CoutHelloWorld(int numThreads)
{
	//вывод с помощью printf и parallel for
#pragma omp parallel num_threads(numThreads)
	{
		cout << "\nHello World! from thread " << omp_get_thread_num() << " of " << omp_get_max_threads() << endl;
	}
}

struct Vector
{
	double* numbers;
	int length;

	Vector()
	{
		numbers = nullptr;
		length = 0;
	}

	Vector(int _length)
	{
		numbers = new double[_length];
		length = _length;
	}
};

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
		int nVectors = 2;
		int nNumbers = 100000000;

		//Последовательное
		//1 этап - Заполнение векторов
		Vector* vectors = new Vector[nVectors];
		//Инициализация
		for (int i = 0; i < nVectors; i++)
		{
			vectors[i] = Vector(nNumbers);
		}

		double time_Start = omp_get_wtime();
		for (int i = 0; i < nVectors; i++)
		{
			for (int j = 0; j < nNumbers; j++)
			{
				int x = 1 + rand() % 10;
				vectors[i].numbers[j] = pow(x, 2/3) * (cos(x) / atan(x));
			}
		}
		cout << "Последовательно. 1 Этап. Среднее время заполнения 1 вектора составило - " << (omp_get_wtime() - time_Start) / nVectors << ". Всего заполнялось " << nVectors << " векторов\n";

		//Последовательное
		//2 этап - Сложение векторов
		Vector resultVector(nNumbers);

		time_Start = omp_get_wtime();
		for (int i = 0; i < nVectors; i++)
		{
			for (int j = 0; j < nNumbers; j++)
			{
				resultVector.numbers[j] += vectors[i].numbers[j];
			}
		}
		cout << "Последовательно. 2 этап. Время сложения векторов составило - " << omp_get_wtime() - time_Start << endl;

		//Последовательное
		//3 этап - Подсчет суммы всех элементов итогового вектора
		double sum = 0;
		time_Start = omp_get_wtime();
		for (int j = 0; j < nNumbers; j++)
		{
			sum += resultVector.numbers[j];
		}
		cout << "Последовательно. 3 этап. Время подсчета суммы всех элементов итогового вектора составило - " << omp_get_wtime() - time_Start << endl;

		//Параллельно
		//1 этап - Заполнение векторов
		vectors = new Vector[nVectors];
		//Инициализация
		for (int i = 0; i < nVectors; i++)
		{
			vectors[i] = Vector(nNumbers);
		}

		time_Start = omp_get_wtime();
#pragma omp parallel
		{
			for (int i = 0; i < nVectors; i++)
			{
				for (int j = 0; j < nNumbers; j++)
				{
					int x = 1 + rand() % 10;
					vectors[i].numbers[j] = pow(x, 2 / 3) * (cos(x) / atan(x));
				}
			}
		}
		cout << "Параллельно. 1 Этап. Среднее время заполнения 1 вектора составило - " << (omp_get_wtime() - time_Start) / nVectors << ". Всего заполнялось " << nVectors << " векторов\n";
	}
	else
	{

	}

	cout << endl << endl;
	string wait;
	cin >> wait;
}
