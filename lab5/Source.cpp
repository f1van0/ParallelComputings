#include<tbb/tbb.h>
#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <String>
#include <ctime>
#include <Windows.h>
#include <Math.h>
#include <fstream>
#include <sstream>
using namespace tbb;
using namespace std;


typedef void(*IntMatrix)(int**&, int**&, int**&, int);
typedef void(*DoubleMatrix)(double**&, double**&, double**&, int);

#pragma region Classes

/* ---- Алгоритм вычисления суммы  ---- */

// TBB Reduction

/// Используемый функцией класс
template<class T>
class reduce_total {
	T* m1;
	T* m2;
	int size;

public:
	double sumRes1, sumRes2;

	void operator() (blocked_range<size_t>& r) {
		T sumRes1_local = sumRes1, sumRes2_local = sumRes2;
		T* m1_local = m1;
		T* m2_local = m2;

		size_t end = r.end();
		for (size_t i = r.begin(); i != end; ++i) {
			sumRes1_local += m1_local[i];
			sumRes2_local += m2_local[i];
		}

		sumRes1 = sumRes1_local;
		sumRes2 = sumRes2_local;
	}

	reduce_total(reduce_total& r, split) : m1(r.m1), m2(r.m2), size(r.size), sumRes1(0.0), sumRes2(0.0) {} // Split
	void join(const reduce_total& r) { sumRes1 += r.sumRes1; sumRes2 += r.sumRes2; } // Join

	reduce_total(T* orig_m1, T* orig_m2, int orig_size) : m1(orig_m1), m2(orig_m2), size(orig_size), sumRes1(0.0), sumRes2(0.0) {} // Основной конструктор
};

/* ---- Алгоритм нахождения максимума  ---- */

// TBB Reduction

// Используемый функцией класс
template<class T>
class reduce_max {
private:
	T* ArrayA;
	T* ArrayB;
public:
	T MaxValue1, MaxValue2;

	void operator() (const blocked_range<size_t>& r) {
		T* a = ArrayA;
		T* b = ArrayB;
		T val1, val2;

		for (size_t i = r.begin(); i != r.end(); ++i) {
			val1 = a[i];
			val2 = b[i];

			if (MaxValue1 < val1) MaxValue1 = val1;
			if (MaxValue2 < val2) MaxValue2 = val2;
		}
	}

	reduce_max(reduce_max& x, tbb::split) : ArrayA(x.ArrayA), ArrayB(x.ArrayB), MaxValue1(-1000), MaxValue2(-1000) {}
	void join(const reduce_max& y) {
		if (y.MaxValue1 > MaxValue1) MaxValue1 = y.MaxValue1;
		if (y.MaxValue2 > MaxValue2) MaxValue2 = y.MaxValue2;
	}

	reduce_max(T* A, T* B) : ArrayA(A), ArrayB(B), MaxValue1(-1000), MaxValue2(-1000) {}
};
#pragma endregion


#pragma region SumMatrix

//Функция суммирует матрицы
template<class T>
void SumMatrix(T**& matrix1, T**& matrix2, T**& result, int size)
{
	//double time_start = omp_get_wtime();

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
		{
			result[i][j] = matrix1[i][j] + matrix2[i][j];
		}

	//double time_stop = omp_get_wtime();
	//return time_stop - time_start;
}

//Функция суммирует матрицы с OMP
template<class T>
void SumMatrixOMP(T**& matrix1, T**& matrix2, T**& result, int size)
{
	//double time_start = omp_get_wtime();

#pragma omp parallel for
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
		{
			result[i][j] = matrix1[i][j] + matrix2[i][j];
		}

	//double time_stop = omp_get_wtime();
	//return time_stop - time_start;
}

//Функция суммирует матрицы с TBB
template<class T>
void SumMatrixTBB(T**& matrix1, T**& matrix2, T**& result, int size)
{
	//double time_start = omp_get_wtime();

	tbb::parallel_for(tbb::blocked_range2d<int>(0, size, 0, size), [&](tbb::blocked_range2d<int> r)
	{
		for (int i = r.rows().begin(); i < r.rows().end(); i++)
			for (int j = r.cols().begin(); j < r.cols().end(); j++)
			{
				//matrix[i][j] = (cos(i) * sin(j / 2.0));
				result[i][j] = matrix1[i][j] + matrix2[i][j];
			}
	});

	//double time_stop = omp_get_wtime();
	//return time_stop - time_start;
}

#pragma endregion

#pragma region ProductMatrix

//Функция вычисляет произведение матриц
template<class T>
void ProductMatrix(T**& matrix1, T**& matrix2, T**& result, int size)
{
	//double time_start = omp_get_wtime();

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
		{
			//matrix[i][j] = (cos(i) * sin(j / 2.0));
			result[i][j] = matrix1[i][j] * matrix2[i][j];
		}

	//double time_stop = omp_get_wtime();
	//return time_stop - time_start;
}

//Функция вычисляет произведение матриц с OMP
template<class T>
void ProductMatrixOMP(T**& matrix1, T**& matrix2, T**& result, int size)
{
	//double time_start = omp_get_wtime();

#pragma omp parallel for
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
		{
			//matrix[i][j] = (cos(i) * sin(j / 2.0));
			result[i][j] = matrix1[i][j] * matrix2[i][j];
		}

	//double time_stop = omp_get_wtime();
	//return time_stop - time_start;
}

//Функция вычисляет произведение матриц с TBB
template<class T>
void ProductMatrixTBB(T**& matrix1, T**& matrix2, T**& result, int size)
{
	//double time_start = omp_get_wtime();

	tbb::parallel_for(tbb::blocked_range2d<int>(0, size, 0, size), [&](tbb::blocked_range2d<int> r)
	{
		for (int i = r.rows().begin(); i < r.rows().end(); i++)
			for (int j = r.cols().begin(); j < r.cols().end(); j++)
			{
				//matrix[i][j] = (cos(i) * sin(j / 2.0));
				result[i][j] = matrix1[i][j] * matrix2[i][j];
			}
	});

	//double time_stop = omp_get_wtime();
	//return time_stop - time_start;
}

#pragma endregion

#pragma region TotalSumMatrix

//Функция суммирует элементы матриц
template<class T>
void TotalSum(T**& matrix1, T**& matrix2, T**& unused, int size) {

	//double time_start = omp_get_wtime();

	T total1 = 0, total2 = 0;

	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			total1 += matrix1[i][j];
			total2 += matrix2[i][j];
		}
	}
	cout << total1 << " " << total2;
	//double time_stop = omp_get_wtime();
	//return time_stop - time_start;
}

//Функция суммирует элементы матриц с OMP reduction
template<class T>
void TotalSumOMP(T**& matrix1, T**& matrix2, T**& unused, int size) {

	//double time_start = omp_get_wtime();

	T total1 = 0, total2 = 0;

	for (int j = 0; j < size; j++)
	{
#pragma omp parallel
#pragma omp for reduction(+:total1)
		for (int i = 0; i < size; i++)
		{
			total1 += matrix1[i][j];
		}

#pragma omp parallel
#pragma omp for reduction(+:total2)
		for (int i = 0; i < size; i++)
		{
			total2 += matrix2[i][j];
		}
	}
	cout << total1 << " " << total2;

	//double time_stop = omp_get_wtime();
	//return time_stop - time_start;
}

//Функция суммирует элементы матриц с TBB
template<class T>
void TotalSumTBB(T**& matrix1, T**& matrix2, T**& unused, int size) {

	//double time_start = omp_get_wtime();

	T total1 = 0, total2 = 0;

	for (int i = 0; i < size; i++) {
		reduce_total<T> r(matrix1[i], matrix2[i], size);
		parallel_reduce(blocked_range<size_t>(0, size), r);
		total1 += r.sumRes1;
		total2 += r.sumRes2;
	}
	cout << total1 << " " << total2;

	//double time_stop = omp_get_wtime();
	//return time_stop - time_start;
}

#pragma endregion

#pragma region MaxElemMatrix

//Функция находит максимальный элемент в 2 матрицах с OMP for
template<class T>
void MaxElem(T**& matrix1, T**& matrix2, T**& unused, int size) {

	//double time_start = omp_get_wtime();

	double max1 = -1000, max2 = -1000;

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (max1 < matrix1[i][j]) max1 = matrix1[i][j];
			if (max2 < matrix1[i][j]) max2 = matrix1[i][j];
		}
	}
	cout << max1 << " " << max2;

	//double time_stop = omp_get_wtime();
	//return time_stop - time_start;
}

//Функция находит максимальный элемент в 2 матрицах с OMP for
template<class T>
void MaxElemOMP(T**& matrix1, T**& matrix2, T**& unused, int size) {

	//double time_start = omp_get_wtime();

	double max1 = -1000, max2 = -1000;

	for (int i = 0; i < size; i++) {
		int j = 0;

#pragma omp parallel for private(j)
		for (j = 0; j < size; j++) {
			if (max1 < matrix1[i][j]) max1 = matrix1[i][j];
			if (max2 < matrix1[i][j]) max2 = matrix1[i][j];
		}
	}
	cout << max1 << " " << max2;
	//double time_stop = omp_get_wtime();
	//return time_stop - time_start;
}

//Функция находит максимальный элемент в 2 матрицах с TBB
template<class T>
void MaxElemTBB(T**& matrix1, T**& matrix2, T**& unused, int size) {

	//double time_start = omp_get_wtime();

	double max1 = -1000, max2 = -1000;

	for (int i = 0; i < size; i++) {
		reduce_max<T> r(matrix1[i], matrix2[i]);
		parallel_reduce(blocked_range<size_t>(0, size), r);
		if (max1 < r.MaxValue1) max1 = r.MaxValue1;
		if (max2 < r.MaxValue2) max2 = r.MaxValue2;
	}
	cout << max1 << " " << max2;
	//double time_stop = omp_get_wtime();
	//return time_stop - time_start;
}

#pragma endregion

#pragma region specFunc
//Заполнение матрицы
template<class T>
void Fill2DMatrix(T** matrix, int size)
{
	tbb::parallel_for(tbb::blocked_range2d<int>(0, size, 0, size), [&](tbb::blocked_range2d<int> r)
	{
		for (int i = r.rows().begin(); i < r.rows().end(); i++)
			for (int j = r.cols().begin(); j < r.cols().end(); j++)
			{
				//matrix[i][j] = (cos(i) * sin(j / 2.0));
				matrix[i][j] = 1000 * pow(i + 1, 3 / 4) * cos(i) / atan(i + 1) + 1000 * pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
			}
	});
}

//Вывести матрицы
template<class T>
void Print2DMatrix(T** matrix, int size)
{
	for (int i = 0; i < size; i++)
			{
				for (int j = 0; j < size; j++)
				{
					//matrix[i][j] = (cos(i) * sin(j / 2.0));
					cout << matrix[i][j] << " ";
				}
				cout << endl;
			}
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

#pragma endregion

void CalcMatrixAvgTimeInt(IntMatrix func, int size, double& time, int iterations)
{
	int** matrix1 = new int*[size];
	int** matrix2 = new int*[size];
	int** result = new int*[size];

	for (int i = 0; i < size; i++)
	{
		matrix1[i] = new int[size];
		matrix2[i] = new int[size];
		result[i] = new int[size];
	}

	Fill2DMatrix<int>(matrix1, size);
	Fill2DMatrix<int>(matrix2, size);

	double avgTime = 0, avgTimeT = 0, correctAVG = 0;
	double startTime, curTime;
	double* Times = new double[iterations];

	std::cout << "[";
	for (int i = 0; i < iterations; i++)
	{
		startTime = omp_get_wtime();
		func(matrix1, matrix2, result, size);
		curTime = omp_get_wtime() - startTime;
		Times[i] = curTime;
		avgTime += curTime;
		std::cout << "#";
	}
	std::cout << "]\n";
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	time = avgTimeT * 1000;

	for (int i = 0; i < size; i++)
	{
		delete[] matrix1[i];
		delete[] matrix2[i];
		delete[] result[i];
	}
	delete matrix1;
	delete matrix2;
	delete result;
}

void CalcMatrixAvgTimeDouble(DoubleMatrix func, int size, double& time, int iterations)
{
	double** matrix1 = new double*[size];
	double** matrix2 = new double*[size];
	double** result = new double*[size];

	for (int i = 0; i < size; i++)
	{
		matrix1[i] = new double[size];
		matrix2[i] = new double[size];
		result[i] = new double[size];
	}

	Fill2DMatrix<double>(matrix1, size);
	Fill2DMatrix<double>(matrix2, size);

	double avgTime = 0, avgTimeT = 0, correctAVG = 0;
	double startTime, curTime;
	double* Times = new double[iterations];

	std::cout << "[";
	for (int i = 0; i < iterations; i++)
	{
		startTime = omp_get_wtime();
		func(matrix1, matrix2, result, size);
		curTime = omp_get_wtime() - startTime;
		Times[i] = curTime;
		avgTime += curTime;
		std::cout << "#";
	}
	std::cout << "]\n";
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	time = avgTimeT * 1000;

	for (int i = 0; i < size; i++)
	{
		delete[] matrix1[i];
		delete[] matrix2[i];
		delete[] result[i];
	}
	delete matrix1;
	delete matrix2;
	delete result;
}

void CalcMatrixFuncs()
{
	std::ofstream resultsFile;
	int* dataAmount = new int[4]{ 3000, 4500, 6000, 8000 };

	double time;
	string* funcsNames = new string[12]{ "Суммирование матриц (последовательно)", "Суммирование матриц (Open MP)", "Суммирование матриц (TBB)", "Перемножение матриц (последовательно)", "Перемножение матриц (Open MP)", "Перемножение матриц (TBB)", "Вычисление суммы элементов 2 матриц (последовательно)", "Вычисление суммы элементов 2 матриц (Open MP)", "Вычисление суммы элементов 2 матриц (TBB)", "Поиск максимального элемента в 2 матрицах (последовательно)", "Поиск максимального элемента в 2 матрицах (Open MP)", "Поиск максимального элемента в 2 матрицах (TBB)" };
	//IntMatrix* funcInt = new IntMatrix[3]{ SumMatrix<int>, SumMatrixOMP<int>, SumMatrixTBB<int> };
	DoubleMatrix* funcDouble = new DoubleMatrix[12]{ SumMatrix<double>, SumMatrixOMP<double>, SumMatrixTBB<double>, ProductMatrix<double>, ProductMatrixOMP<double>, ProductMatrixTBB<double>, TotalSum<double>, TotalSumOMP<double>, TotalSumTBB<double>, MaxElem<double>, MaxElemOMP<double>, MaxElemTBB<double> };
	IntMatrix* funcInt = new IntMatrix[12]{ SumMatrix<int>, SumMatrixOMP<int>, SumMatrixTBB<int>, ProductMatrix<int>, ProductMatrixOMP<int>, ProductMatrixTBB<int>, TotalSum<int>, TotalSumOMP<int>, TotalSumTBB<int>, MaxElem<int>, MaxElemOMP<int>, MaxElemTBB<int> };
	double* T1 = new double[3];
	resultsFile.open("Task1Results.csv", std::ios_base::app);
	resultsFile << "Int;;;;;;;;;;Double;\n";
	resultsFile << "Функция сортировки;Потоки;Время;Sp(n);Время;Sp(n);Время;Sp(n);Время;Sp(n);Время;Sp(n);Время;Sp(n);Время;Sp(n);Время;Sp(n);\n";

	for (int j = 0; j < 2; j++)
	{
		if (j == 0)
		{
			cout << "Тип данных Int";
		}
		else
		{
			cout << "Тип данных Double";
		}

		resultsFile << ";;";
		for (int cr = 0; cr < 4; cr++)
			resultsFile << "НД" << cr + 1 << ": " << dataAmount[cr] << ";;";

		for (int cr = 0; cr < 4; cr++)
			resultsFile << "НД" << cr + 1 << ": " << dataAmount[cr] << ";;";

		resultsFile << endl;

		for (int i = 0; i < 12; i++)
		{
			resultsFile << funcsNames[i];
			std::cout << funcsNames[i] << endl;
			for (int t = 2; t < 5; t++)
			{
				if (i == 0 || i == 3 || i == 6 || i == 9)
					t = 1;

				std::cout << "Потоков: " << t << endl;
				resultsFile << ";" << t << ";";
				omp_set_num_threads(t);

				for (int d = 0; d < 4; d++) //4
				{
					std::cout << "Количество элементов: " << dataAmount[d] << endl;

					if (j == 0)
						CalcMatrixAvgTimeInt(funcInt[i], dataAmount[d], time, 20);
					else
						CalcMatrixAvgTimeDouble(funcDouble[i], dataAmount[d], time, 20);

					if (i == 0 || i == 3 || i == 6 || i == 9)
						T1[d] = time;
					resultsFile << time << ";" << T1[d] / time << "^;";
					std::cout << " - Длительность: " << time << endl;
				}
				std::cout << endl;
				resultsFile << endl;

				if (i == 0 || i == 3 || i == 6 || i == 9)
					break;
			}
		}
		resultsFile << endl;
	}
}

int main()
{
	setlocale(LC_ALL, "Russian");
	int wait;
	//task_group tg;
	//tg.run(MyTask("1")); /* Генерация задачи */
	//tg.run(MyTask("2")); /* Генерация задачи */
	//tg.wait(); /* Ожидание завершения задачи */

	CalcMatrixFuncs();
	
	std::cin >> wait;
	return 0;
}