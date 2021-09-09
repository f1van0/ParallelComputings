// Experiments_Demo.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <string>
#include <omp.h>

using namespace std;

// Определение типа для функций тестирования (имеют общуу структуру но вызывают соответствующие функции)
typedef double(*TestFunctTempl)(double*&, double*&, double*&, int&);
// Функция рассчета среднеарифметического значения в доверительном интервале
double AvgTrustedInterval(double& avg, double*& times, int& cnt);
// Функция для замера времени Функции работы с данными 
double TestIter(void* Funct, double*& arr_in1, double*& arr_in2, double*& arr_res, int size, int iterations);

/// Функции работы с данными
// заполнение массивов два варианта
double Fillarr1(double* arr, int size);
double Fillarr2(double* arr, int size);
// Сложение массивов 
double SumArrays(double*& arr_in1, double*& arr_in2, double*& arr_res, int& size);
// Получение суммы элементов массива 
double ElArrSum(double*& arr, int& size, double& sum);

/// Функции тестирования 
double TestFillarr1(double*& arr, double*& empty1, double*& empty2, int& size)
{
	return Fillarr1(arr, size);
}

double TestFillarr2(double*& empty1, double*& arr, double*& empty2, int& size)
{
	return Fillarr2(arr, size);
}

double TestSum(double*& arr_in1, double*& arr_in2, double*& arr_res, int& size)
{
	return SumArrays(arr_in1, arr_in2, arr_res, size);
}

double TestElArrSum(double*& arr, double*& empty1, double*& empty2, int& size)
{
	double summel = 0;
	return ElArrSum(arr, size, summel);
}
// Реализация функций
double Fillarr1(double* arr, int size)
{
	double t_start = omp_get_wtime();
	for (int i = 0; i < size; i++)
		arr[i] = (sin(i) + 0.5 * cos(i)) * i;
	double t_end = omp_get_wtime();
	return t_end - t_start;
}

double Fillarr2(double* arr, int size)
{
	double t_start = omp_get_wtime();
	for (int i = 0; i < size; i++)
		arr[i] = (cos(i) * sin(i / 2.0)) * i;
	double t_end = omp_get_wtime();
	return t_end - t_start;
}

double SumArrays(double*& arr_in1, double*& arr_in2, double*& arr_res, int& size)
{
	double t_start = omp_get_wtime();
	for (int i = 0; i < size; i++)
		arr_res[i] = arr_in1[i] + arr_in2[i];
	double t_end = omp_get_wtime();
	return t_end - t_start;
}

double ElArrSum(double*& arr, int& size, double& sum)
{
	double t_start = omp_get_wtime();
	double temp_summ = 0;
	for (int i = 0; i < size; i++)
		temp_summ += arr[i];
	double t_end = omp_get_wtime();
	return t_end - t_start;
}

double TestIter(void* Funct, double*& arr_in1, double*& arr_in2, double*& arr_res, int size, int iterations)
{
	double curtime = 0, avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	double* Times = new double[iterations];
	for (int i = 0; i < iterations; i++)
	{
		// Запуск функции и получение времени в миллисекундах
		curtime = ((*(TestFunctTempl)Funct)(arr_in1, arr_in2, arr_res, size)) * 1000;
		// запись времени в массив для определения среднеарифметического значения в доверительном интервале
		Times[i] = curtime;
		avgTime += curtime;
		//       std::cout << curtime << std::endl;
	}
	// Вычисление среднеарифметического по всем итерациям и вывод значения на экран
	avgTime /= iterations;
	std::cout << "AvgTime:" << avgTime << std::endl;
	// Определения среднеарифметического значения в доверительном интервале по всем итерациям и вывод значения на экран
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	std::cout << "AvgTimeTrusted:" << avgTimeT << std::endl;
	return avgTimeT;
}

// Функция рассчета среднеарифметического значения в доверительном интервале
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


void test_functions(void** Functions, string(&function_names)[4], int size, int iterations)
{
	double* a, *b, *c;
	a = new double[size];
	b = new double[size];
	c = new double[size];
	double Functions_time_ms[4];
	// проведение замеров времени работы по каждой функции
	for (int i = 0; i < 4; i++)
		Functions_time_ms[i] = TestIter(Functions[i], a, b, c, size, iterations);

	// Вывод результатов замера (можно организовать вывод в файл)
	for (int i = 0; i < 4; i++)
		cout << function_names[i] << "\t" << Functions_time_ms[i] << " ms." << endl;
	delete[]a;
	delete[]b;
	delete[]c;
}

int main()
{
	void** Functions = new void*[4]{ TestFillarr1 ,TestFillarr2,TestSum,TestElArrSum };
	string  function_names[4]{ "Заполнение вар 1","Заполнение вар 2","Сумма массивов A и B вар 1","Сумма элементов массива С вар 1" };
	setlocale(LC_ALL, "Russian");
	test_functions(Functions, function_names, 100000, 2000);
}