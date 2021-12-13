#include<tbb/tbb.h>
#include <omp.h>
#include<iostream>
using namespace tbb;
using namespace std;

class MyTask
{
public:MyTask(const char* name) : name_(name) {}
	   void operator()() const
	   {
		   std::cout << "Hello from task " << name_ << std::endl;
	   }
private:
	const char* name_;
};



/* ---- Алгоритм вычисления суммы  ---- */

// TBB Reduction

/// Используемый функцией класс
class reduce_total {
	double* m1;
	double* m2;
	int size;

public:
	double sumRes1, sumRes2;

	void operator() (blocked_range<size_t>& r) {
		double sumRes1_local = sumRes1, sumRes2_local = sumRes2;
		double* m1_local = m1;
		double* m2_local = m2;

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

	reduce_total(double* orig_m1, double* orig_m2, int orig_size) : m1(orig_m1), m2(orig_m2), size(orig_size), sumRes1(0.0), sumRes2(0.0) {} // Основной конструктор
};

/// Непосредственно функция
double matrTotalSum_tbbReduct(double**& matrix1, double**& matrix2, double**& unused, int size) {

	double time_start = omp_get_wtime();

	double total1 = 0.0, total2 = 0.0;

	for (int i = 0; i < size; i++) {
		reduce_total r(matrix1[i], matrix2[i], size);
		parallel_reduce(blocked_range<size_t>(0, size), r);
		total1 += r.sumRes1;
		total2 += r.sumRes2;
	}

	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}

/* ---- Алгоритм нахождения максимума  ---- */

// TBB Reduction

// Используемый функцией класс
class reduce_max {
private:
	double* ArrayA;
	double* ArrayB;
public:
	double MaxValue1, MaxValue2;

	void operator() (const blocked_range<size_t>& r) {
		double* a = ArrayA;
		double* b = ArrayB;
		double val1, val2;

		for (size_t i = r.begin(); i != r.end(); ++i) {
			val1 = a[i];
			val2 = b[i];

			if (MaxValue1 < val1) MaxValue1 = val1;
			if (MaxValue2 < val2) MaxValue2 = val2;
		}
	}

	reduce_max(reduce_max& x, tbb::split) : ArrayA(x.ArrayA), ArrayB(x.ArrayB), MaxValue1(-DBL_MAX), MaxValue2(-DBL_MAX) {}
	void join(const reduce_max& y) {
		if (y.MaxValue1 > MaxValue1) MaxValue1 = y.MaxValue1;
		if (y.MaxValue2 > MaxValue2) MaxValue2 = y.MaxValue2;
	}

	reduce_max(double* A, double* B) : ArrayA(A), ArrayB(B), MaxValue1(-DBL_MAX), MaxValue2(-DBL_MAX) {}
};

/// Непосредственно функция
double matrMax_tbbReduct(double**& matrix1, double**& matrix2, double**& unused, int size) {

	double time_start = omp_get_wtime();

	double max1 = -DBL_MAX, max2 = -DBL_MAX;

	for (int i = 0; i < size; i++) {
		reduce_max r(matrix1[i], matrix2[i]);
		parallel_reduce(blocked_range<size_t>(0, size), r);
		if (max1 < r.MaxValue1) max1 = r.MaxValue1;
		if (max2 < r.MaxValue2) max2 = r.MaxValue2;
	}

	double time_stop = omp_get_wtime();
	return time_stop - time_start;
}


void FillMatrix(double**& matrix, int size)
{
	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < size; i++)
		{
			matrix[i][j] = pow(i + 1, 3 / 4)* cos(i) / atan(i + 1) + pow(j + 1, 3 / 4)* cos(j) / atan(j + 1);
		}
	}
}

int main()
{
	int wait;
	//task_group tg;
	//tg.run(MyTask("1")); /* Генерация задачи */
	//tg.run(MyTask("2")); /* Генерация задачи */
	//tg.wait(); /* Ожидание завершения задачи */
	
	int size = 500;

	double** m1 = new double*[size];
	double** m2 = new double*[size];
	double** m;

	for (int i = 0; i < size; i++)
	{
		m1[i] = new double[size];
		m2[i] = new double[size];
	}

	FillMatrix(m1, size);
	FillMatrix(m2, size);

	double time = matrMax_tbbReduct(m1, m2, m, size);
	cout << time;
	std::cin >> wait;
	return 0;
}