#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <String>
#include <ctime>
#include <Windows.h>
#include <Math.h>
#include <fstream>
#include <sstream>
#include "BMPFileRW.hpp"

#define Random(a, b) (double)rand() / (RAND_MAX + 1)*((b)-(a)) + (a)

using namespace std;

typedef void(*IntSortFunc)(int*, int);
typedef void(*DoubleSortFunc)(double*, int);
typedef void(*ByteSortFunc)(BYTE*, int);

typedef void(*FilteringFunc)(RGBQUAD**, int, int, int, RGBQUAD**, void*);

//Обмен значениями
template<class T>
void Swap(T& one, T& two)
{
	T temp = one;
	one = two;
	two = temp;
}

template<class T>
void Print10(T* arr)
{
	for (int i = 0; i < 10; i++)
	{
		cout << arr[i] << " ";
	}
	cout << endl;
}

//Сортировка пузырьком
template<class T>
void BubbleSortConsistently(T* arr, int length)
{
	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < length - i - 1; j++)
			if (arr[j] > arr[j + 1])
			{
				Swap(arr[j], arr[j + 1]);
			}
	}
}

//Сортировка чет-нечет
template<class T>
void BubbleEvenSortConsistently(T* arr, int length)
{
	int res;
	do {
		res = 0;
		for (int i = 0; i < length - 1; i += 2)
		{
			if (arr[i] > arr[i + 1])
			{
				Swap(arr[i], arr[i + 1]);
				res = 1;
			}
		}
		for (int i = 1; i < length - 1; i += 2)
		{
			if (arr[i] > arr[i + 1])
			{
				Swap(arr[i], arr[i + 1]);
				res = 1;
			}
		}
	} while (res != 0);
}

//Сортировка чет-нечет параллельная с for
template<class T>
void BubbleEvenSortParallelFor(T* arr, int length)
{
	int res;
	do {
		res = 0;
#pragma omp parallel for reduction(+:res)
		for (int i = 0; i < length - 1; i += 2)
		{
			if (arr[i] > arr[i + 1])
			{
				Swap(arr[i], arr[i + 1]);
				res++;
			}
		}
#pragma omp parallel for reduction(+:res)
		for (int i = 1; i < length - 1; i += 2)
		{
			if (arr[i] > arr[i + 1])
			{
				Swap(arr[i], arr[i + 1]);
				res++;
			}
		}
	} while (res != 0);
}

//Сортировка Шелла
template<class T>
void ShellSortConsistently(T *arr, int length)
{
	int step, i, j;

	for (step = length / 2; step > 0; step /= 2)
	{
		// Перечисление элементов, которые сортируются на определённом шаге
		for (i = step; i < length; i++)
		{
			// Перестановка элементов внутри подсписка, пока i-тый не будет отсортирован
			for (j = i - step; j >= 0 && arr[j] > arr[j + step]; j -= step)
			{
				Swap(arr[j], arr[j + step]);
			}
		}
	}
}

//Сортировка Шелла параллельная с for
template<class T>
void ShellSortParallelFor(T *arr, int length)
{
	int step, i, j;

	for (step = length / 2; step > 0; step /= 2)
	{
		// Перечисление элементов, которые сортируются на определённом шаге
#pragma omp parallel for shared(step) private(i, j)
		for (i = step; i < length; i++)
		{
			// Перестановка элементов внутри подсписка, пока i-тый не будет отсортирован
			for (j = i - step; j >= 0 && arr[j] > arr[j + step]; j -= step)
			{
				Swap(arr[j], arr[j + step]);
			}
		}
	}
}

//Быстрая сортировка
template<class T>
void QuickSortConsistently(T* arr, int length) {
	long i = 0, j = length - 1;
	T  p;
	p = arr[length >> 1];
	do {
		while (arr[i] < p) i++;
		while (arr[j] > p) j--;

		if (i <= j) {
			Swap(arr[i], arr[j]);
			i++; j--;
		}
	} while (i < j);
	if (j > 0)
		QuickSortConsistently(arr, j + 1);
	if (length > i)
		QuickSortConsistently(arr + i, length - i);
}

//Быстрая сортировка распараллеленая секциями
template<class T>
void QuickSortParallelSectionsPart(T* arr, int length) {

	int i = 0, j = length - 1;
	T temp, p;
	p = arr[length >> 1];

	do {
		while (arr[i] < p) i++;
		while (arr[j] > p) j--;
		if (i <= j) {
			Swap(arr[i], arr[j]);
			i++; j--;
		}
	} while (i < j);

#pragma omp parallel sections
	{
#pragma omp section
		{
			if (j > 0) QuickSortParallelSectionsPart(arr, j + 1);
		}
#pragma omp section
		{
			if (length > i) QuickSortParallelSectionsPart(arr + i, length - i);
		}
	}
}

//Быстрая сортировка распараллеленая секциями
template<class T>
void QuickSortParallelSections(T* arr, int length) {

	int i = 0, j = length - 1;
	T temp, p;
	p = arr[length >> 1];

	do {
#pragma omp parallel sections
		{
#pragma omp section
			{ while (arr[i] < p) i++; }
#pragma omp section
			{ while (arr[j] > p) j--; }

		}
		if (i <= j) {
			Swap(arr[i], arr[j]);
			i++; j--;
		}
	} while (i < j);

#pragma omp parallel sections
	{
#pragma omp section
		{
			if (j > 0) QuickSortParallelSectionsPart(arr, j + 1);
		}
#pragma omp section
		{
			if (length > i) QuickSortParallelSectionsPart(arr + i, length - i);
		}
	}
}

//Заполняет последовательность случайными числами
template<class T>
void GetValues(T *array, int length)
{
	for (int i = 0; i < length; i++)
	{
		array[i] = length * pow(i + 1, 3 / 4)* cos(i) / atan(i + 1);
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

void SortIntFuncAverageTime(IntSortFunc func, int dataAmount, double& time, int iterations)
{
	int* arr = new int[dataAmount];
	GetValues(arr, dataAmount);
	double avgTime = 0, avgTimeT = 0, correctAVG = 0;
	double startTime, curTime;
	double* Times = new double[iterations];

	std::cout << "[";
	for (int i = 0; i < iterations; i++)
	{
		GetValues(arr, dataAmount);
		startTime = omp_get_wtime();
		func(arr, dataAmount);
		curTime = omp_get_wtime() - startTime;
		Times[i] = curTime;
		avgTime += curTime;
		std::cout << "#";
	}
	std::cout << "]\n";
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	time = avgTimeT * 1000;
	delete[] arr;
}

void SortDoubleFuncAverageTime(DoubleSortFunc func, int dataAmount, double& time, int iterations)
{
	double* arr = new double[dataAmount];
	GetValues(arr, dataAmount);
	double avgTime = 0, avgTimeT = 0, correctAVG = 0;
	double startTime, curTime;
	double* Times = new double[iterations];

	std::cout << "[";
	for (int i = 0; i < iterations; i++)
	{
		startTime = omp_get_wtime();
		func(arr, dataAmount);
		curTime = omp_get_wtime() - startTime;
		Times[i] = curTime;
		avgTime += curTime;
		std::cout << "#";
	}
	std::cout << "]\n";
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	time = avgTimeT * 1000;
	delete[] arr;
}

void TaskSort()
{
	std::ofstream resultsFile;
	int* dataAmount = new int[4]{5000, 15000, 35000, 50000};

	double time;
	string* funcsNames = new string[7]{ "Классический алгоритм пузырька", "Чет-нечетной перестановки (последовательный)", "Чет-нечет перестановки (параллельный с for)", "Алгоритм Шелла (последовательный)", "Алгоритм Шелла (параллельный с for)", "Алгоритм быстрой сортировки (последовательный)", "Алгоритм быстрой сортировки (параллельный с sections)" };
	IntSortFunc* sortFuncs = new IntSortFunc[7]{ BubbleSortConsistently<int>, BubbleEvenSortConsistently<int>, BubbleEvenSortParallelFor<int>, ShellSortConsistently<int>, ShellSortParallelFor<int>, QuickSortConsistently<int>, QuickSortParallelSections<int> };
	DoubleSortFunc* sortFuncsDouble = new DoubleSortFunc[7]{ BubbleSortConsistently<double>, BubbleEvenSortConsistently<double>, BubbleEvenSortParallelFor<double>, ShellSortConsistently<double>, ShellSortParallelFor<double>, QuickSortConsistently<double>, QuickSortParallelSections<double> };
	double* T1 = new double[3];
	resultsFile.open("Task1Results.csv", std::ios_base::app);
	resultsFile << "Функция сортировки;Потоки;Время;Sp(n);Ep(n);Время;Sp(n);Ep(n);Время;Sp(n);Ep(n);\n";

	for (int j = 0; j < 2; j++)
	{
		if (j == 0)
		{
			resultsFile << "Тип данных Int;;";
			std::cout << "Сортировка массивов с типом данных Int" << endl;
		}
		else
		{
			resultsFile << "Тип данных Double;;";
			std::cout << "Сортировка массивов с типом данных Double" << endl;
		}
		
		for (int cr = 0; cr < 4; cr++)
			resultsFile << "НД" << cr + 1 << ": " << dataAmount[cr] << ";;;";

		resultsFile << endl;

		for (int i = 0; i < 7; i++)
		{
			resultsFile << funcsNames[i];
			std::cout << funcsNames[i] << endl;
			for (int t = 2; t < 5; t++)
			{
				if (i == 0 || i == 1 || i == 3 || i == 5)
					t = 1;
					
				std::cout << "Потоков: " << t << endl;
				resultsFile << ";" << t << ";";
				omp_set_num_threads(t);

				for (int d = 0; d < 4; d++) //4
				{
					std::cout << "Количество элементов: " << dataAmount[d] << endl;

					if (j == 0)
						SortIntFuncAverageTime(sortFuncs[i], dataAmount[d], time, 20);
					else
						SortDoubleFuncAverageTime(sortFuncsDouble[i], dataAmount[d], time, 20);

					if (i == 0 || i == 1 || i == 3 || i == 5)
						T1[d] = time;
					resultsFile << time << ";" << T1[d] / time << "^;" << T1[d] / (t * time) << "^;";
					std::cout << " - Длительность сортировки: " << time << endl;
				}
				std::cout << endl;
				resultsFile << endl;

				if (i == 0 || i == 1 || i == 3 || i == 5)
					break;
			}
		}
		resultsFile << endl;
	}
}

//заполнение медиального массива
RGBQUAD* getMedial(RGBQUAD **&image, int width, int height, int x, int y, int kSize)
{
	int index = 0;
	RGBQUAD* barray = new RGBQUAD[(2 * kSize + 1) * (2 * kSize + 1)];
	int coordX;
	int coordY;
	for (int dy = -kSize; dy <= kSize; dy++)
	{
		coordY = y + dy;
		for (int dx = -kSize; dx <= kSize; dx++)
		{
			coordX = x + dx;
			if (coordX < 0)
				coordX = 0;

			if (coordX >= width)
				coordX = width - 1;

			if (coordY < 0)
				coordY = 0;

			if (coordY >= height)
				coordY = height - 1;
			barray[index] = image[coordY][coordX];
			index++;
		}
	}
	return barray;
}

//Сортировка массива РГБ
RGBQUAD* sortRGB(RGBQUAD* arr, long length, void* sortFunc)
{
	BYTE *red = new BYTE[length];
	BYTE *green = new BYTE[length];
	BYTE *blue = new BYTE[length];

	for (int i = 0; i < length; i++)
	{
		red[i] = arr[i].rgbRed;
		green[i] = arr[i].rgbGreen;
		blue[i] = arr[i].rgbBlue;
	}

	(*(ByteSortFunc)sortFunc)(red, length);
	(*(ByteSortFunc)sortFunc)(green, length);
	(*(ByteSortFunc)sortFunc)(blue, length);
	RGBQUAD* resultRGBArr = new RGBQUAD[length];
	for (int i = 0; i < length; i++)
		resultRGBArr[i] = { blue[i], green[i],red[i], 0 };
	delete[] red;
	delete[] green;
	delete[] blue;
	return resultRGBArr;
}

//сортировка массива РГБ с Omp sections и for
RGBQUAD* sortRGB_OMP(RGBQUAD* arr, long length, void* sortFunc)
{
	BYTE *red = new BYTE[length];
	BYTE *green = new BYTE[length];
	BYTE *blue = new BYTE[length];
#pragma omp parallel for
	for (int i = 0; i < length; i++)
	{
		red[i] = arr[i].rgbRed;
		green[i] = arr[i].rgbGreen;
		blue[i] = arr[i].rgbBlue;
	}

#pragma omp parallel sections
	{
#pragma omp section
		{
			(*(ByteSortFunc)sortFunc)(red, length);
		}
#pragma omp section
		{
			(*(ByteSortFunc)sortFunc)(green, length);
		}
#pragma omp section
		{
			(*(ByteSortFunc)sortFunc)(blue, length);
		}
	}
	RGBQUAD* narr = new RGBQUAD[length];
	for (int i = 0; i < length; i++)
		narr[i] = { blue[i], green[i],red[i],  0 };
	delete[] red;
	delete[] green;
	delete[] blue;
	return narr;
}

//медианная фильтрация
void MedianFiltering(RGBQUAD** &RGB, int height, int width, int kSize, RGBQUAD** &RGBresult, void* sortFunc)
{
	RGBQUAD *temp1, *temp2;
	int size = (2 * kSize + 1) * (2 * kSize + 1);
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			//в окне H x W кладу пиксели в массив temp
			temp1 = getMedial(RGB, width, height, x, y, kSize); //заполняю медиальный массив
			temp2 = sortRGB(temp1, size, sortFunc); // сортирую каждую из компонент
			RGBresult[y][x] = temp2[size / 2]; // вытаскиваю срединный элемент
			delete[] temp1;
			delete[] temp2;
		}
	}
}

//медианная фильтрация (параллельная)
void MedianFilteringOMP(RGBQUAD** &RGB, int height, int width, int kSize, RGBQUAD** &RGBresult, void* sortFunc)
{
	int size = (2 * kSize + 1) * (2 * kSize + 1);
#pragma omp parallel for
	for (int y = 0; y < height; y++)
	{
		RGBQUAD *temp1, *temp2;
		for (int x = 0; x < width; x++)
		{
			//в окне H x W ложу пиксели в массив temp
			temp1 = getMedial(RGB, width, height, x, y, kSize); //заполняю медиальный массив
			temp2 = sortRGB_OMP(temp1, size, sortFunc);
			RGBresult[y][x] = temp2[size / 2]; // вытаскиваю срединный элемент
			delete[] temp1;
			delete[] temp2;
		}
	}
}

void FilteringFuncAverageTime(int filtheringMethod, void* sortFunc, string fileName, string methodFullName, int kSize, double& time, int iterations)
{

	double avgTime = 0, avgTimeT = 0, correctAVG = 0;
	double startTime, curTime;
	double* Times = new double[iterations];

	std::cout << "[";
	for (int i = 0; i < iterations; i++)
	{
		RGBQUAD** sourceImage;
		RGBQUAD** resultImage;
		BITMAPFILEHEADER head;
		BITMAPINFOHEADER info;
		BMPRead(sourceImage, head, info, fileName.c_str());
		resultImage = new RGBQUAD*[info.biHeight];
		for (int i = 0; i < info.biHeight; i++)
			resultImage[i] = new RGBQUAD[info.biWidth];
		startTime = omp_get_wtime();

		if (filtheringMethod == 0)
			MedianFiltering(sourceImage, info.biHeight, info.biWidth, kSize, resultImage, sortFunc);
		else
			MedianFilteringOMP(sourceImage, info.biHeight, info.biWidth, kSize, resultImage, sortFunc);

		curTime = omp_get_wtime() - startTime;
		Times[i] = curTime;
		avgTime += curTime;
		std::cout << "#";

		if (i == 0)
		{
			BMPWrite(resultImage, head, info, methodFullName.c_str());
		}

		for (int i = 0; i < info.biHeight; i++)
		{
			delete[] sourceImage[i];
		}
		delete[] sourceImage;

		for (int i = 0; i < info.biHeight; i++)
		{
			delete[] resultImage[i];
		}
		delete[] resultImage;
	}
	std::cout << "]\n";
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	time = avgTimeT;
}

void TaskMedianFiltering()
{
	std::ofstream resultsFile;

	void** sortFuncs = new void*[2]{ QuickSortConsistently<BYTE>, ShellSortConsistently<BYTE>};
	string* sortFuncsNames = new string[2]{ "Быстрая сортировка (последовательная)", "Сортировка Шелла (последовательная)" };
	string* filteringFuncsNames = new string[2]{"Последовательная медианная фильтрация", "Параллельная медианная фильтрация"};
	string* inputFiles = new string[4]{"500x500.bmp", "640x420.bmp", "840x480.bmp", "1280x720.bmp"};
	int* kSize = new int[3]{ 3, 9, 12 };
	int iterations = 20;
	stringstream ss;
	double time;

	double** T1 = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		T1[i] = new double[4];
	}
	resultsFile.open("Task2Results.csv", std::ios_base::app);
	resultsFile << "Функция;Потоки;Ksize;НД1;;НД2;;НД3;;НД4\n";
	resultsFile << ";;;Время;Sp(n);Время;Sp(n);Время;Sp(n);Время;Sp(n);\n";

	for (int i = 1; i < 2; i++)
	{
		std::cout << sortFuncsNames[i] << endl;
		for (int j = 1; j < 2; j++)
		{
			resultsFile << sortFuncsNames[i] << " [" << filteringFuncsNames[j] << "];";
			cout << "___" << filteringFuncsNames[j] << "___" << endl;
			for (int t = 2; t < 5; t++)
			{
				if (j == 0)
					t = 1;

				omp_set_num_threads(t);

				if (t > 2)
					resultsFile << ";";

				resultsFile << t << ";";
				std::cout << "Потоков: " << t << endl;
				for (int k = 0; k < 3; k++)
				{
					if (k != 0)
						resultsFile << ";;";

					resultsFile << kSize[k] << ";";
					std::cout << "Ksize = " << kSize[k] << endl;
					for (int d = 0; d < 4; d++)
					{
						ss = stringstream();
						ss << inputFiles[d] << "_ouput_" << filteringFuncsNames[j] << "__" << sortFuncsNames[i] << "_k" << kSize[k] << ".bmp";
						std::cout << "Input file: " << inputFiles[d] << endl;
						FilteringFuncAverageTime(j, sortFuncs[i], inputFiles[d], ss.str(), kSize[k], time, 20);
						if (j == 0)
							T1[k][d] = time;
						std::cout << ". Длительность: " << time << " с." << endl;
						resultsFile << time << "^;" << (double)T1[k][d] / time << "^;";
					}
					resultsFile << endl;
				}

				if (j == 0)
					break;
			}
		}
	}
}

void main()
{
	srand(time(0));
	setlocale(LC_ALL, "Russian");
	int choice;

	std::cout << "[1]: Задание 1\n"
		<< "[2]: Задание 2\n";
	cin >> choice;
	if (choice == 1)
	{
		TaskSort();
	}
	else
	{
		TaskMedianFiltering();
	}
}