#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <String>
#include <ctime>
#include <Windows.h>
#include <Math.h>
#include <fstream>
#include "BMPFileRW.hpp"

#define Random(a, b) (double)rand() / (RAND_MAX + 1)*((b)-(a)) + (a)

using namespace std;

typedef void(*IntSortFunc)(int*, int);
typedef void(*DoubleSortFunc)(double*, int);
typedef void(*ByteSortFunc)(BYTE*, int);

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
	for (int i = 0; i < length; i++)
	{
		for (int j = (i + 1) % 2; j < length - 1; j += 2)
			if (arr[j] > arr[j + 1])
			{
				Swap(arr[j], arr[j + 1]);
			}
	}
}

//Сортировка чет-нечет параллельная с for
template<class T>
void BubbleEvenSortParallelFor(T* arr, int length)
{
	for (int i = 0; i < length; i++)
	{
#pragma omp parallel for shared(i, sorted, narr,compare) 
		for (int j = (i + 1) % 2; j < length - 1; j += 2)
			if (arr[j] > arr[j + 1])
			{
				Swap(arr[j], arr[j + 1]);
			}
	}
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
#pragma omp parallel for shared(narr, compare, step) private(j)
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
void QuickSortParallelSectionsPart(T* arr, int length, int numThreads) {

	int i = 0, j = length - 1;
	T temp, p;
	p = arr[length >> 1];

	do {
#pragma omp parallel sections shared(a, i, j) num_threads(threads)
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

#pragma omp parallel sections shared(j, i, a) num_threads(threads)
	{
#pragma omp section
		{
			if (j > 0) QuickSortParallelSectionsPart(arr, j + 1, numThreads - 1);
		}
#pragma omp section
		{
			if (length > i) QuickSortParallelSectionsPart(arr + i, length - i, numThreads - 2);
		}
	}
}

//Быстрая сортировка распараллеленая секциями
template<class T>
void QuickSortParallelSections(T* arr, int length) {

	int i = 0, j = length - 1;
	T temp, p;
	p = arr[length >> 1];
	int numThreads = 1;

#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
	}

	do {
#pragma omp parallel sections shared(a, i, j) num_threads(threads)
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

#pragma omp parallel sections shared(j, i, a) num_threads(threads)
	{
#pragma omp section
		{
			if (j > 0) QuickSortParallelSectionsPart(arr, j + 1, numThreads - 1);
		}
#pragma omp section
		{
			if (length > i) QuickSortParallelSectionsPart(arr + i, length - i, numThreads - 2);
		}
	}
}

//Заполняет последовательность случайными числами
template<class T>
void GetRandomValues(T *array, int length, int minValue, int maxValue)
{
	for (int i = 0; i < length; i++)
		array[i] = Random(minValue, maxValue);
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

void SortIntFuncAverageTime(void* func, int dataAmount, double& time, int iterations)
{
	int* arr = new int[dataAmount];
	GetRandomValues(arr, dataAmount, -dataAmount, dataAmount);
	double avgTime = 0, avgTimeT = 0, correctAVG = 0;
	double startTime, curTime;
	double* Times = new double[iterations];

	std::cout << "[";
	for (int i = 0; i < iterations; i++)
	{
		GetRandomValues(arr, dataAmount, -dataAmount, dataAmount);
		startTime = omp_get_wtime();
		(*(IntSortFunc)func)(arr, dataAmount);
		curTime = omp_get_wtime() - startTime;
		Times[i] = curTime;
		avgTime += curTime;
		std::cout << "#";
	}
	std::cout << "]\n";
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	time = avgTimeT;
	delete[] arr;
}

void SortDoubleFuncAverageTime(void* func, int dataAmount, double& time, int iterations)
{
	double* arr = new double[dataAmount];
	GetRandomValues(arr, dataAmount, -dataAmount, dataAmount);
	double avgTime = 0, avgTimeT = 0, correctAVG = 0;
	double startTime, curTime;
	double* Times = new double[iterations];

	std::cout << "[";
	for (int i = 0; i < iterations; i++)
	{
		startTime = omp_get_wtime();
		(*(DoubleSortFunc)func)(arr, dataAmount);
		curTime = omp_get_wtime() - startTime;
		Times[i] = curTime;
		avgTime += curTime;
		std::cout << "#";
	}
	std::cout << "]\n";
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	time = avgTimeT;
	delete[] arr;
}

template<class T>
void TaskSort()
{
	int* dataAmount = new int[3]{10000, 25000, 50000};

	double time;
	string* funcsNames = new string[7]{ "Классический алгоритм пузырька", "Чет-нечетной перестановки (последовательный)", "Чет-нечет перестановки (параллельный с for)", "Алгоритм Шелла (последовательный)", "Алгоритм Шелла (параллельный с for)", "Алгоритм быстрой сортировки (последовательный)", "Алгоритм быстрой сортировки (параллельный с sections)" };
	void** sortFuncs = new void*[7]{ BubbleSortConsistently<T>, BubbleEvenSortConsistently<T>, BubbleEvenSortParallelFor<T>, ShellSortConsistently<T>, ShellSortParallelFor<T>, QuickSortConsistently<T>, QuickSortParallelSections<T> };

	std::cout << "Сортировка массивов с типом данных Int" << endl;
	for (int d = 0; d < 4; d++)
	{
		std::cout << "- Количество элементов в массиве: " << dataAmount[d] << endl;
		for (int i = 0; i < 7; i++)
		{
			std::cout << funcsNames[i] << endl;
			for (int t = 2; t < 5; t++)
			{
				omp_set_num_threads(t);
				SortIntFuncAverageTime(sortFuncs[i], dataAmount[d], time, 20);
				std::cout << "- - Потоков: " << t << ". " << time << " с. длилась сортировка" << endl;
				if (i == 0 || i == 1 || i == 3 || i == 5)
					break;
			}
			std::cout << endl;
		}
	}

	std::cout << "Сортировка массивов с типом данных Double" << endl;
	for (int d = 0; d < 4; d++)
	{
		std::cout << "- Количество элементов в массиве: " << dataAmount[d] << endl;
		for (int i = 0; i < 7; i++)
		{
			cout << funcsNames[i] << endl;
			for (int t = 2; t < 5; t++)
			{
				omp_set_num_threads(t);
				SortDoubleFuncAverageTime(sortFuncs[i], dataAmount[d], time, 20);
				std::cout << "- - Потоков: " << t << ". " << time << " с. длилась сортировка" << endl;
				if (i == 0 || i == 1 || i == 3 || i == 5)
					break;
			}
			std::cout << endl;
		}
	}
}

//заполнение медиального массива
//image - исходная картинка
//(x,y) - центр рамки
//RH, RW - радиусы рамки по высоте(height) и ширине(width)
RGBQUAD* getMedial(RGBQUAD **&image, int width, int height, int x, int y, int RH, int RW)
{
	int index = 0;
	RGBQUAD* barray = new RGBQUAD[(2 * RW + 1) * (2 * RH + 1)];
	int coordX;
	int coordY;
	for (int dy = -RH; dy <= RH; dy++)
	{
		coordY = y + dy;
		for (int dx = -RW; dx <= RW; dx++)
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
	//BYTE *red1;
	//BYTE *green1;
	//BYTE *blue1;
	for (int i = 0; i < length; i++)
	{
		red[i] = arr[i].rgbRed;
		green[i] = arr[i].rgbGreen;
		blue[i] = arr[i].rgbBlue;
	}

	(*(ByteSortFunc)sortFunc)(red, length);
	(*(ByteSortFunc)sortFunc)(green, length);
	(*(ByteSortFunc)sortFunc)(blue, length);
	RGBQUAD* narr = new RGBQUAD[length];
	for (int i = 0; i < length; i++)
		narr[i] = { blue[i], green[i],red[i], 0 };
	//delete[] red1;
	//delete[] green1;
	//delete[] blue1;
	delete[] red;
	delete[] green;
	delete[] blue;
	return narr;
}

//сортировка массива РГБ с Omp sections+for
RGBQUAD* sortRGBAsync(RGBQUAD* arr, long length, void* sortFunc)
{
	BYTE *red = new BYTE[length];
	BYTE *green = new BYTE[length];
	BYTE *blue = new BYTE[length];
	//BYTE *red1;
	//BYTE *green1;
	//BYTE *blue1;
#pragma omp parallel for shared(arr, red, blue, green)
	for (int i = 0; i < length; i++)
	{
		red[i] = arr[i].rgbRed;
		green[i] = arr[i].rgbGreen;
		blue[i] = arr[i].rgbBlue;
	}
#pragma omp parallel sections shared(red, blue, green, length, red1, blue1, green1)
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
	//delete[] red1;
	//delete[] green1;
	//delete[] blue1;
	return narr;
}

//медианная фильтрация
//RGB - исходное изображение
//RH, RW - радиусы рамки по вертикали и горизонтали
//method - метод сортировки байтового массива
//RGBresult должен быть инициализирован, в него возвращается результат
void medianFiltering(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult, void* sortFunc)
{
	RGBQUAD *temp1, *temp2;
	int size = (2 * RH + 1) * (2 * RW + 1);
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			//в окне H x W ложу пиксели в массив temp
			temp1 = getMedial(RGB, width, height, x, y, RH, RW); //заполняю медиальный массив
			temp2 = sortRGB(temp1, size, sortFunc); // сортирую каждую из компонент
			RGBresult[y][x] = temp2[size / 2]; // вытаскиваю срединный элемент
			delete[] temp1;
			delete[] temp2;
		}
	}
}

//медианная фильтрация c распараллеливанием сортировки по компонентам
//RGB - исходное изображение
//RH, RW - радиусы рамки по вертикали и горизонтали
//method - метод сортировки байтового массива
//RGBresult должен быть инициализирован, в него возвращается результат
void medianFilteringParallel(RGBQUAD** &RGB, int height, int width, int RH, int RW, RGBQUAD** &RGBresult, void* sortFunc)
{
	RGBQUAD *temp1, *temp2;
	int size = (2 * RH + 1) * (2 * RW + 1);
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			//в окне H x W ложу пиксели в массив temp
			temp1 = getMedial(RGB, width, height, x, y, RH, RW); //заполняю медиальный массив
			temp2 = sortRGBAsync(temp1, size, sortFunc);
			RGBresult[y][x] = temp2[size / 2]; // вытаскиваю срединный элемент
			delete[] temp1;
			delete[] temp2;
		}
	}
}

void TaskMedianFiltering()
{
	void** sortFuncs = new void*[2]{QuickSortParallelSections<BYTE>, ShellSortParallelFor<BYTE>};
	string* sortFuncsNames = new string[2]{ "Быстрая сортировка", "Сортировка Шелла" };
	int* kSize = new int[3]{ 9, 15, 21 };
	RGBQUAD** RGB;
	RGBQUAD** result;
	BITMAPFILEHEADER head;
	BITMAPINFOHEADER info;
	string str = "";
	double time, timestart;

	//
	BMPRead(RGB, head, info, "C:\\Users\\vanya\\Source\\Repos\\ParallelnieVichisleniya\\lab3\\input.bmp");
	result = new RGBQUAD*[info.biHeight];
	for (int i = 0; i < info.biHeight; i++)
		result[i] = new RGBQUAD[info.biWidth];
	medianFiltering(RGB, info.biHeight, info.biWidth, kSize[0], kSize[0], result, sortFuncs[0]);
	BMPWrite(result, head, info, "output.bmp");
	//
	/*
	for (int d = 0; d < 3; d++)
	{
		std::cout << "K size = " << kSize[d] << endl;
		for (int i = 0; i < 2; i++)
		{
			std::cout << sortFuncsNames[i] << endl;
			for (int t = 2; t < 5; t++)
			{
				BMPRead(RGB, head, info, "input.bmp");
				omp_set_num_threads(t);
				timestart = omp_get_wtime();
				medianFilteringParallel(RGB, info.biHeight, info.biWidth, kSize[d], kSize[d], result, sortFuncs[i]);
				time = omp_get_wtime() - timestart;
				str = "output_" + sortFuncsNames[i]; //+ "_" + kSize[d] + "_" + t + ".bmp";
				BMPWrite(result, head, info, str.c_str());
				std::cout << "- - Потоков: " << t << ". " << time << " с. длилась сортировка" << endl;
				//if (i == 0 || i == 1 || i == 3 || i == 5)
				//	break;
			}
			std::cout << endl;
		}
	}
	*/
	//Беру алгоритмы
	//1 реализация - ShellSort.
	//                          Последовательный вариант с medianFiltering + sortRGB.
	//                          Параллельный вариант с medianFilteringAsyncSort + sortRGBAsync
	//2 реализация - QuickSort.
	//                          Последовательный вариант с medianFiltering + sortRGB.
	//                          Параллельный вариант с medianFilteringAsyncSort + sortRGBAsync
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
		TaskSort<int>();
	}
	else
	{
		TaskMedianFiltering();
	}
}