#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <String>
#include <ctime>
#include <Windows.h>
#include <Math.h>
#include <fstream>
#include <sstream>
#include "BMPFileRW.h"

# define PI 3.14159265358979323846

using namespace std;

typedef void(*FilteringFunc)(RGBQUAD**, int, int, int, RGBQUAD**);

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


//Линейный среднеарифметический фильтр последовательный
void LineMiddleFiltering(RGBQUAD** &RGB, int height, int width, int kSize, RGBQUAD** &RGBresult)
{
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			int rgbBlue = 0, rgbGreen = 0, rgbRed = 0;
			for (int DY = -kSize; DY <= kSize; DY++)
			{
				int KY = Y + DY;
				if (KY < 0)
					KY = 0;
				if (KY > height - 1)
					KY = height - 1;
				for (int DX = -kSize; DX <= kSize; DX++)
				{
					int KX = X + DX;
					if (KX < 0)
						KX = 0;
					if (KX > width - 1)
						KX = width - 1;
					rgbBlue += RGB[KY][KX].rgbBlue;
					rgbGreen += RGB[KY][KX].rgbGreen;
					rgbRed += RGB[KY][KX].rgbRed;
				}
			}
			RGBresult[Y][X].rgbBlue = rgbBlue / ((kSize * 2 + 1) * (kSize * 2 + 1));
			RGBresult[Y][X].rgbGreen = rgbGreen / ((kSize * 2 + 1) * (kSize * 2 + 1));
			RGBresult[Y][X].rgbRed = rgbRed / ((kSize * 2 + 1) * (kSize * 2 + 1));
		}
}

//Линейный среднеарифметический фильтр параллельный
void LineMiddleFilteringParallel(RGBQUAD** &RGB, int height, int width, int kSize, RGBQUAD** &RGBresult)
{
#pragma omp parallel for firstprivate(kSize, height, width) shared(RGB, RGBresult)
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			int rgbBlue = 0, rgbGreen = 0, rgbRed = 0;
			for (int DY = -kSize; DY <= kSize; DY++)
			{
				int KY = Y + DY;
				if (KY < 0)
					KY = 0;
				if (KY > height - 1)
					KY = height - 1;
				for (int DX = -kSize; DX <= kSize; DX++)
				{
					int KX = X + DX;
					if (KX < 0)
						KX = 0;
					if (KX > width - 1)
						KX = width - 1;
					rgbBlue += RGB[KY][KX].rgbBlue;
					rgbGreen += RGB[KY][KX].rgbGreen;
					rgbRed += RGB[KY][KX].rgbRed;
				}
			}
			RGBresult[Y][X].rgbBlue = rgbBlue / ((kSize * 2 + 1) * (kSize * 2 + 1));
			RGBresult[Y][X].rgbGreen = rgbGreen / ((kSize * 2 + 1) * (kSize * 2 + 1));
			RGBresult[Y][X].rgbRed = rgbRed / ((kSize * 2 + 1) * (kSize * 2 + 1));
		}
	/*
#pragma omp parallel for
	for (int Y = 0; Y < height; Y++)
		for (int X = 0; X < width; X++)
		{
			int rgbBlue = 0, rgbGreen = 0, rgbRed = 0;
			for (int DY = -kSize; DY <= kSize; DY++)
			{
				int KY = Y + DY;
				if (KY < 0)
					KY = 0;
				if (KY > height - 1)
					KY = height - 1;
				for (int DX = -kSize; DX <= kSize; DX++)
				{
					int KX = X + DX;
					if (KX < 0)
						KX = 0;
					if (KX > width - 1)
						KX = width - 1;
					rgbBlue += RGB[KY][KX].rgbBlue;
					rgbGreen += RGB[KY][KX].rgbGreen;
					rgbRed += RGB[KY][KX].rgbRed;
				}
			}
			RGBresult[Y][X].rgbBlue = rgbBlue / ((kSize * 2 + 1) * (kSize * 2 + 1));
			RGBresult[Y][X].rgbGreen = rgbGreen / ((kSize * 2 + 1) * (kSize * 2 + 1));
			RGBresult[Y][X].rgbRed = rgbRed / ((kSize * 2 + 1) * (kSize * 2 + 1));
		}
	*/
}

//Формирование матрицы коэффициентов для фильтрации Гаусса
double** GaussMatrixCoefficients(int kSize, double q) {
	double** Result = new double*[kSize * 2 + 1];
	for (int i = 0; i < kSize * 2 + 1; i++)
		Result[i] = new double[kSize * 2 + 1];
	double SUM = 0;
	for (int Y = -kSize; Y <= kSize; Y++)
		for (int X = -kSize; X <= kSize; X++) {
			double CF = (1 / (2 * PI * q * q)) * exp(-1 * (X * X + Y * Y) / (2 * q * q));
			Result[Y + kSize][X + kSize] = CF;
			SUM += CF;
		}
	for (int Y = -kSize; Y <= kSize; Y++)
		for (int X = -kSize; X <= kSize; X++)
			Result[Y + kSize][X + kSize] /= SUM;
	return Result;
}

//Формирование матрицы коэффициентов для фильтрации Гаусса параллельное
double** GaussMatrixCoefficientsOMP(int kSize, double q) {
	double** Result = new double*[kSize * 2 + 1];
	for (int i = 0; i < kSize * 2 + 1; i++)
		Result[i] = new double[kSize * 2 + 1];
	double SUM = 0;
#pragma omp parallel for reduction(+:SUM) shared(Result)
	for (int Y = -kSize; Y <= kSize; Y++)
		for (int X = -kSize; X <= kSize; X++) {
			double CF = (1 / (2 * PI * q * q)) * exp(-1 * (X * X + Y * Y) / (2 * q * q));
			Result[Y + kSize][X + kSize] = CF;
			SUM += CF;
		}
#pragma omp parallel for shared(Result)
	for (int Y = -kSize; Y <= kSize; Y++)
		for (int X = -kSize; X <= kSize; X++) {
			Result[Y + kSize][X + kSize] /= SUM;
		}
	return Result;
}

//Линейный фильтр Гаусса последовательный;
//Возвращает RGBresult указывающий на выходную картинку
void LineGaussFiltering(RGBQUAD** &RGB, int height, int width, int kSize, RGBQUAD** &RGBresult)
{
	double** CoefMatrix = GaussMatrixCoefficients(kSize, kSize / 3.0); //Сигма тут
	for (int Y = 0; Y < height; Y++)
	{
		for (int X = 0; X < width; X++)
		{
			double rgbBlue = 0, rgbGreen = 0, rgbRed = 0;
			for (int DY = -kSize; DY <= kSize; DY++)
			{
				int KY = Y + DY;
				if (KY < 0)
					KY = 0;
				if (KY > height - 1)
					KY = height - 1;
				for (int DX = -kSize; DX <= kSize; DX++)
				{
					int KX = X + DX;
					if (KX < 0)
						KX = 0;
					if (KX > width - 1)
						KX = width - 1;
					double tmp = CoefMatrix[DY + kSize][DX + kSize];
					rgbBlue += RGB[KY][KX].rgbBlue * tmp;
					rgbGreen += RGB[KY][KX].rgbGreen * tmp;
					rgbRed += RGB[KY][KX].rgbRed * tmp;
				}
			}
			if (rgbBlue < 0)	rgbBlue = 0;
			if (rgbBlue > 255)	rgbBlue = 255;
			if (rgbGreen < 0)	rgbGreen = 0;
			if (rgbGreen > 255)	rgbGreen = 255;
			if (rgbRed < 0)		rgbRed = 0;
			if (rgbRed > 255)	rgbRed = 255;
			RGBresult[Y][X].rgbBlue = rgbBlue;
			RGBresult[Y][X].rgbGreen = rgbGreen;
			RGBresult[Y][X].rgbRed = rgbRed;
		}
	}
	for (int i = 0; i < kSize; i++)
		delete[] CoefMatrix[i];
	delete[] CoefMatrix;
}

//Линейный фильтр Гаусса параллельеный;
//Возвращает RGBresult указывающий на выходную картинку
void LineGaussFilteringOMP(RGBQUAD** &RGB, int height, int width, int kSize, RGBQUAD** &RGBresult)
{
	double** CoefMatrix = GaussMatrixCoefficientsOMP(kSize, kSize / 3.0); //Сигма тут
#pragma omp parallel for
	for (int Y = 0; Y < height; Y++)
	{
		for (int X = 0; X < width; X++)
		{
			double rgbBlue = 0, rgbGreen = 0, rgbRed = 0;
			for (int DY = -kSize; DY <= kSize; DY++)
			{
				int KY = Y + DY;
				if (KY < 0)
					KY = 0;
				if (KY > height - 1)
					KY = height - 1;
				for (int DX = -kSize; DX <= kSize; DX++)
				{
					int KX = X + DX;
					if (KX < 0)
						KX = 0;
					if (KX > width - 1)
						KX = width - 1;
					double tmp = CoefMatrix[DY + kSize][DX + kSize];
					rgbBlue += RGB[KY][KX].rgbBlue * tmp;
					rgbGreen += RGB[KY][KX].rgbGreen * tmp;
					rgbRed += RGB[KY][KX].rgbRed * tmp;
				}
			}
			if (rgbBlue < 0)	rgbBlue = 0;
			if (rgbBlue > 255)	rgbBlue = 255;
			if (rgbGreen < 0)	rgbGreen = 0;
			if (rgbGreen > 255)	rgbGreen = 255;
			if (rgbRed < 0)		rgbRed = 0;
			if (rgbRed > 255)	rgbRed = 255;
			RGBresult[Y][X].rgbBlue = rgbBlue;
			RGBresult[Y][X].rgbGreen = rgbGreen;
			RGBresult[Y][X].rgbRed = rgbRed;
		}
	}
	for (int i = 0; i < kSize; i++)
		delete[] CoefMatrix[i];
	delete[] CoefMatrix;
}

void FilteringFuncAverageTime(int filtheringMethodIndex, string fileName, string methodFullName, int kSize, double& time, int iterations)
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

	
		switch (filtheringMethodIndex)
		{
		case 0:
		{
			LineMiddleFiltering(sourceImage, info.biHeight, info.biWidth, kSize, resultImage);
			break;
		}
		case 1:
		{
			LineMiddleFilteringParallel(sourceImage, info.biHeight, info.biWidth, kSize, resultImage);
			break;
		}
		case 2:
		{
			LineGaussFiltering(sourceImage, info.biHeight, info.biWidth, kSize, resultImage);
			break;
		}
		default:
		{
			LineGaussFilteringOMP(sourceImage, info.biHeight, info.biWidth, kSize, resultImage);
			break;
		}
		}
		

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

void TaskFilteringMethods()
{
	std::ofstream resultsFile;

	string* filteringFuncsNames = new string[4]{ "Среднеарифметический фильтр (последовательный)", "Среднеарифметический фильтр (параллельный)", "Фильтр Гаусса (последовательный)", "Фильтр Гаусса (параллельный)" };
	void** filteringFuncs = new void*[4]{ LineMiddleFiltering, LineMiddleFilteringParallel, LineGaussFiltering, LineGaussFilteringOMP };
	string* inputFiles = new string[4]{ "500x500.bmp", "720x480.bmp", "1600x1200.bmp", "1920x1080.bmp"};
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

	for (int j = 0; j < 4; j++)
	{
		resultsFile << filteringFuncsNames[j] << ";";
		cout << "___" << filteringFuncsNames[j] << "___" << endl;
		for (int t = 2; t < 5; t++)
		{
			if (j == 0 || j == 2)
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

				resultsFile << kSize[k] * 2 + 1 << ";";
				std::cout << "Ksize = " << kSize[k] * 2 + 1 << endl;
				for (int d = 0; d < 4; d++)
				{
					ss = stringstream();
					ss << inputFiles[d] << "_ouput_" << filteringFuncsNames[j] << "_k" << kSize[k] << ".bmp";
					std::cout << "Input file: " << inputFiles[d] << endl;
					FilteringFuncAverageTime(j, inputFiles[d], ss.str(), kSize[k], time, 20);
					if (j == 0 || j == 2)
						T1[k][d] = time;
					std::cout << "Длительность: " << time << " с." << endl;
					resultsFile << time << "^;" << (double)T1[k][d] / time << "^;";
				}
				resultsFile << endl;
			}

			if (j == 0 || j == 2)
				break;
		}
	}
}

void main()
{
	srand(time(0));
	setlocale(LC_ALL, "Russian");
	int choice;
	TaskFilteringMethods();
}