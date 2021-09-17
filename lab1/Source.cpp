#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <String>
#include <ctime>
#include <Math.h>
#include <fstream>
#include <list>

using namespace std;

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

void Test(int*& a)
{
	a = a + 1;
}

void FillConsistently(Vector& vector, double& time)
{
	double time_Start = omp_get_wtime();
	for (int i = 0; i < vector.length; i++)
	{
		int x = 1 + rand() % 10;
		vector.numbers[i] = pow(x, 2 / 3) * (cos(x) / atan(x));
	}
	time = omp_get_wtime() - time_Start;
}

void FillParallel(Vector& vector, double& time)
{
	double time_Start = omp_get_wtime();
	
#pragma omp parallel
	{
		int currentThread = omp_get_thread_num();
		int startIndex = 0;
		int numThreads = omp_get_num_threads();
		int part = (float)vector.length / numThreads;
		if (currentThread != 0)
		{
			startIndex = part * currentThread;
		}

		if (currentThread == numThreads - 1)
		{
			part = vector.length - (part * currentThread);
		}

		for (int i = 0; i < part; i++)
		{
			int x = 1 + rand() % 10;
			vector.numbers[startIndex + i] = pow(x, 2 / 3) * (cos(x) / atan(x));
			//cout << "Заполение " << startIndex + i << " ячейки\n";
		}
	}

	time = omp_get_wtime() - time_Start;
}

void VectorsSumConsistently(Vector vector1, Vector vector2, Vector& resultVector, double& time)
{
	double time_Start = omp_get_wtime();
	for (int i = 0; i < vector1.length; i++)
	{
		resultVector.numbers[i] = vector1.numbers[i] + vector2.numbers[i];
	}
	time = omp_get_wtime() - time_Start;
}

void VectorsSumParallelFor(Vector vector1, Vector vector2, Vector& resultVector, double& time)
{
	double time_Start = omp_get_wtime();
#pragma omp for
	for (int i = 0; i < vector1.length; i++)
	{
		resultVector.numbers[i] = vector1.numbers[i] + vector2.numbers[i];
	}
	time = omp_get_wtime() - time_Start;
}

void VectorsSumParallelSections(Vector vector1, Vector vector2, Vector& resultVector, double& time, int numThreads)
{
	int s1, s2, s3;
	s1 = vector1.length / numThreads;
	s2 = vector1.length;
	s3 = vector1.length;

	if (numThreads == 3)
	{
		s2 = vector1.length / numThreads * 2;
	}
	else if (numThreads == 4)
	{
		s3 = vector1.length / numThreads * 3;
	}
	
	double time_Start = omp_get_wtime();
#pragma omp sections
	{
#pragma omp section
		{
			for (int i = 0; i < s1; i++)
			{
				resultVector.numbers[i] = vector1.numbers[i] + vector2.numbers[i];
			}
		}
#pragma omp section
		{
			for (int i = s1; i < s2; i++)
			{
				resultVector.numbers[i] = vector1.numbers[i] + vector2.numbers[i];
			}
		}
#pragma omp section
		{
			for (int i = s2; i < s3; i++)
			{
				resultVector.numbers[i] = vector1.numbers[i] + vector1.numbers[i];
			}
		}
#pragma omp section
		{
			for (int i = s3; i < vector1.length; i++)
			{
				resultVector.numbers[i] = vector1.numbers[i] + vector2.numbers[i];
			}
		}
	}
	time = omp_get_wtime() - time_Start;
}

void MagnitudeVectorConsistently(Vector vector, double& sum, double& time)
{
	sum = 0;
	double time_Start = omp_get_wtime();
	for (int i = 0; i < vector.length; i++)
	{
		sum += vector.numbers[i];
	}
	time = omp_get_wtime() - time_Start;
}

void MagnitudeVectorReductor(Vector vector, double& sum, double& time)
{
	double result = 0;
	double time_Start = omp_get_wtime();

#pragma omp parallel
#pragma omp for reduction(+:result)
		for (int i = 0; i < vector.length; i++)
		{
			result += vector.numbers[i];
		}

	sum = result;
	time = omp_get_wtime() - time_Start;
}

void MagnitudeVectorCritical(Vector vector, double& sum, double& time)
{
	sum = 0;
	double time_Start = omp_get_wtime();
#pragma omp critical
	{

		for (int i = 0; i < vector.length; i++)
		{
			sum += vector.numbers[i];
		}
	}
	time = omp_get_wtime() - time_Start;
}

void Task2()
{
	int nVectors = 2;
	int nNumbers = 1000000;
	int iterations = 20;
	double time = 0;
	Vector* vectors = new Vector[nVectors];
	//Инициализация
	for (int i = 0; i < nVectors; i++)
	{
		vectors[i] = Vector(nNumbers);
	}

	double avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	double* Times = new double[iterations];
	for (int i = 0; i < iterations; i++)
	{
		FillConsistently(vectors[0], time);
		Times[i] = time;
		avgTime += time;
	}
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	cout << "Заполнение вектора. Последовательно. " << (int)(avgTimeT * 1000) << endl;


	avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	Times = new double[iterations];
	for (int i = 0; i < iterations; i++)
	{
		FillParallel(vectors[1], time);
		Times[i] = time;
		avgTime += time;
	}
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	cout << "Заполнение вектора. Параллельно. " << (int)(avgTimeT * 1000) << endl;

	Vector vectorC1(nNumbers);
	Vector vectorC2(nNumbers);
	Vector vectorC3(nNumbers);

	avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	Times = new double[iterations];
	for (int i = 0; i < iterations; i++)
	{
		VectorsSumConsistently(vectors[0], vectors[1], vectorC1, time);
		Times[i] = time;
		avgTime += time;
	}
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	cout << "Сложение векторов. Последовательно. " << (int)(avgTimeT * 1000) << endl;
	
	avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	Times = new double[iterations];
	for (int i = 0; i < iterations; i++)
	{
		VectorsSumParallelFor(vectors[0], vectors[1], vectorC2, time);
		Times[i] = time;
		avgTime += time;
	}
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	cout << "Сложение векторов. Параллельный FOR. " << (int)(avgTimeT * 1000) << endl;
	
	avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	Times = new double[iterations];
	for (int i = 0; i < iterations; i++)
	{
		VectorsSumParallelSections(vectors[0], vectors[1], vectorC2, time, omp_get_num_threads());
		Times[i] = time;
		avgTime += time;
	}
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	cout << "Сложение векторов. Параллельный Sections. " << (int)(avgTimeT * 1000) << endl;

	VectorsSumParallelSections(vectors[0], vectors[1], vectorC3, time, omp_get_num_threads());

	double sumC1 = 0;
	double sumC2 = 0;
	double sumC3 = 0;


	avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	Times = new double[iterations];
	for (int i = 0; i < iterations; i++)
	{
		MagnitudeVectorConsistently(vectorC1, sumC1, time);
		Times[i] = time;
		avgTime += time;
	}
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	cout << "Суммирование итогового вектора. Последовательно. " << (int)(avgTimeT * 1000) << endl;

	avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	Times = new double[iterations];
	for (int i = 0; i < iterations; i++)
	{
		MagnitudeVectorReductor(vectorC2, sumC2, time);
		Times[i] = time;
		avgTime += time;
	}
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	cout << "Суммирование итогового вектора. Параллельно с редукторами. " << (int)(avgTimeT * 1000) << endl;

	avgTime = 0, avgTimeT = 0, correctAVG = 0;;
	Times = new double[iterations];
	for (int i = 0; i < iterations; i++)
	{
		MagnitudeVectorCritical(vectorC3, sumC3, time);
		Times[i] = time;
		avgTime += time;
	}
	avgTime /= iterations;
	avgTimeT = AvgTrustedInterval(avgTime, Times, iterations);
	cout << "Суммирование итогового вектора. Параллельно с критическими секциями. " << (int)(avgTimeT * 1000) << endl;
}

struct CalcInfo
{
	int dataAmount;
	double time;

	CalcInfo()
	{
		dataAmount = 0;
		time = 0;
	}

	CalcInfo(double _time, int _dataAmount)
	{
		dataAmount = _dataAmount;
		time = _time;
	}
};

class CalcPerThread
{
	list<CalcInfo> calcInfos;

public:

	int threadsNum;

	CalcPerThread()
	{
		threadsNum = -1;
		calcInfos = list<CalcInfo>();
	}

	CalcPerThread(int _threadsNum, double _time, int _dataAmount)
	{
		threadsNum = _threadsNum;
		calcInfos = list<CalcInfo>();
		Add(_time, _dataAmount);
	}

	bool Find(int _dataAmount)
	{
		for (CalcInfo elem : calcInfos)
		{
			if (elem.dataAmount == _dataAmount)
				return true;
		}

		return false;
	}

	void Add(double _time, int _dataAmount)
	{
		if (!Find(_dataAmount))
		{
			CalcInfo newElem(_time, _dataAmount);
			calcInfos.push_back(newElem);
		}
	}

	void PrintElem(string fileName)
	{
		std::ofstream resultsFile;
		resultsFile.open(fileName, std::ios_base::app);
		resultsFile << ";" << threadsNum;
		for (CalcInfo elem : calcInfos)
		{
			int cTime = (int)(elem.time * 1000);
			resultsFile << ";" << cTime;// << "Ж" << elem.dataAmount;
		}

		resultsFile << endl;
		resultsFile.close();
	}
};

class UnitMethod
{
public:
	string name;
	list<CalcPerThread> threadCalcs;

public:
	UnitMethod()
	{
		name = "";
		threadCalcs = list<CalcPerThread>();
	}

	UnitMethod(string _name)
	{
		name = _name;
		threadCalcs = list<CalcPerThread>();
	}

	UnitMethod(string _name, int _threadsNum, double _time, int _dataAmount)
	{
		name = _name;
		threadCalcs = list<CalcPerThread>();
		Add(_threadsNum, _time, _dataAmount);
	}

	int Find(int _threadsNum)
	{
		int i = 0;
		for (CalcPerThread elem : threadCalcs)
		{
			if (elem.threadsNum == _threadsNum)
			{
				return i;
			}

			i++;
		}

		return -1;
	}

	void Add(int _threadsNum, double _time, int _dataAmount)
	{
		int i = Find(_threadsNum);
		if (i == -1)
		{
			CalcPerThread newCalcPerThread(_threadsNum, _time, _dataAmount);
			threadCalcs.push_back(newCalcPerThread);
		}
		else
		{
			list <CalcPerThread> ::iterator it = threadCalcs.begin();
			advance(it, i);
			it->Add(_time, _dataAmount);
		}
	}

	void PrintElem(string fileName)
	{
		std::ofstream resultsFile;
		resultsFile.open(fileName, std::ios_base::app);
		resultsFile << name;
		resultsFile.close();

		for (CalcPerThread elem : threadCalcs)
		{
			elem.PrintElem(fileName);
		}
	}
};

class CalcData
{
	list<UnitMethod> methods;

public:
	CalcData()
	{
		methods = list<UnitMethod>();
	}

	int Find(string name)
	{
		int i = 0;
		for(UnitMethod elem : methods)
		{
			if (elem.name == name)
				return i;
			i++;
		}

		return -1;
	}

	void Add(string name, int _numThreads, double _time, int _dataAmount)
	{
		int i = Find(name);
		if (i == -1)
		{
			UnitMethod newMethod(name, _numThreads, _time, _dataAmount);
			methods.push_back(newMethod);
		}
		else
		{
			list <UnitMethod> ::iterator it = methods.begin();
			advance(it, i);
			it->Add(_numThreads, _time, _dataAmount);
		}
	}

	void PrintResults(int* numbersAmount, string fileName)
	{
		std::ofstream resultsFile;
		resultsFile.open(fileName, std::ios_base::out);
		resultsFile << "функция;Кол-во потоков;НД " << numbersAmount[0] << ";НД2 " << numbersAmount[1] << ";НД3 " << numbersAmount[2] << ";НД4 "<< numbersAmount[3] << endl;
		resultsFile.close();
		for (UnitMethod elem : methods)
		{
			elem.PrintElem(fileName);
		}
	}
};

void PrintResult(string name1, string name2, string name3, double time, int nA, int nT, string fileName)
{
	cout << "[" << name1 << "] [" << name2 << "] [" << name3 << "]. Потоков: " << nT << " Времени потрачено: " << (int)(time * 1000) << endl;
	
	std::ofstream resultsFile;
	resultsFile.open(fileName, std::ios_base::app);
	resultsFile << name1 << ";" << name2 << ";" << name3 << ";" << nT << ";" << (int)(time * 1000) << endl;
	resultsFile.close();
}

void Task3()
{
	CalcData calcData;
	string fileName = "results.csv";
	std::ofstream resultsFile;
	resultsFile.open(fileName, std::ios_base::out);
	resultsFile << "Способ заполнения;Способ суммирования;Способ подсчета итогового вектора;Количество потоков;Время работы\n";
	resultsFile.close();
	int* numbersAmount = new int[4]{750000, 1000000, 1250000, 1500000};
	double time;
	double time1, time2;
	string name1, name2, name3;
	int iterations = 20; // минимум 20
	double trustedAverageTime1 = 0;
	double trustedAverageTime2 = 0;
	double trustedAverageTime3 = 0;

	double* timeCalculations1 = new double[iterations];
	double* timeCalculations2 = new double[iterations];
	double* timeCalculations3 = new double[iterations];

	for (int nA = 0; nA < 4; nA++)
	{
		for (int nT = 2; nT < 5; nT++)
		{
			omp_set_num_threads(nT);

			for (int j = 0; j < 2; j++)
			{
				Vector vector1(numbersAmount[nA]);
				Vector vector2(numbersAmount[nA]);
				double averageTime = 0;
				for (int iter = 0; iter < iterations; iter++)
				{
					vector1 = Vector(numbersAmount[nA]);
					vector2 = Vector(numbersAmount[nA]);

					switch (j)
					{
						case 0:
						{
							name1 = "Последоватьное заполнение";
							FillConsistently(vector1, time);
							FillConsistently(vector2, time);

							break;
						}
						default:
						{
							name1 = "Параллельное заполнение";
							FillParallel(vector1, time);
							FillParallel(vector2, time);

							break;
						}
					}

					timeCalculations1[iter] = time;
					averageTime += time;
				}
				averageTime = averageTime / iterations;
				trustedAverageTime1 = AvgTrustedInterval(averageTime, timeCalculations1, iterations);

				time1 = time;
				if (j == 0)
					calcData.Add(name1, 0, time, numbersAmount[nA]);
				else
					calcData.Add(name1, nT, time, numbersAmount[nA]);
				//Почему там аргументы превращаются в мусор?

				for (int t = 0; t < 3; t++)
				{
					averageTime = 0;
					Vector resultVector(numbersAmount[nA]);
					for (int iter = 0; iter < iterations; iter++)
					{
						resultVector = Vector(numbersAmount[nA]);

						switch (t)
						{
						case 0:
						{
							name2 = "Последоватьное суммирование";
							VectorsSumConsistently(vector1, vector2, resultVector, time);
							break;
						}
						case 1:
						{
							name2 = "Параллельное суммирование с for";
							VectorsSumParallelFor(vector1, vector2, resultVector, time);
							break;
						}
						default:
						{
							name2 = "Параллельное суммирование с sections";
							VectorsSumParallelSections(vector1, vector2, resultVector, time, omp_get_num_threads());
							break;
						}
						}

						timeCalculations2[iter] = time;
						averageTime += time;
					}
					averageTime = averageTime / iterations;
					trustedAverageTime2 = AvgTrustedInterval(averageTime, timeCalculations2, iterations);

					time2 = time;

					if (t == 0)
						calcData.Add(name2, 0, time, numbersAmount[nA]);
					else
						calcData.Add(name2, nT, time, numbersAmount[nA]);

					for (int i = 0; i < 3; i++)
					{
						averageTime = 0;
						for (int iter = 0; iter < iterations; iter++)
						{
							double sum = 0;

							switch (i)
							{
							case 0:
							{
								name3 = "Последоватьный подсчет суммы";
								MagnitudeVectorConsistently(resultVector, sum, time);
								break;
							}
							case 1:
							{
								name3 = "Параллельный (с редукторами) подсчет суммы";
								MagnitudeVectorReductor(resultVector, sum, time);
								break;
							}
							default:
							{
								name3 = "Параллельный (с крит. секциями) подсчет суммы";
								MagnitudeVectorCritical(resultVector, sum, time);
								break;
							}
							}

							timeCalculations3[iter] = time;
							averageTime += time;
						}
						averageTime = averageTime / iterations;
						trustedAverageTime3 = AvgTrustedInterval(averageTime, timeCalculations3, iterations);

						if (i == 0)
							calcData.Add(name3, 0, time, numbersAmount[nA]);
						else
							calcData.Add(name3, nT, time, numbersAmount[nA]);

						if (nA == 1) PrintResult(name1, name2, name3, trustedAverageTime1 + trustedAverageTime2 + trustedAverageTime3, nA, nT, fileName);
					}
				}
			}
		}
	}
	calcData.PrintResults(numbersAmount, "resultTable.csv");
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

	cout << endl << endl;
	string wait;
	cin >> wait;
}
