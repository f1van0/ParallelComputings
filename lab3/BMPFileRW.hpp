#pragma once
#include <Windows.h>
#include <iostream>
#include <fstream>

void BMPWrite(RGBQUAD **&, BITMAPFILEHEADER&, BITMAPINFOHEADER&, const char *);
void BMPRead(RGBQUAD **&, BITMAPFILEHEADER&, BITMAPINFOHEADER&, const char *);


void BMPRead(RGBQUAD** &rgb, BITMAPFILEHEADER &header, \
	BITMAPINFOHEADER &bmiHeader, const char* fin)
{
    // ��������� ����
	std::ifstream InFile(fin, std::ios::binary);
	// ��������� ��������� �����
	InFile.read((char*)(&header), sizeof(BITMAPFILEHEADER));
	// ��������� ������������ ����� �����������
	InFile.read((char*)(&bmiHeader), sizeof(BITMAPINFOHEADER));
    // �������� ������ ��� ������ RGB �������� ��������� RGBQUAD
    rgb = new RGBQUAD*[bmiHeader.biHeight];
    for (int i = 0; i < bmiHeader.biHeight; i++)
    {
        rgb[i] = new RGBQUAD[bmiHeader.biWidth];
    }
    // ��������� ������ ����������� � ������ �������� RGB 
    for (int i = 0; i < bmiHeader.biHeight; i++)
    {
        for (int j = 0; j < bmiHeader.biWidth; j++)
        {
			InFile.read((char*)(&rgb[i][j]), 3); // .rgbBlue .rgbGreen .rgbRed;
        }
    }
    // ��������� ����
	InFile.close();
	
}

void BMPWrite(RGBQUAD** &rgb, BITMAPFILEHEADER &header, \
	BITMAPINFOHEADER &bmiHeader, const char *fout)
{
    // ��������� ���� ��� ������ ����������� � ������ BMP
	std::ofstream OutFile(fout, std::ios::binary);
	//// ���������� ��������� �����
	OutFile.write((char*)(&header), sizeof(BITMAPFILEHEADER));
	//// ���������� ������������ ����� �����������
	OutFile.write((char*)(&bmiHeader), sizeof(BITMAPINFOHEADER));
    // ���������� ������ ����������� �� ������ �������� RGB � ���� 
    for (int i = 0; i < bmiHeader.biHeight; i++)
        for (int j = 0; j < bmiHeader.biWidth; j++)
        {
			OutFile.write((char*)&(rgb[i][j]), 3);// .rgbBlue .rgbGreen .rgbRed;
        }
    // ��������� ����
		OutFile.close();
}
