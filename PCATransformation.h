#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "gdal_priv.h"
#include "cpl_conv.h"
#include "gdal_alg.h"
#include "ogrsf_frmts.h"

//�Գ���ң��Ӱ�����PCA�任
class CPCATransformation
{
public:
	CPCATransformation();
	~CPCATransformation();
	//PCA�任ִ�к���
	string Execute(string sImageFile, string sOutputFile);

private:
	//����Ӱ���ֵ
	double* CalculateMean();
	//����Ӱ��Э�������
	double* CalculateCovMatrix(double *pMean);
	//�����Ÿ�ȣ�Jacobi��������ʵ�Գƾ����ȫ������ֵ����������
	bool Eejcb(double *a, double *v, double eps, int nMaxIterNum);
	//������ֵ��С������������
	void SortEigenVector(double *pEigenVector, double *pCovAfterEejcb);
	//��������������Ӱ�����PCA�任
	string GetPCAResult(string sOutputFile, double *pEigenVector);
	//��Ӱ��
	string OpenDataset(string sImageFile);
	//�ر�Ӱ��
	void CloseDataset();

	//���Ժ������������Э������������������
	void OutputMatrix(double *pMatrix, string sOutputFile);

	void OutputVector(double *pVector, string sOutputFile);


private:
	GDALDataset *m_pDataset;                //Ӱ�����ݼ�
	int m_nRasterXSize;                     //Ӱ����
	int m_nRasterYSize;                     //Ӱ��߶�
	int m_nRasterCount;                     //Ӱ�񲨶���
};

