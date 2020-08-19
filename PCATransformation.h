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

//对超大遥感影像进行PCA变换
class CPCATransformation
{
public:
	CPCATransformation();
	~CPCATransformation();
	//PCA变换执行函数
	string Execute(string sImageFile, string sOutputFile);

private:
	//计算影像均值
	double* CalculateMean();
	//计算影像协方差矩阵
	double* CalculateCovMatrix(double *pMean);
	//利用雅格比（Jacobi）方法求实对称矩阵的全部特征值及特征向量
	bool Eejcb(double *a, double *v, double eps, int nMaxIterNum);
	//按特征值大小排列特征向量
	void SortEigenVector(double *pEigenVector, double *pCovAfterEejcb);
	//根据特征向量对影像进行PCA变换
	string GetPCAResult(string sOutputFile, double *pEigenVector);
	//打开影像
	string OpenDataset(string sImageFile);
	//关闭影像
	void CloseDataset();

	//测试函数，用于输出协方差矩阵或者特征向量
	void OutputMatrix(double *pMatrix, string sOutputFile);

	void OutputVector(double *pVector, string sOutputFile);


private:
	GDALDataset *m_pDataset;                //影像数据集
	int m_nRasterXSize;                     //影像宽度
	int m_nRasterYSize;                     //影像高度
	int m_nRasterCount;                     //影像波段数
};

