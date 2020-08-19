#include "PCATransformation.h"



CPCATransformation::CPCATransformation()
{
	GDALAllRegister();
	OGRRegisterAll();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	CPLSetConfigOption("SHAPE_ENCODING", "");
	m_pDataset = NULL;
	m_nRasterXSize = 0;
	m_nRasterYSize = 0;
	m_nRasterCount = 0;
}


CPCATransformation::~CPCATransformation()
{
	CloseDataset();
}

//执行PCA变换
string CPCATransformation::Execute(string sImageFile, string sOutputFile)
{
	string sResult = "";
	//打开数据集
	sResult = OpenDataset(sImageFile);
	if (sResult != "")
		return sResult;

	//计算影像均值
	double *pMean = CalculateMean();
	OutputVector(pMean, "D:\\Mean.txt");
	//计算影像协方差矩阵
	double *pCovMatrix = CalculateCovMatrix(pMean);

	OutputMatrix(pCovMatrix, "D:\\Matrix1.txt");

	double eps = 0.0001;  //控制精度要求
	double *pEigenVector = new double[m_nRasterCount * m_nRasterCount];

	//求解特征值及特征向量
	bool b = Eejcb(pCovMatrix, pEigenVector, eps, 100000);
	OutputMatrix(pCovMatrix, "D:\\Matrix2.txt");
	OutputMatrix(pCovMatrix, "D:\\Vector1.txt");
	if (b == false)
	{
		CloseDataset();
		delete[]pMean;
		delete[]pCovMatrix;
		delete[]pEigenVector;
		return "Calculate Eigen Vector Failed!";
	}

	//按特征值大小对特征向量排序
	SortEigenVector(pEigenVector, pCovMatrix);

	OutputMatrix(pCovMatrix, "D:\\Matrix3.txt");
	OutputMatrix(pCovMatrix, "D:\\Vector2.txt");

	//求解PCA变换结果
	sResult = GetPCAResult(sOutputFile, pEigenVector);
	if(sResult != "")
	{
		CloseDataset();
		delete[]pMean;
		delete[]pCovMatrix;
		delete[]pEigenVector;
		return sResult;
	}

	CloseDataset();
	delete[]pMean;
	delete[]pCovMatrix;
	delete[]pEigenVector;
	return "";
}

void CPCATransformation::CloseDataset()
{
	if (m_pDataset)
	{
		GDALClose(m_pDataset);
		m_pDataset = NULL;
	}
}

string CPCATransformation::OpenDataset(string sImageFile)
{
	m_pDataset = (GDALDataset*)GDALOpen(sImageFile.c_str(), GA_ReadOnly);
	if (m_pDataset == NULL)
	{
		return "Open File " + sImageFile + " Failed!";
	}
	m_nRasterXSize = m_pDataset->GetRasterXSize();
	m_nRasterYSize = m_pDataset->GetRasterYSize();
	m_nRasterCount = m_pDataset->GetRasterCount();
	return "";
}


//计算均值（考虑NoData）
double* CPCATransformation::CalculateMean()
{
	cout << "Calculate Mean Begin!" << endl;
	double *pMean = new double[m_nRasterCount];
	double *pNoDataValue = new double[m_nRasterCount];
	for(int i = 0; i < m_nRasterCount; i++)
	{
		pMean[i] = 0;
		pNoDataValue[i] = m_pDataset->GetRasterBand(i + 1)->GetNoDataValue();
	}

	double *pData = new double[m_nRasterXSize * m_nRasterCount];
	int nCount = 0;
	bool bIsNoData = false;

	int nOldPercent = 0, nNewPercent = 0;
	cout << "0%";
	for(int y = 0; y < m_nRasterYSize; y++)
	{
		m_pDataset->RasterIO(GF_Read, 0, y, m_nRasterXSize, 1, pData, m_nRasterXSize, 1,
			GDT_Float64, m_nRasterCount, NULL, 0, 0, 0);
		for(int x = 0; x < m_nRasterXSize; x++)
		{
			bIsNoData = false;
			for(int n = 0; n < m_nRasterCount; n++)
			{
				if(pData[n * m_nRasterXSize + x] == pNoDataValue[n])
				{
					bIsNoData = true;
					break;
				}
			}
			if(bIsNoData == false)
			{
				nCount++;
				for(int n = 0; n < m_nRasterCount; n++)
				{
					pMean[n] += pData[n * m_nRasterXSize + x];
				}
			}
		}

		nNewPercent = y * 100 / m_nRasterYSize;
		if(nOldPercent != nNewPercent)
		{
			nOldPercent = nNewPercent;
			cout << "\r" << nOldPercent << "%";
		}

	}
	cout << endl;

	delete[]pData;
	delete[]pNoDataValue;
	for(int n = 0; n < m_nRasterCount; n++)
	{
		pMean[n] /= nCount;
	}
	cout << "Calculate Mean Finished!" << endl;
	
	return pMean;
}

//计算协方差矩阵
double* CPCATransformation::CalculateCovMatrix(double *pMean)
{
	cout << "Calculate CovMatrix Begin!" << endl;
	double *pCovMatrix = new double[m_nRasterCount * m_nRasterCount];
	for (int i = 0; i < m_nRasterCount * m_nRasterCount; i++)
		pCovMatrix[i] = 0;

	double *pNoDataValue = new double[m_nRasterCount];
	for (int i = 0; i < m_nRasterCount; i++)
	{
		pNoDataValue[i] = m_pDataset->GetRasterBand(i + 1)->GetNoDataValue();
	}

	double *pData = new double[m_nRasterXSize * m_nRasterCount];
	int nCount = 0;
	bool bIsNoData = false;

	int nOldPercent = 0, nNewPercent = 0;
	cout << "0%";
	for (int y = 0; y < m_nRasterYSize; y++)
	{
		m_pDataset->RasterIO(GF_Read, 0, y, m_nRasterXSize, 1, pData, m_nRasterXSize, 1,
			GDT_Float64, m_nRasterCount, NULL, 0, 0, 0);
		for (int x = 0; x < m_nRasterXSize; x++)
		{
			bIsNoData = false;
			for (int n = 0; n < m_nRasterCount; n++)
			{
				if (pData[n * m_nRasterXSize + x] == pNoDataValue[n])
				{
					bIsNoData = true;
					break;
				}
			}
			if (bIsNoData == false)
			{
				nCount++;
				for(int i = 0; i < m_nRasterCount; i++)
				{
					for(int j = 0; j < m_nRasterCount; j++)
					{
						pCovMatrix[i * m_nRasterCount + j] += ((pData[i * m_nRasterXSize + x] - pMean[i]) * 
							(pData[j * m_nRasterXSize + x] - pMean[j]));
					}
				}
			}
		}

		nNewPercent = y * 100 / m_nRasterYSize;
		if (nOldPercent != nNewPercent)
		{
			nOldPercent = nNewPercent;
			cout << "\r" << nOldPercent << "%";
		}

	}
	cout << endl;

	delete[]pData;
	delete[]pNoDataValue;

	for (int i = 0; i < m_nRasterCount * m_nRasterCount; i++)
		pCovMatrix[i] /= (nCount - 1);
	cout << "Calculate CovMatrix Finished!" << endl;
	return pCovMatrix;
}

bool CPCATransformation::Eejcb(double *a, double *v, double eps, int nMaxIterNum)
{
	cout << "Calculate Eigen Vector Begin!" << endl;
	int i, j, p, q, u, w, t, s, l;
	double fm, cn, sn, omega, x, y, d;

	int n = m_nRasterCount;

	l = 1;
	//初始化特征向量矩阵使其全为0
	for (i = 0; i <= n - 1; i++)
	{
		v[i * n + i] = 1.0;
		for (j = 0; j <= n - 1; j++)
		{
			if (i != j)
				v[i * n + j] = 0.0;
		}
	}
	cout << "\rThe 0-th Iteration";
	while (true) //循环
	{
		fm = 0.0;
		for (i = 0; i <= n - 1; i++)   // 出, 矩阵a( 特征值 ), 中除对角线外其他元素的最大绝对值
		{
			//这个最大值是位于a[p][q] ,等于fm
			for (j = 0; j <= n - 1; j++)
			{
				d = fabs(a[i * n + j]);

				if ((i != j) && (d > fm))
				{
					fm = d;
					p = i;
					q = j;
				}
			}
		}

		if (fm < eps)   //精度复合要求
		{
			cout << "Calculate Eigen Vector Finished!" << endl;
			return true; //正常返回
		}

		if (l > nMaxIterNum)     //迭代次数太多
			return false;//失败返回

		l++;       //   迭代计数器
		u = p * n + q;
		w = p * n + p;
		t = q * n + p;
		s = q * n + q;
		x = -a[u];
		y = (a[s] - a[w]) / 2.0;		//x y的求法不同
		omega = x / sqrt(x * x + y * y);	//sin2θ

											//tan2θ=x/y = -2.0*a[u]/(a[s]-a[w])
		if (y < 0.0)
			omega = -omega;

		sn = 1.0 + sqrt(1.0 - omega * omega);
		sn = omega / sqrt(2.0 * sn);		//sinθ
		cn = sqrt(1.0 - sn * sn);			//cosθ

		fm = a[w];   //   变换前的a[w]   a[p][p]
		a[w] = fm * cn * cn + a[s] * sn * sn + a[u] * omega;
		a[s] = fm * sn * sn + a[s] * cn * cn - a[u] * omega;
		a[u] = 0.0;
		a[t] = 0.0;

		//   以下是旋转矩阵,旋转了了p行,q行,p列,q列
		//   但是四个特殊点没有旋转(这四个点在上述语句中发生了变化)
		//   其他不在这些行和列的点也没变
		//   旋转矩阵,旋转p行和q行
		for (j = 0; j <= n - 1; j++)
		{
			if ((j != p) && (j != q))
			{
				u = p * n + j;
				w = q * n + j;
				fm = a[u];
				a[u] = a[w] * sn + fm * cn;
				a[w] = a[w] * cn - fm * sn;
			}
		}

		//旋转矩阵,旋转p列和q列
		for (i = 0; i <= n - 1; i++)
		{
			if ((i != p) && (i != q))
			{
				u = i * n + p;
				w = i * n + q;
				fm = a[u];
				a[u] = a[w] * sn + fm * cn;
				a[w] = a[w] * cn - fm * sn;
			}
		}

		//记录旋转矩阵特征向量
		for (i = 0; i <= n - 1; i++)
		{
			u = i * n + p;
			w = i * n + q;
			fm = v[u];
			v[u] = v[w] * sn + fm * cn;
			v[w] = v[w] * cn - fm * sn;
		}

		cout << "\rThe" << l << "-th Iteration";
	}
	cout << "Calculate Eigen Vector Finished!" << endl;
	return true;
}

void CPCATransformation::SortEigenVector(double *pEigenVector, double *pCovAfterEejcb)
{
	cout << "Sort Eigen Vector Begin!" << endl;
	for (int i = 0; i < m_nRasterCount - 1; i++)
	{
		for (int j = i + 1; j < m_nRasterCount; j++)
		{
			if (pCovAfterEejcb[j * m_nRasterCount + j] > pCovAfterEejcb[i * m_nRasterCount + i])
			{
				double temp = 0;
				temp = pCovAfterEejcb[j * m_nRasterCount + j];
				pCovAfterEejcb[j * m_nRasterCount + j] = pCovAfterEejcb[i * m_nRasterCount + i];
				pCovAfterEejcb[i * m_nRasterCount + i] = temp;
				for (int k = 0; k < m_nRasterCount; k++)
				{
					temp = 0;
					temp = pEigenVector[k * m_nRasterCount + j];
					pEigenVector[k * m_nRasterCount + j] = pEigenVector[k * m_nRasterCount + i];
					pEigenVector[k * m_nRasterCount + i] = temp;
				}
			}
		}
	}
	cout << "Sort Eigen Vector Finished!" << endl;
}

string CPCATransformation::GetPCAResult(string sOutputFile, double *pEigenVector)
{
	cout << "Calculate PAC Result Begin!" << endl;
	double **pVector = new double*[m_nRasterCount];
	for(int i = 0; i < m_nRasterCount; i++)
	{
		pVector[i] = new double[m_nRasterCount];
		for(int j = 0; j < m_nRasterCount; j++)
		{
			pVector[i][j] = pEigenVector[i + j * m_nRasterCount];
		}
	}

	GDALDriver *pDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	GDALDataset *pDataset = pDriver->Create(sOutputFile.c_str(), m_nRasterXSize, m_nRasterYSize, m_nRasterCount, GDT_Float32, NULL);
	if (pDataset == NULL)
		return "Create File " + sOutputFile + " Failed!";

	double *pGeoTrans = new double[6];
	m_pDataset->GetGeoTransform(pGeoTrans);
	pDataset->SetGeoTransform(pGeoTrans);
	delete[]pGeoTrans;
	pDataset->SetProjection(m_pDataset->GetProjectionRef());

	double *pNoDataValue = new double[m_nRasterCount];
	for (int i = 0; i < m_nRasterCount; i++)
	{
		pNoDataValue[i] = m_pDataset->GetRasterBand(i + 1)->GetNoDataValue();
		pDataset->GetRasterBand(i + 1)->SetNoDataValue(pNoDataValue[i]);
	}

	double *pData = new double[m_nRasterXSize * m_nRasterCount];
	float *pOutputData = new float[m_nRasterXSize * m_nRasterCount];
	bool bIsNoData = false;

	int nOldPercent = 0, nNewPercent = 0;
	cout << "0%";
	for (int y = 0; y < m_nRasterYSize; y++)
	{
		m_pDataset->RasterIO(GF_Read, 0, y, m_nRasterXSize, 1, pData, m_nRasterXSize, 1,
			GDT_Float64, m_nRasterCount, NULL, 0, 0, 0);
		for (int x = 0; x < m_nRasterXSize; x++)
		{
			bIsNoData = false;
			for (int n = 0; n < m_nRasterCount; n++)
			{
				if (pData[n * m_nRasterXSize + x] == pNoDataValue[n])
				{
					bIsNoData = true;
					break;
				}
			}
			if (bIsNoData == false)
			{
				for(int i = 0; i < m_nRasterCount; i++)
				{
					pOutputData[i * m_nRasterXSize + x] = 0;
					for(int k = 0; k < m_nRasterCount; k++)
					{
						pOutputData[i * m_nRasterXSize + x] += (float)(pVector[i][k] * pData[k * m_nRasterXSize + x]);
					}
				}
			}
			else
			{
				for (int n = 0; n < m_nRasterCount; n++)
					pOutputData[n * m_nRasterXSize + x] = (float)pNoDataValue[n];
			}
		}

		pDataset->RasterIO(GF_Write, 0, y, m_nRasterXSize, 1, pOutputData, m_nRasterXSize, 1,
			GDT_Float32, m_nRasterCount, NULL, 0, 0, 0);

		nNewPercent = y * 100 / m_nRasterYSize;
		if (nOldPercent != nNewPercent)
		{
			nOldPercent = nNewPercent;
			cout << "\r" << nOldPercent << "%";
		}

	}
	cout << endl;
	
	for (int i = 0; i < m_nRasterCount; i++)
		delete[]pVector[i];

	delete[]pVector;
	delete[]pNoDataValue;
	delete[]pData;
	delete[]pOutputData;
	GDALClose(pDataset);
	cout << "Calculate PAC Result Begin!" << endl;
	return "";
}

void CPCATransformation::OutputMatrix(double *pMatrix, string sOutputFile)
{
	ofstream out(sOutputFile.c_str());
	out << setiosflags(ios::fixed) << setprecision(5);
	for(int i = 0; i < m_nRasterCount; i++)
	{
		for(int j = 0; j < m_nRasterCount; j++)
		{
			out << pMatrix[i * m_nRasterCount + j] << " ";
		}
		out << endl;
	}
	out.close();
}

void CPCATransformation::OutputVector(double *pVector, string sOutputFile)
{
	ofstream out(sOutputFile.c_str());
	out << setiosflags(ios::fixed) << setprecision(5);
	for (int i = 0; i < m_nRasterCount; i++)
	{
		out << pVector[i] << endl;
	}
	out.close();
}