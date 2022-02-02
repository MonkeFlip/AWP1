#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
#include <immintrin.h>

void VectorizedAdd(int x_matrix, int y_matrix, float** matrixToAdd, float** result);
void NonVectorizedAdd(int x_matrix, int y_matrix, float** matrixToAdd, float** result);

void MultiplicationWithVectorization(int M, int N, int K, float** __restrict A, float** __restrict B, float** __restrict C);
void MultiplicationWithoutVectorization(int M, int N, int K, float** A, float** B, float** C);

void DisplayMatrix(float**** matrix, int x_matrix, int y_matrix);
void FillMatrix(float**** matrix, int x_matrix, int y_matrix);

void ManualVectorization(int M, int N, int K, float** __restrict A, float** __restrict B, float** __restrict C);

void ClearMatrix(float**** matrix, int x_matrix, int y_matrix);

const int size = 200;

int x_matrix1 = 1;
int y_matrix1 = 4;
int x_matrix2 = y_matrix1;
int y_matrix2 = x_matrix1;
int x_result = y_matrix1;
int y_result = x_matrix2;

int main()
{
	using namespace std::chrono;
	float**** matrix1;
	float**** matrix2;
	float**** result;

	matrix1 = new float***[size];
	matrix2 = new float*** [size];
	result = new float*** [size];
	float** temp;

	for (int i = 0; i < size; i++)
	{
		matrix1[i] = new float** [size];
		matrix2[i] = new float** [size];
		result[i] = new float** [size];
	}
	

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			matrix1[i][j] = new float* [y_matrix1];
			matrix2[i][j] = new float* [y_matrix2];
			result[i][j] = new float* [y_result];
		}
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int k = 0; k < y_matrix1; k++)
			{
				matrix1[i][j][k] = new float [x_matrix1];
			}
		}
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int k = 0; k < y_matrix2; k++)
			{
				matrix2[i][j][k] = new float[x_matrix2];
			}
		}
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int k = 0; k < y_result; k++)
			{
				result[i][j][k] = new float[x_result];
			}
		}
	}

	temp = new float* [y_matrix1];
	for (int i = 0; i < y_matrix1; i++)
	{
		temp[i] = new float[y_matrix1];
	}

	//end of memory allocation 

	FillMatrix(matrix1, x_matrix1, y_matrix1);
	FillMatrix(matrix2, x_matrix2, y_matrix2);
	//std::cout << "Matrix1:" << std::endl;
	//DisplayMatrix(matrix1, x_matrix1, y_matrix1);
	//std::cout << "Matrix2:" << std::endl;
	//DisplayMatrix(matrix2, x_matrix2, y_matrix2);

	////////////////////////////////////

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int r = 0; r < size; r++)
			{
				MultiplicationWithVectorization(y_matrix1, x_matrix2, x_matrix1, matrix1[i][r], matrix2[r][j], temp);
				VectorizedAdd(y_matrix1, x_matrix2, temp, result[i][j]);
			}
		}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "Vectorized method time: " << time_span.count() << " seconds." << std::endl;
	//DisplayMatrix(result, x_result, y_result);
	ClearMatrix(result, x_result, y_result);

	t1 = high_resolution_clock::now();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int r = 0; r < size; r++)
			{
				MultiplicationWithoutVectorization(y_matrix1, x_matrix2, x_matrix1, matrix1[i][r], matrix2[r][j], temp);
				NonVectorizedAdd(y_matrix1, x_matrix2, temp, result[i][j]);
			}
		}
	}
	t2 = high_resolution_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "Not vectorized method time: " << time_span.count() << " seconds." << std::endl;
	//DisplayMatrix(result, x_result, y_result);
	ClearMatrix(result, x_result, y_result);

	t1 = high_resolution_clock::now();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int r = 0; r < size; r++)
			{
				ManualVectorization(y_matrix1, x_matrix2, x_matrix1, matrix1[i][r], matrix2[r][j], temp);
				NonVectorizedAdd(y_matrix1, x_matrix2, temp, result[i][j]);
			}
		}
	}
	t2 = high_resolution_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "Manually vectorized method time: " << time_span.count() << " seconds." << std::endl;
	//DisplayMatrix(result, x_result, y_result);
	ClearMatrix(result, x_result, y_result);
}

void MultiplicationWithVectorization(int M, int N, int K, float** __restrict A, float** __restrict B, float** __restrict C)
{
	for (int i = 0; i < M; ++i)
	{
		float* c = C[i];
		for (int j = 0; j < N; ++j)
			c[j] = 0;
		for (int k = 0; k < K; ++k)
		{
			const float* b = B[k];
			float a = A[i][k];
			for (int j = 0; j < N; ++j)
				c[j] += a * b[j];
		}
	}
}

void MultiplicationWithoutVectorization(int M, int N, int K, float** A, float** B, float** C)
{
	for (int i = 0; i < M; ++i)
	{
		float* c = C[i];
		for (int j = 0; j < N; ++j)
			c[j] = 0;
		for (int k = 0; k < K; ++k)
		{
			const float* b = B[k];
			float a = A[i][k];
#pragma loop(no_vector)
			for (int j = 0; j < N; ++j)
				c[j] += a * b[j];
		}
	}
}

void ManualVectorization(int M, int N, int K, float** __restrict A, float** __restrict B, float** __restrict C)
{
	for (int i = 0; i < M; ++i)
	{
		float* c = C[i];
		for (int j = 0; j < N; j += 8)
			_mm256_storeu_ps(c + j + 0, _mm256_setzero_ps());
		for (int k = 0; k < K; ++k)
		{
			const float* b = B[k];
			__m256 a = _mm256_set1_ps(A[i][k]);
			for (int j = 0; j < N; j += 8)
			{
				_mm256_storeu_ps(c + j + 0, _mm256_fmadd_ps(a,
					_mm256_loadu_ps(b + j + 0), _mm256_loadu_ps(c + j + 0)));
			}
		}
	}
}

void VectorizedAdd(int x_matrix, int y_matrix, float** matrixToAdd, float** result)
{
	for (int i = 0; i < y_matrix; i++)
	{
		for (int j = 0; j < x_matrix; j++)
		{
			result[i][j] += matrixToAdd[i][j];
		}
	}
}

void NonVectorizedAdd(int x_matrix, int y_matrix, float** matrixToAdd, float** result)
{
	for (int i = 0; i < y_matrix; i++)
	{
#pragma loop(no_vector)
		for (int j = 0; j < x_matrix; j++)
		{
			result[i][j] += matrixToAdd[i][j];
		}
	}
}

void DisplayMatrix(float**** matrix, int x_matrix, int y_matrix)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int r = 0; r < y_matrix; r++)
			{
				for (int k = 0; k < x_matrix; k++)
				{
					std::cout << matrix[i][j][r][k]<<" ";
				}
			}
		}
	}
	std::cout << std::endl;
}

void FillMatrix(float**** matrix, int x_matrix, int y_matrix)
{
	int counter = 0;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int r = 0; r < y_matrix; r++)
			{
				for (int k = 0; k < x_matrix; k++)
				{
					matrix[i][j][r][k] = counter++;
				}
			}
		}
	}
}

void ClearMatrix(float**** matrix, int x_matrix, int y_matrix)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int r = 0; r < y_matrix; r++)
			{
				for (int k = 0; k < x_matrix; k++)
				{
					matrix[i][j][r][k] = 0;
				}
			}
		}
	}
}