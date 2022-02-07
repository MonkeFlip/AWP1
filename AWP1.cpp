#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
#include <immintrin.h>
#include <cmath>

bool Compare(float**** m1, float**** m2, float**** m3, int x);

void VectorizedAdd(int x_matrix_size, int y_matrix_size, float** matrixToAdd, float** result);
void NonVectorizedAdd(int x_matrix_size, int y_matrix_size, float** matrixToAdd, float** result);

void MultiplicationWithVectorization(int y_matrix1_size, int x_matrix2_size, int x_matrix1_size, float**  A, float**  B, float**  C);
void MultiplicationWithoutVectorization(int y_matrix1_size, int x_matrix2_size, int x_matrix1_size, float** A, float** B, float** C);

void DisplayMatrix(float**** matrix, int x_matrix_size, int y_matrix_size);
void FillMatrix(float**** matrix, int x_matrix_size, int y_matrix_size);

void ManuallyVectorizedMultiplication(int y_matrix1_size, int x_matrix2_size, int x_matrix1_size, float**  A, float**  B, float**  C);

void ClearMatrix(float**** matrix, int x_matrix_size, int y_matrix_size);

const int size = 100;

int x_matrix1 = 10;
int y_matrix1 = 40;
int x_matrix2 = y_matrix1;
int y_matrix2 = x_matrix1;
int x_result = y_matrix1;
int y_result = x_matrix2;

int main()
{
	using namespace std::chrono;
	float**** matrix1;
	float**** matrix2;
	float**** result3;
	float**** result1;
	float**** result2;

	matrix1 = new float***[size];
	matrix2 = new float*** [size];
	result3 = new float*** [size];
	result1 = new float*** [size];
	result2 = new float*** [size];
	float** temp;

	for (int i = 0; i < size; i++)
	{
		matrix1[i] = new float** [size];
		matrix2[i] = new float** [size];
		result3[i] = new float** [size];
		result1[i] = new float** [size];
		result2[i] = new float** [size];
	}
	

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			matrix1[i][j] = new float* [y_matrix1];
			matrix2[i][j] = new float* [y_matrix2];
			result3[i][j] = new float* [y_result];
			result1[i][j] = new float* [y_result];
			result2[i][j] = new float* [y_result];
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
				result3[i][j][k] = new float[x_result];
				result1[i][j][k] = new float[x_result];
				result2[i][j][k] = new float[x_result];
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
	ClearMatrix(result1, x_matrix2, x_matrix2);
	ClearMatrix(result2, x_matrix2, x_matrix2);
	ClearMatrix(result3, x_matrix2, x_matrix2);
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int r = 0; r < size; r++)
			{
				MultiplicationWithVectorization(y_matrix1, x_matrix2, x_matrix1, matrix1[i][r], matrix2[r][j], temp);
				VectorizedAdd(y_matrix1, x_matrix2, temp, result1[i][j]);
			}
		}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "Vectorized method time: " << time_span.count() << " seconds." << std::endl;
	//DisplayMatrix(result1, x_result, y_result);

	//ClearMatrix(result1, x_result, y_result);

	t1 = high_resolution_clock::now();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int r = 0; r < size; r++)
			{
				MultiplicationWithoutVectorization(y_matrix1, x_matrix2, x_matrix1, matrix1[i][r], matrix2[r][j], temp);
				NonVectorizedAdd(y_matrix1, x_matrix2, temp, result2[i][j]);
			}
		}
	}
	t2 = high_resolution_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "Not vectorized method time: " << time_span.count() << " seconds." << std::endl;
	//DisplayMatrix(result2, x_result, y_result);
	//ClearMatrix(result, x_result, y_result);

	t1 = high_resolution_clock::now();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int r = 0; r < size; r++)
			{
				ManuallyVectorizedMultiplication(y_matrix1, x_matrix2, x_matrix1, matrix1[i][r], matrix2[r][j], temp);
				NonVectorizedAdd(y_matrix1, x_matrix2, temp, result3[i][j]);
			}
		}
	}
	t2 = high_resolution_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "Manually vectorized method time: " << time_span.count() << " seconds." << std::endl;
	//DisplayMatrix(result3, x_matrix2, x_matrix2);
	if (Compare(result1, result2, result3, x_matrix2))
	{
		std::cout << "Matrices are equal." << std::endl;;
	}
	else
	{
		std::cout << "Matrices aren't equal." << std::endl;;
	}
}

void MultiplicationWithVectorization(int y_matrix1_size, int x_matrix2_size, int x_matrix1_size, float**  A, float**  B, float**  C)
{
	for (int i = 0; i < y_matrix1_size; ++i)
	{
		float* c = C[i];
		for (int j = 0; j < x_matrix2_size; ++j)
			c[j] = 0;
		for (int k = 0; k < x_matrix1_size; ++k)
		{
			const float* b = B[k];
			float a = A[i][k];
			for (int j = 0; j < x_matrix2_size; ++j)
				c[j] += a * b[j];
		}
	}
}

void MultiplicationWithoutVectorization(int y_matrix1_size, int x_matrix2_size, int x_matrix1_size, float** A, float** B, float** C)
{
	for (int i = 0; i < y_matrix1_size; ++i)
	{
		float* c = C[i];
		for (int j = 0; j < x_matrix2_size; ++j)
			c[j] = 0;
		for (int k = 0; k < x_matrix1_size; ++k)
		{
			const float* b = B[k];
			float a = A[i][k];
#pragma loop(no_vector)
			for (int j = 0; j < x_matrix2_size; ++j)
				c[j] += a * b[j];
		}
	}
}

void ManuallyVectorizedMultiplication(int y_matrix1_size, int x_matrix2_size, int x_matrix1_size, float**  A, float**  B, float**  C)
{
	for (int i = 0; i < y_matrix1_size; ++i)
	{
		float* c = C[i];
		for (int j = 0; j < x_matrix2_size; j += 8)
			_mm256_storeu_ps(c + j + 0, _mm256_setzero_ps());
		for (int k = 0; k < x_matrix1_size; ++k)
		{
			const float* b = B[k];
			__m256 a = _mm256_set1_ps(A[i][k]);
			for (int j = 0; j < x_matrix2_size; j += 8)
			{
				/*_mm256_store_ps(c + j + 0, _mm256_fmadd_ps(a,
					_mm256_load_ps(b + j + 0), _mm256_load_ps(c + j + 0)));*/
				//////////////////////////
				__m256 temp = _mm256_mul_ps(a, _mm256_load_ps(b + j + 0));
				_mm256_store_ps(c + j + 0, _mm256_add_ps(temp, _mm256_load_ps(c + j + 0)));
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

void DisplayMatrix(float**** matrix, int x_matrix_size, int y_matrix)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int r = 0; r < y_matrix; r++)
			{
				for (int k = 0; k < x_matrix_size; k++)
				{
					std::cout << matrix[i][j][r][k]<<" ";
				}
			}
		}
	}
	std::cout << std::endl;
}

bool Compare(float**** m1, float**** m2, float**** m3, int x)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int k = 0; k < x; k++)
			{
				for (int r = 0; r < x; r++)
				{
					float epsilon = 1;
					/*if (abs(m1[i][j][k][r] - m2[i][j][k][r]) >= epsilon || abs(m2[i][j][k][r] - m3[i][j][k][r]) >= epsilon)
					{
						std::cout << "m1: " << m1[i][j][k][r] << std::endl;
						std::cout << "m2: " << m2[i][j][k][r] << std::endl;
						std::cout << "m3: " << m3[i][j][k][r] - m2[i][j][k][r] << std::endl;
						return false;
					}*/
					if ((m1[i][j][k][r] != m2[i][j][k][r]) || (m2[i][j][k][r] != m3[i][j][k][r]))
					{
						std::cout << "i: "<< i << " j: " << j << " k: " <<k << " r: " << r << std::endl;
						std::cout << "m1: " << m1[i][j][k][r] - m2[i][j][k][r] << std::endl;
						std::cout << "m2: " << m2[i][j][k][r] << std::endl;
						std::cout << "m3: " << m3[i][j][k][r] - m2[i][j][k][r] << std::endl;
						return false;
					}
				}
			}
		}
	}

	return true;
}

void FillMatrix(float**** matrix, int x_matrix_size, int y_matrix_size)
{
	float counter = 0;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int r = 0; r < y_matrix_size; r++)
			{
				for (int k = 0; k < x_matrix_size; k++)
				{
					matrix[i][j][r][k] = counter;
					counter += 0.1;
				}
			}
		}
	}
}

void ClearMatrix(float**** matrix, int x_matrix_size, int y_matrix_size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int r = 0; r < y_matrix_size; r++)
			{
				for (int k = 0; k < x_matrix_size; k++)
				{
					matrix[i][j][r][k] = 0;
				}
			}
		}
	}
}