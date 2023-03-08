#pragma once
#include<iostream>
#include<vector>
#include<math.h>
using namespace std;
class Matrix {
private:
	int rows, columns;
	vector<double>* elements;
public:
	Matrix(int rows, int columns);
	Matrix(Matrix* m);
	~Matrix();
	double GetRows() const;
	double GetColumns() const;
	double GetElement(int i, int j) const;
	void SetElement(int i, int j, double value);
	void Print();
	Matrix* GetRow(int i);
	Matrix* GetColumn(int j);
	double DotProduct(Matrix* b);
	double Norm();
	Matrix* Ortho(Matrix* b);
	Matrix* Gauss(Matrix* b);
	Matrix* LU(Matrix* b);
	Matrix* JordanGauss(Matrix* b);
	Matrix* SquareRoots(Matrix* b);
	Matrix* Multiply(Matrix* b);
	Matrix* MultByNumber(double x);
	Matrix* Plus(Matrix* B);
	Matrix* Minus(Matrix* B);
	static Matrix* Transpose(Matrix* S);
	//double Sign(double a);
	double TDet();
	Matrix* operator*(Matrix* b);
	Matrix* operator*(double b);
	Matrix* operator-(Matrix* b);
	Matrix* operator+(Matrix* b);
};