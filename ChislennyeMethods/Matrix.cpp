#include"Matrix.h"


Matrix::Matrix(int rows, int columns) {
	this->rows = rows;
	this->columns = columns;
	this->elements = new vector<double>(rows * columns);

}
Matrix::~Matrix() {
	delete elements;
}
double Matrix::GetRows() const {
	return rows;
}
double Matrix::GetColumns() const {
	return columns;
}
double Matrix::GetElement(int i, int j) const {
	return (*elements)[i * columns + j];
}
double Matrix::TDet() {
	double det = 1;
	for (int i = 0; i < columns; i++) {
		det *= (*elements)[i * columns + i];
	}
	return det;
}

Matrix* Matrix::GetRow(int i) {
	Matrix* b = new Matrix(1, columns);
	for (int j = 0; j < columns; j++)
		b->SetElement(0, j, this->GetElement(i, j));
	return b;
}

Matrix* Matrix::GetColumn(int j) {
	Matrix* b = new Matrix(rows, 1);
	for (int i = 0; i < rows; i++)
		b->SetElement(i, 0, this->GetElement(i, j));
	return b;
}

double Matrix::DotProduct(Matrix* b) {
	return this->Multiply(b)->GetElement(0, 0);
}

double Matrix::Norm() {
	Matrix* B = Transpose(this);
	//B->Print();
	return sqrt(B->DotProduct(this));
}
Matrix::Matrix(Matrix* m) {
	rows = m->GetRows();
	columns = m->GetColumns();
	elements = new vector<double>(rows * columns);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			(*this->elements)[i * columns + j] = m->GetElement(i, j);
		}
	}
}
void Matrix::SetElement(int i, int j, double value) {
	(*elements)[i * columns + j] = value;
}

void Matrix::Print() {
	cout << "rows" << ' ' << rows << endl;
	cout << "columns" << ' ' << columns << endl;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			cout << (*elements)[i * columns + j] << ' ';
		}
		cout << endl;
	}
}

Matrix* Matrix::Gauss(Matrix* b) {
	Matrix* AForward = new Matrix(this);
	AForward->Print();
	Matrix* x = new Matrix(b);

	for (int k = 0; k < columns; k++) {
		for (int i = k+1; i < rows; i++) {
			double a = AForward->GetElement(i, k) / AForward->GetElement(k, k);
			for (int j =k; j < columns; j++) {
				if (AForward->GetElement(k, k) != 0 && i!=k) {
					AForward->SetElement(i, j, AForward->GetElement(i, j) - a * AForward->GetElement(k, j));

				}
			}
			x->SetElement(i, 0, x->GetElement(i, 0) - a * x->GetElement(k, 0));
		}
	}
	for (int k = rows - 1; k >= 0; k--) {
		double sum = 0;
		for (int j = k + 1; j < rows; j++) {
			sum += (AForward->GetElement(k, j) * x->GetElement(j, 0));
		}
		x->SetElement(k, 0, (x->GetElement(k, 0) - sum) / AForward->GetElement(k, k));
	}
	AForward->Print();
	cout << AForward->TDet() << endl;
	return x;
}

Matrix* Matrix::JordanGauss(Matrix* b) {
	Matrix* AForward = new Matrix(this);
	AForward->Print();
	Matrix* x = new Matrix(b);

	for (int k = 0; k < columns; k++) {
		for (int i = k + 1; i < rows; i++) {
			double a = AForward->GetElement(i, k) / AForward->GetElement(k, k);
			for (int j = k; j < columns; j++) {
				if (AForward->GetElement(k, k) != 0 && i != k) {
					AForward->SetElement(i, j, AForward->GetElement(i, j) - a * AForward->GetElement(k, j));

				}
			}
			x->SetElement(i, 0, x->GetElement(i, 0) - a * x->GetElement(k, 0));
		}
	}
	for (int k = columns - 1; k >=0; k--) {
		double sum = 0;
		for (int j = k-1; j >=0 ; j--) {
			x->SetElement(j, 0, x->GetElement(j, 0) - x->GetElement(k, 0) * AForward->GetElement(j, k) / AForward->GetElement(k, k));
			AForward->SetElement(j, k,0);
		}
		x->SetElement(k, 0, x->GetElement(k, 0) / AForward->GetElement(k, k));
	}
	AForward->Print();
	cout << AForward->TDet() << endl;
	return x;
}

double Sign(double a) {
	if (a > 0)
		return 1;
	if (a < 0)
		return -1;
	return 0;
}

double Abs(double x) {
	if (x < 0)
		x = -x;
	return x;
}

Matrix* Matrix::Transpose(Matrix* S) {
	Matrix* T = new Matrix(S->GetColumns(), S->GetRows());
	for (int i = 0; i < S->GetRows(); i++) {
		for (int j = 0; j < S->GetColumns(); j++) {
			//double s = S->GetElement(j, i);
			T->SetElement(j, i, S->GetElement(i,j));
		}
	}
	return T;
}

Matrix* Matrix::SquareRoots(Matrix* b) {
	Matrix* S = new Matrix(this);
	Matrix* d = new Matrix(this);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			d->SetElement(i, j, 0);
			S->SetElement(i, j, 0);
		}
	}
	for (int i = 0; i < rows; i++) {
		double sum = 0;
		for (int p = 0; p < i; p++) {
			double s = S->GetElement(p, i);
			sum += (s * s * d->GetElement(p, p));
		}
		d->SetElement(i, i, Sign(this->GetElement(i, i) - sum));
		S->SetElement(i, i, sqrt(Abs(this->GetElement(i, i) - sum)));
		for (int j = i+1; j < columns; j++) {
			double sum = 0;
			for (int p = 0; p < i; p++) {
				sum += S->GetElement(p, i) * d->GetElement(p, p) * S->GetElement(p, j);
			}
			S->SetElement(i, j, (this->GetElement(i, j) - sum) / (d->GetElement(i,i) * S->GetElement(i, i)));
		}
	}
	cout << "d:" << endl;
	d->Print();
	cout << "S:" << endl;
	S->Print();
	Matrix* STD = Transpose(S)->Multiply(d);
	cout << "STD" << endl;
	STD->Print();
	Matrix* y = new Matrix(b);
	for (int k = 0; k < rows; k++) {
		for (int j = 0; j < k; j++) {
			y->SetElement(k, 0, y->GetElement(k, 0) - STD->GetElement(k, j) * y->GetElement(j, 0));
		}
		y->SetElement(k, 0, y->GetElement(k, 0) / STD->GetElement(k, k));
	}
	cout << "y:" << endl;
	y->Print();
	cout << "STD * y" << endl;
	(STD->Multiply(y))->Print();
	Matrix* x = new Matrix(y);
	for (int k = rows - 1; k >= 0; k--) {
		x->SetElement(k, 0, y->GetElement(k, 0));
		for (int j = k + 1; j < rows; j++) {
			x->SetElement(k, 0, x->GetElement(k, 0) - S->GetElement(k, j) * x->GetElement(j, 0));
		}
		x->SetElement(k, 0, x->GetElement(k, 0) / S->GetElement(k, k));
	}
	return x;
}


Matrix* Matrix::LU(Matrix* b) {
	Matrix* L = new Matrix(this);
	Matrix* U = new Matrix(this);
	this->Print();
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			if (i > j)
				U->SetElement(i, j, 0);
			if (i < j)
				L->SetElement(i, j, 0);
			if (i == j)
				L->SetElement(i, j, 1);
		}
	}
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			double a = this->GetElement(i, j);
			if (i <= j) {
				U->SetElement(i, j, a);
				for (int k = 0; k < i; k++) {
					U->SetElement(i, j, U->GetElement(i, j) - L->GetElement(i, k) * U->GetElement(k, j));
				}
			}
			else {
				L->SetElement(i, j, a);
				for (int k = 0; k < j; k++) {
					L->SetElement(i, j, L->GetElement(i, j) - L->GetElement(i, k) * U->GetElement(k, j));
				}

				L->SetElement(i, j, L->GetElement(i, j)/ U->GetElement(j,j));
			}

		}
	}
	cout << "L:" << endl;
	L->Print();
	cout << "Det L:" << ' ' << L->TDet() << endl;
	cout << "U:" << endl;
	U->Print();
	cout << "Det U:" << ' ' << U->TDet() << endl;
	Matrix* y = new Matrix(b);
	for (int k = 0; k < rows; k++) {
		for (int j = 0; j < k; j++) {
			y->SetElement(k, 0, y->GetElement(k, 0) - L->GetElement(k, j) * y->GetElement(j, 0));
		}
	}

	Matrix* x = new Matrix(y);
	for (int k = rows - 1; k >= 0; k--) {
		x->SetElement(k, 0, y->GetElement(k, 0));
		for (int j = k+1; j < rows; j++) {
			x->SetElement(k, 0, x->GetElement(k, 0) - U->GetElement(k, j) * x->GetElement(j, 0));
		}
		x->SetElement(k, 0, x->GetElement(k, 0) / U->GetElement(k, k));
	}
	return x;
}
Matrix* Matrix::MultByNumber(double x) {
	Matrix* A = new Matrix(this);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			A->SetElement(i, j, this->GetElement(i, j) * x);
		}
	}
	return A;
}

Matrix* Matrix::Minus(Matrix* B) {
	Matrix* A = new Matrix(rows, columns);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			A->SetElement(i, j, this->GetElement(i, j) - B->GetElement(i, j));
		}
	}
	return A;
}

Matrix* Matrix::Plus(Matrix* B) {
	Matrix* A = new Matrix(rows, columns);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			A->SetElement(i, j, this->GetElement(i, j) + B->GetElement(i, j));
		}
	}
	return A;
}
Matrix* Matrix::operator-(Matrix* b) {
	Matrix* A = new Matrix(rows, columns);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			A->SetElement(i, j, this->GetElement(i, j) - b->GetElement(i,j));
		}
	}
	return A;
}
Matrix* Matrix::operator+(Matrix* b) {
	Matrix* A = new Matrix(rows, columns);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			A->SetElement(i, j, this->GetElement(i, j) + b->GetElement(i, j));
		}
	}
	return A;
}
Matrix* Matrix::operator*(Matrix* b) {
	return this->Multiply(b);
}
Matrix* Matrix::operator*(double b) {
	return this->MultByNumber(b);
}
Matrix* Matrix::Ortho(Matrix* b) {
	Matrix* A = new Matrix(rows + 1, columns + 1);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			A->SetElement(i, j, this->GetElement(i, j));
		}
	}
	for (int i = 0; i < rows; i++) {
		A->SetElement(i, columns, -b->GetElement(i, 0));
	}
	for (int j = 0; j < columns; j++) {
		A->SetElement(rows, j, 0);
	}
	A->SetElement(rows, columns, 1);

	vector<Matrix*> u(6);
	vector<Matrix*> v(6);
	u[0] = Transpose(A->GetRow(0));
	v[0] = u[0]->MultByNumber(1 / u[0]->Norm());
	//v[0]->Print();
	for (int i = 1; i <= rows; i++) {
		u[i] = Transpose(A->GetRow(i));
		//u[i]->Print();
		for (int j = 0; j < i; j++) {
			cout << A->GetRow(i)->DotProduct(v[j]) << endl;
			u[i]->Print();
			u[i]= u[i]->Minus(v[j]->MultByNumber (A->GetRow(i)->DotProduct(v[j])));
		}
		//cout << 1;
		
		v[i] = u[i]->MultByNumber(1 / u[i]->Norm());
		/*cout << "u:" << i << endl;
		u[i]->Print();
		cout << "v:" << i << endl;
		v[i]->Print();
		*/
		}
	Matrix* x = new Matrix(rows, 1);
	for (int i = 0; i < rows; i++) {
		x->SetElement(i, 0, v[rows]->GetElement(i, 0)/v[rows]->GetElement(rows,0));
	}
	return x;
}
Matrix* Matrix::Multiply(Matrix* b) {
	if (columns == b->GetRows()) {
		Matrix* result = new Matrix(rows, b->GetColumns());
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < b->GetColumns(); j++) {
				double sum = 0;
				for (int k = 0; k < columns; k++) {
					sum += (*elements)[i * columns + k] * b->GetElement(k, j);
				}
				result->SetElement(i, j, sum);
			}
		}
		return result;
	}
	return new Matrix(0, 0);
}