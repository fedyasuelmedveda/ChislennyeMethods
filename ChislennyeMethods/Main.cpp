#include"Matrix.h"

int main() {
	double NN = 10;
	double j = 1.5 + 0.1 * NN;
	double k = NN;
	double l = k;
	double m = 1;

	Matrix* A = new Matrix(5, 5);

	A->SetElement(0, 0, j * m);
	A->SetElement(0, 1, j * 0.5);
	A->SetElement(0, 2, 0);
	A->SetElement(0, 3, 0.2 * l);
	A->SetElement(0, 4, 0);

	A->SetElement(1, 0, j * 0.5);
	A->SetElement(1, 1, j);
	A->SetElement(1, 2, 0.3 * j);
	A->SetElement(1, 3, 0);
	A->SetElement(1, 4, 0.1 * l);

	A->SetElement(2, 0, 0);
	A->SetElement(2, 1, j * 0.3);
	A->SetElement(2, 2, 10);
	A->SetElement(2, 3, -0.3 * j);
	A->SetElement(2, 4, 0.5 * l);

	A->SetElement(3, 0, 0.2 * k);
	A->SetElement(3, 1, 0);
	A->SetElement(3, 2, -0.3 * j);
	A->SetElement(3, 3, j);
	A->SetElement(3, 4, -0.1 * j);

	A->SetElement(4, 0, 0);
	A->SetElement(4, 1, 0.1 * k);
	A->SetElement(4, 2, 0.5 * k);
	A->SetElement(4, 3, -0.1 * j);
	A->SetElement(4, 4, j * m);
	Matrix* b = new Matrix(5, 1);
	b->SetElement(0, 0, -j + 0.05 * j * j);
	b->SetElement(1, 0, -0.8 * j + 0.1 * j * j - 0.02 * l * j);
	b->SetElement(2, 0, -10 + 0.03 * j * j - 0.1 * l * j);
	b->SetElement(3, 0, -0.2 * k + 0.3 * j + 0.02 * j * j);
	b->SetElement(4, 0, 0.01 * j * k - 0.5 * k - 0.2 * j * j);
	
	A->Print();
	b->Print();
	
	Matrix* x = A->Ortho(b);
	cout << 'x' << endl;
	x->Print();
	
	Matrix* C = A->Multiply(x);
	cout << 'C' << endl;

	C->Print();
	cout << 'b' << endl;

	b->Print();
	/*
	Matrix* D = Matrix::Transpose(A);
	Matrix* r2 = A->GetRow(2);
	Matrix* c3 = A->GetColumn(3);
	cout << "Dot product:" << endl;
	cout << r2->DotProduct(c3) << endl;
	cout << "Norm c3:" << endl << c3->Norm() << endl;
	//D->Print();
	*/
	cout << "det:" << endl << A->DetRecursive(1) << endl;
}