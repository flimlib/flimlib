#pragma once

template <class T>
class ParamMatrix {
public:
	// the internal array
	T **const arr;

	// size of the matrix
	const int nrow;
	const int ncol;

	ParamMatrix(int row, int col);

	ParamMatrix(T **arr, int row, int col);

	~ParamMatrix();

	static T **new_mat(int row, int col);
};

template <class T>
ParamMatrix<T>::ParamMatrix(int row, int col) 
	: arr(new_mat(row, col)), nrow(row), ncol(col) { }

template <class T>
ParamMatrix<T>::ParamMatrix(T **arr, int row, int col) 
	: arr(arr), nrow(row), ncol(col) { }

template <class T>
ParamMatrix<T>::~ParamMatrix() {
	if (this->arr)
		delete[] this->arr;
}

/* copied from GCI_ecf_Matrix() */
template <class T>
T **ParamMatrix<T>::new_mat(int row, int col) {
	if (row == 0 || col == 0)
		return NULL;
	
	size_t row_size = (size_t) row * sizeof(T *);
	size_t data_size = (size_t) row * col * sizeof(T);
	void *raw = new T[row_size + data_size];
	T **head = (T **) raw;
	T *data = (T *) (head + row);

	if (NULL == raw) {
		return NULL;
	}

	for (int i = 0; i < row; ++i) {
		head[i] = data;
		data += col;
	}
	return head;
}
