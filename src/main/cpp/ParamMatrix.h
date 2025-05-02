/*-
 * #%L
 * FLIMLib package for exponential curve fitting of fluorescence lifetime data.
 * %%
 * Copyright (C) 2010 - 2025 University of Oxford and Board of Regents of the
 * University of Wisconsin-Madison.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
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
