#include "Matrix.hpp"
#include <iostream>

namespace cie
{

// default constructor that creates an empty matrix
DoubleMatrix::DoubleMatrix() : 
	length_(0), rows_(0), columns_(0), matrix_(new double[0])
{}

// constructor that creates a matrix filled with zeroes
DoubleMatrix::DoubleMatrix(size_t rows, size_t columns):
	length_(rows * columns), rows_(rows), columns_(columns), matrix_(new double[rows * columns])
{
	// fills matrix with zeros
	for (size_t i = 0; i < length_; ++i)
	{
		matrix_[i] = 0.0;
	}
}

// constructor that creates a square identity matrix
DoubleMatrix::DoubleMatrix(size_t size) :
	length_(size * size), rows_(size), columns_(size), matrix_(new double[size * size])
{
	// fills matrix with zeros
	for (size_t i = 0; i < length_; ++i)
	{
		matrix_[i] = 0.0;
	}
	// fills diagonals with 1
	for (size_t i = 0; i < length_; i = i + size + 1)
	{
		matrix_[i] = 1.0;
	}
}

// copy constructor
DoubleMatrix::DoubleMatrix(DoubleMatrix & m):
	length_(m.numberOfRows() * m.numberOfColumns()), rows_(m.numberOfRows()), columns_(m.numberOfColumns()),
	matrix_(new double[m.numberOfRows() * m.numberOfColumns()])
{
	for (size_t i = 0; i < length_; ++i)
	{
		matrix_[i] = m.matrix_[i];
	}
}

// destructor
DoubleMatrix::~DoubleMatrix()
{
	delete[] matrix_;
}

// linear to 2D transformation
double & DoubleMatrix::at(size_t row, size_t columns)
{
	return matrix_[row * columns_ + columns];
}

// public function to give the number of rows
size_t DoubleMatrix::numberOfRows() const
{
	return rows_;
}

// public function to give the number of columns
size_t DoubleMatrix::numberOfColumns() const
{
	return columns_;
}

// resize function - works for all combinations
void DoubleMatrix::resize(size_t newRows, size_t newColumns)
{
	// assigns new length of the matrix
	size_t newLength_ = newRows * newColumns;
	// creates a temporary matrix 
	double* temp_matrix = new double[newLength_];

	// different conditions of size change
	if (newRows >= rows_ && newColumns >= columns_) // both rows and columns increase
	{
		// fill temporary matrix with zeroes
		for (size_t i = 0; i < newLength_; ++i)
		{
			temp_matrix[i] = 0;
		}
		// overwrites temporary matrix with old data from matrix_
		for (size_t i = 0; i < columns_; ++i)
		{
			for (size_t j = 0; j < rows_; ++j)
			{
				temp_matrix[j * newColumns + i] = at(j, i);
			}
		}
	}
	else if(newRows < rows_ && newColumns < columns_) // both rows and columns decrease
	{
		// overwrites temporary matrix with old data from matrix_
		for (size_t i = 0; i < newColumns; ++i)
		{
			for (size_t j = 0; j < newRows; ++j)
			{
				temp_matrix[j * newColumns + i] = at(j, i);
			}
		}
	}
	else if (newRows < rows_ && newColumns >= columns_) // rows decrease, columns increase
	{
		// fill temporary matrix with zeroes
		for (size_t i = 0; i < newLength_; ++i)
		{
			temp_matrix[i] = 0;
		}
		// overwrites temporary matrix with old data from matrix_
		for (size_t i = 0; i < columns_; ++i)
		{
			for (size_t j = 0; j < newRows; ++j)
			{
				temp_matrix[j * newColumns + i] = at(j, i);
			}
		}
	}
	else // rows increase, columns decrease
	{
		// fill temporary matrix with zeroes
		for (size_t i = 0; i < newLength_; ++i)
		{
			temp_matrix[i] = 0;
		}
		// overwrites temporary matrix with old data from matrix_
		for (size_t i = 0; i < newColumns; ++i)
		{
			for (size_t j = 0; j < rows_; ++j)
			{
				temp_matrix[j * newColumns + i] = at(j, i);
			}
		}
	}

	// deletes original matrix_
	delete[] matrix_;
	// creates new matrix_ with new size
	matrix_ = new double[newLength_];
	// assignes temp_matrix to new matrix_
	for (size_t i = 0; i < newLength_; ++i)
	{
		matrix_[i] = temp_matrix[i];
	}
	
	// deletes temp_matrix
	delete[] temp_matrix;
	// updates length, columns and rows
	length_ = newLength_;
	columns_ = newColumns;
	rows_ = newRows;
}

// push new line to the right
void DoubleMatrix::push_right(double value)
{
	resize(rows_, columns_ + 1);
	size_t i = columns_ - 1;
	for (size_t j = 0; j < rows_; ++j)
	{
		matrix_[j * columns_ + i] = value;
	}
}

// pust new line to the bottom
void DoubleMatrix::push_down(double value)
{
	resize(rows_ + 1, columns_);
	size_t j = rows_ - 1;
	for (size_t i = 0; i < columns_; ++i)
	{
		matrix_[j * columns_ + i] = value;
	}
}

// push new line to the left
void DoubleMatrix::push_left(double value)
{
	// length icreases by one column which is n-rows
	size_t newColumns_ = columns_ + 1;
	size_t newLength_ = newColumns_ * rows_;
	// creates a temporary matrix 
	double* temp_matrix = new double[newLength_];

	// fill temporary matrix with values
	for (size_t j = 0; j < rows_; ++j)
	{
		temp_matrix[j * newColumns_] = value;

	}
	// overwrites temporary matrix with old data from matrix_
	for (size_t i = 0; i < columns_; ++i)
	{
		for (size_t j = 0; j < rows_; ++j)
		{
			temp_matrix[j * newColumns_ + i + 1] = at(j, i);
		}
	}

	// deletes original matrix_
	delete[] matrix_;
	// creates new matrix_ with new size
	matrix_ = new double[newLength_];
	// assignes temp_matrix to new matrix_
	for (size_t i = 0; i < newLength_; ++i)
	{
		matrix_[i] = temp_matrix[i];
	}

	// deletes temp_matrix
	delete[] temp_matrix;
	// updates length, columns and rows
	length_ = newLength_;
	columns_ = newColumns_;
}

// push new line to the top
void DoubleMatrix::push_up(double value)
{
	// length icreases by one column which is n-rows
	size_t newRows_ = rows_ + 1;
	size_t newLength_ = newRows_ * columns_;
	// creates a temporary matrix 
	double* temp_matrix = new double[newLength_];

	// fill temporary matrix with values
	for (size_t i = 0; i < columns_; ++i)
	{
		temp_matrix[i] = value;
	}
	// overwrites temporary matrix with old data from matrix_
	for (size_t i = 0; i < columns_; ++i)
	{
		for (size_t j = 0; j < rows_; ++j)
		{
			temp_matrix[(j + 1) * columns_ + i] = at(j, i);
		}
	}

	// deletes original matrix_
	delete[] matrix_;
	// creates new matrix_ with new size
	matrix_ = new double[newLength_];
	// assignes temp_matrix to new matrix_
	for (size_t i = 0; i < newLength_; ++i)
	{
		matrix_[i] = temp_matrix[i];
	}

	// deletes temp_matrix
	delete[] temp_matrix;
	// updates length, columns and rows
	length_ = newLength_;
	rows_ = newRows_;
}

// assignment operator
DoubleMatrix & DoubleMatrix::operator=(DoubleMatrix & other)
{
	// assigns new number of rows and columns to the matrix
	rows_ = other.numberOfRows();
	columns_ = other.numberOfColumns();
	length_ = rows_ * columns_;

	// deletes old matrix and creates new one with new size
	delete[] matrix_;
	matrix_ = new double[length_];

	// copies the data to the matrix
	for (size_t i = 0; i < length_; ++i)
	{
		matrix_[i] = other.matrix_[i];
	}
	return *this; // I don't know why we have to put this here
}

// addition operator
DoubleMatrix & DoubleMatrix::operator+(DoubleMatrix & other)
{
	// checks if two matrices are the same size
	if (rows_ == other.numberOfRows() && columns_ == other.numberOfColumns())
	{
		// adds one matrix to another
		for (size_t i = 0; i < length_; ++i)
		{
			matrix_[i] += other.matrix_[i];
		}
	}
	else
	{
		std::cout << "\nMatrices are not right dimensions, addition not possible!" << std::endl;
	}
	return *this; // I don't know why we have to put this here
}

DoubleMatrix & DoubleMatrix::operator-(DoubleMatrix & other)
{
	// checks if two matrices are the same size
	if (rows_ == other.numberOfRows() && columns_ == other.numberOfColumns())
	{
		// substracts one matrix from another
		for (size_t i = 0; i < length_; ++i)
		{
			matrix_[i] -= other.matrix_[i];
		}
	}
	else
	{
		std::cout << "\nMatrices are not right dimensions, addition not possible!" << std::endl;
	}
	return *this; // I don't know why we have to put this here
}

// print function
const void print(DoubleMatrix& m)
{
	std::cout << "\nYour matrix is:" << std::endl;
	for (size_t i = 0; i < m.numberOfRows(); ++i)
	{
		for (size_t j = 0; j < m.numberOfColumns(); ++j)
		{
			std::cout << m.at(i, j) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

const void print_dots(DoubleMatrix & m)
{
	{
		std::cout << "\nYour matrix is:" << std::endl;
		for (size_t i = 0; i < m.numberOfRows(); ++i)
		{
			for (size_t j = 0; j < m.numberOfColumns(); ++j)
			{
				if (m.at(i, j) != 0.0)
				{
					std::cout << "X ";
				}
				else
				{
					std::cout << "- ";
				}
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}



// ask for user input with some fancy GI
void userInput(DoubleMatrix & m)
{
	// print out the grid of the matrix
	std::cout << "Please enter values!\n\n     ";	
	for (size_t i = 0; i < m.numberOfColumns(); ++i)
	{
		std::cout << i << " ";
	}
	std::cout << std::endl;
	std::cout << "    ";
	for (size_t i = 0; i < m.numberOfColumns()*2; ++i)
	{
		std::cout <<"_";
	}
	std::cout << std::endl;
	for (size_t i = 0; i < m.numberOfRows(); ++i)
	{
		std::cout << i << "  | ";
		for (size_t j = 0; j < m.numberOfColumns(); ++j)
		{
			std::cout << ". ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	// ask for values
	for (size_t i = 0; i < m.numberOfRows(); ++i)
	{
		for (size_t j = 0; j < m.numberOfColumns(); ++j)
		{
			size_t a;
			std::cout << "(" << i << "," << j << ") = ";
			std::cin >> a;
			m.at(i, j) = a;
		}
	}
}





VEKTOR::VEKTOR() : 
	x_(0.0), y_(0.0), length_(0.0)
{}

VEKTOR::VEKTOR(double x, double y) :
	x_(x), y_(y)
{
	length_ = sqrt(x * x + y * y);
}

double VEKTOR::x()
{
	return x_;
}

double VEKTOR::y()
{
	return y_;
}

double VEKTOR::length()
{
	return length_;
}

void VEKTOR::norm()
{
	x_ = x_ / length_;
	y_ = y_ / length_;

	length_ = 1.0;
}

VEKTOR & VEKTOR::operator=(VEKTOR & const other)
{
	// assigns new variables
	x_ = other.x_;
	y_ = other.y_;
	length_ = other.length_;

	return *this; // I don't know why we have to put this here
}

VEKTOR & VEKTOR::operator+( VEKTOR & const other)
{
	// assigns new variables
	VEKTOR result;
	result.x_ = x_ + other.x_;
	result.y_ = y_ + other.y_;
	
	result.length_ = sqrt(result.x_ * result.x_ + result.y_ * result.y_);

	return result;
}

VEKTOR & VEKTOR::operator-(VEKTOR & const other)
{
	// assigns new variables
	VEKTOR result;
	result.x_ = x_ - other.x_;
	result.y_ = y_ - other.y_;

	result.length_ = sqrt(result.x_ * result.x_ + result.y_ * result.y_);

	return result;
}


const void print(VEKTOR & v)
{
	std::cout << "(" << v.x() << ", " << v.y() << ")" << std::endl;
	return void();
}

}