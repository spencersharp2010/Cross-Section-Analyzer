#ifndef MATRIX_HPP
#define MATRIX_HPP

namespace cie
{
	class DoubleMatrix
	{
	public:
		// default constructor - size 0
		DoubleMatrix();

		// default constructor - 0 initialized
		DoubleMatrix(size_t rows, size_t columns);

		// default unity square matrix constructor
		DoubleMatrix(size_t size);

		// copy constructor
		DoubleMatrix(DoubleMatrix& m);

		// destructor
		~DoubleMatrix();

		// linear to 2d transformation
		double& at(size_t row, size_t columns);

		// return matrix dimensions
		size_t numberOfRows() const;
		size_t numberOfColumns() const;
		
		// matrix resize function - work for any new size
		void resize(size_t newRows, size_t newColumns);

		// matrix push functions - easy
		void push_right(double value);
		void push_down(double value);
		// matrix push functions - hard
		void push_left(double value);
		void push_up(double value);

		// assignment operator
		DoubleMatrix& operator= (DoubleMatrix& other);
		// addition operator
		DoubleMatrix& operator+ (DoubleMatrix& other);
		// substraction operator
		DoubleMatrix& operator- (DoubleMatrix& other);

	private:
		size_t length_;
		size_t rows_;
		size_t columns_;
		double* matrix_;
	};
	// free floating functions
	const void print(DoubleMatrix& m);
	const void print_dots(DoubleMatrix& m);
	void userInput(DoubleMatrix& m);
	
	
	class VEKTOR
	{
	public:
		// default constructor - zero initialized
		VEKTOR();
		// constructor
		VEKTOR(double x, double y);


		// public access to private members
		double x();
		double y();
		double length();

		// normalization
		void norm();

		// assignment operator
		VEKTOR& operator= (VEKTOR& const other);
		// addition operator
		VEKTOR& operator+ (VEKTOR& const other);
		// substraction operator
		VEKTOR& operator- (VEKTOR& const other);

	private:
		double x_;
		double y_;
		double length_;
	};

	const void print(VEKTOR& v);

} //namespace cie

#endif // !MATRIX_HPP
