#include "Matrix.hpp"
#include <iostream>
#include <vector>
#include <array>
#include <math.h>
#include <algorithm>
#include <stack>
#include <cstdlib>


// search function // 0-> code for false, anything else-> (position + 1) of the data points of interest
// this function finds the point only if it's group number is 0 (zero)
int Search(std::vector<std::vector<int>> points, int i, int j)
{
	// loops over all the data in points of interest
	for (int n = 0; n < points[0].size(); ++n)
	{
		if (points[0][n] == i && points[1][n] == j && points[2][n] == 0)
		{
			return n + 1;
		}
	}
	return 0;
}

// recursive find gropus function
void FindGroups(std::vector<std::vector<int>>& points, int i, int j, int group_number)
{
	
	// finds the position of the point and changes it's group number
	int position = Search(points, i, j);
	if (position != 0)
	{
		points[2][position - 1] = group_number;
		
	}

	//go right
	if (Search(points, i, j + 1) != 0)
	{
		FindGroups(points, i, j + 1, group_number);
	}
	//go down
	if (Search(points, i + 1, j) != 0)
	{
		FindGroups(points, i + 1, j, group_number);
	}
	//go left
	if (Search(points, i, j - 1) != 0)
	{
		FindGroups(points, i, j - 1, group_number);
	}
	//go up
	if (Search(points, i - 1, j) != 0)
	{
		FindGroups(points, i - 1, j, group_number);
	}
}

// main find groups function
void FindGroups(std::vector<std::vector<int>>& points)
{
	// initialize global group number
	int group_number = 1;
	// loop through all entries
	for (int n = 0; n < points[0].size(); ++n)
	{
		// check if the point has already been assigned a group number
		if (points[2][n] == 0)
		{
			// find the coresponding group members for that point
			FindGroups(points, points[0][n], points[1][n], group_number);
			// increment the group number
			++group_number;
		}
	}
}








// calculate spatial moment
double m_ij(cie::DoubleMatrix& picture, int i, int j)
{
	double moment = 0.0;
	// preform loop over double sum
	for (int x = 0; x < picture.numberOfColumns(); ++x)
	{
		for (int y = 0; y < picture.numberOfRows(); ++y)
		{
			moment += std::pow(x,i) * std::pow(y,j) * picture.at(x, y);
		}
	}
	return moment;
}

// calculate central moment
double mu_ij(cie::DoubleMatrix& picture, int i, int j)
{
	// spatial moments
	double m00 = m_ij(picture, 0, 0); 
	double m01 = m_ij(picture, 0, 1);
	double m10 = m_ij(picture, 1, 0);
	// mass center
	double x_center = m10 / m00;
	double y_center = m01 / m00;


	double moment = 0.0;
	// preform loop over double sum
	for (int x = 0; x < picture.numberOfColumns(); ++x)
	{
		for (int y = 0; y < picture.numberOfRows(); ++y)
		{
			moment += std::pow(x - x_center, i) * std::pow(y - y_center, j) * picture.at(x, y);
		}
	}
	return moment;
}

// normalized central moment
double nu_ij(cie::DoubleMatrix& picture, int i, int j)
{
	double m00 = m_ij(picture, 0, 0);
	double mu = mu_ij(picture, i, j);

	double exponent = 1 + 0.5 * (i + j);

	return mu / std::pow(m00, exponent);
}

// First Invariant
std::array<double, 7 > HU_moments(cie::DoubleMatrix& picture)
{
	// calculate normalized central moments
	double nu02 = nu_ij(picture, 0, 2);
	double nu03 = nu_ij(picture, 0, 3);
	double nu11 = nu_ij(picture, 1, 1);
	double nu12 = nu_ij(picture, 1, 2);
	double nu20 = nu_ij(picture, 2, 0);
	double nu21 = nu_ij(picture, 2, 1);
	double nu30 = nu_ij(picture, 3, 0);
	
	// Calculate individual invariants
	double HU_1 = nu20 + nu02;
	double HU_2 = (nu20 - nu02) * (nu20 - nu02) + 4 * nu11 * nu11;
	double HU_3 = (nu30 - 3 * nu12) * (nu30 - 3 * nu12) + (3 * nu21 - nu03) * (3 * nu21 - nu03);
	double HU_4 = (nu30 + nu12) * (nu30 + nu12) + (nu21 + nu03) * (nu21 + nu03);
	double HU_5 = (nu30 - 3 * nu12) * (nu30 + nu12) * ((nu30 + nu12) * (nu30 + nu12) - 3 * (nu21 + nu03) * (nu21 + nu03)) + (3 * nu21 - nu03) * (nu21 + nu03) * (3 * (nu30 + nu12) * (nu30 + nu12) - (nu21 + nu03) * (nu21 + nu03));
	double HU_6 = (nu20 - nu02) * ((nu30 + nu12) * (nu30 + nu12) - (nu21 + nu03) * (nu21 + nu03)) + 4 * nu11 * (nu30 + nu12) * (nu21 + nu03);
	double HU_7 = (3 * nu21 - nu03) * (nu30 + nu12) * ((nu30 + nu12) * (nu30 + nu12) - 3 * (nu21 + nu03) * (nu21 + nu03)) - (nu30 - 3 * nu12) * (nu21 + nu03) * (3 * (nu30 + nu12) * (nu30 + nu12) - (nu21 + nu03) * (nu21 + nu03));

	// return result
	return { HU_1, HU_2, HU_3, HU_4, HU_5, HU_6, HU_7 };
}



// function that returns normalized log10 value of the number
double log10_normalization(double value)
{
	if (value > 0.0)
	{
		return std::log10(abs(value));
	}
	if (value < 0.0)
	{
		return -std::log10(abs(value));
	}
	if (value == 0.0)
	{
		return 0.0;
	}
}


// comparison method 2 from the sources - sum over the whole vector in log10 scale
// not sure what is the point, could be simpler and just comapare raw values
// use with extra caution when values are cose to zero
double Comparison_2(const std::array<double, 7>& HU_A, const std::array<double, 7>& HU_B)
{
	double result = 0.0;
	// loop over the whole HU_invariants array
	for (int i = 0; i < 7; ++i)
	{
		// normalize the values and add the together
		double m_A = log10_normalization(HU_A[i]);
		double m_B = log10_normalization(HU_B[i]);

		result += std::abs(m_A - m_B);
	}
	return result;
}

// comparison method 1 - just compares raw values
double Comparison_1(const std::array<double, 7>& HU_A, const std::array<double, 7>& HU_B)
{
	double result = 0.0;
	// loop over the whole HU_invariants array
	for (int i = 0; i < 7; ++i)
	{
		result += std::abs(HU_A[i] - HU_B[i]);
	}
	return result;
}



// evaluation method 2 from the sources 
double Comparison_2(const std::array<double, 7>& HU)
{
	double result = 0.0;
	// loop over the whole HU_invariants array
	for (int i = 0; i < 7; ++i)
	{
		// normalize the values and add the together
		double m_A = log10_normalization(HU[i]);

		result += std::abs(m_A);
	}
	return result;
}


// print function for HU moments array
void PrintMoments(const std::array<double, 7>& moments)
{
	std::cout << "(";
	for (int i = 0; i < 6; ++i)
	{
		std::cout << moments[i] << ", ";
	}
	std::cout << moments[6] << ")" << std::endl;
}





// scaning functions

// function that finds if the point lies in the rectangle give by any
// combination of two edge point
bool IsInside(cie::VEKTOR test_point, cie::VEKTOR v1, cie::VEKTOR v2)
{
	// find min -> bottom left point
	double x1 = std::min(v1.x(), v2.x());
	double y1 = std::min(v1.y(), v2.y());
	// find max -> upper right point
	double x2 = std::max(v1.x(), v2.x());
	double y2 = std::max(v1.y(), v2.y());


	// check if it lies between x bounds
	if (x1 < test_point.x() && test_point.x() < x2)
	{
		// check if it lies between y bounds
		if (y1 < test_point.y() && test_point.y() < y2)
		{
			return true;
		}
		return false;
	}
	return false;
}








// function that checks if the points from the list of coordinates
// input list of point to check and two vectors of bounding box
bool IsInside(std::vector<double>& x_coordinates, std::vector<double> y_coordinates, cie::VEKTOR v1, cie::VEKTOR v2)
{
	// loop over the whole list of points
	for (int n = 0; n < x_coordinates.size(); ++n)
	{
		// create a test point
		cie::VEKTOR test_point(x_coordinates[n], y_coordinates[n]);

		// check each point if it is inside
		if (IsInside(test_point, v1, v2) == true)
		{
			// function returns true when first such point is found
			return true;
		}
	}
	return false;
}




// function that "scans" a sqare area of points around a specified origin and returns a matrix of 1s and 0s
cie::DoubleMatrix Scan(std::vector<double> x_coordinates, std::vector<double> y_coordinates, double size, int resolution, cie::VEKTOR center)
{
	// initialize a square matrix of results - zero intialized
	cie::DoubleMatrix result(resolution, resolution);

	// calculate the width of one bin
	double bin_width = size / resolution;

	// calculate the null vector of the matrix - position of (0,0)
	cie::VEKTOR center_offset(size / 2.0, size / 2.0);
	cie::VEKTOR null = center - center_offset;

	// WARNING - I am not sure about mirroring and orientation of the picture with
	// respect to the global coordiante system

	// scan along x-coordinate -> j
	for (int x = 0; x < resolution; ++x)
	{
		// scan along y-coordinate -> i
		for (int y = 0; y < resolution; ++y)
		{
			// calculate the v1 and v2 vector - min and max of data bin
			// scaning is row-wise along x-coordinate
			cie::VEKTOR v1(bin_width * y, bin_width * x);
			cie::VEKTOR v2(bin_width * (y + 1), bin_width * (x + 1));
			//cie::VEKTOR v1(bin_width * x, bin_width * y);
			//cie::VEKTOR v2(bin_width * (x + 1), bin_width * (y + 1));
			v1 = v1 + null;
			v2 = v2 + null;

			// call inside funcition that check the entire list if the point is inside
			if (IsInside(x_coordinates, y_coordinates, v1, v2) == true)
			{
				// invert the result as printing occurs in reverse order
				result.at(resolution - 1 - x, y) = 1;
			}
		}
	}
	return result;
}

// returns 1 if test point defined by (x,y) is inside rectangle defined by (x1_, x2_, y1_, y2_). Otherwise, returns 0
bool insideRectangle(double x, double y, double x1_, double x2_, double y1_, double y2_)
{
	return x >= x1_ && x <= x2_ && y >= y1_ && y <= y2_;
}

//// ------------------------------------------------------------------------------------------------------------
//// function that determines if the line segment defined by 2 points cuts the rectangle defined by 2 points
//// combination of two edge point
bool doesItCut(cie::VEKTOR vktLine1, cie::VEKTOR vktLine2, cie::VEKTOR vktRect1, cie::VEKTOR vktRect2)
{
	int numberOfDiscretization = 100;
	// find min -> bottom left point of bounding box
	double x1_rect = std::min(vktRect1.x(), vktRect2.x());
	double y1_rect = std::min(vktRect1.y(), vktRect2.y());
	// find max -> upper right point of bounding box
	double x2_rect = std::max(vktRect1.x(), vktRect2.x());
	double y2_rect = std::max(vktRect1.y(), vktRect2.y());

	// find min -> bottom left point of line
	double x1_line = std::min(vktLine1.x(), vktLine2.x());
	double y1_line = std::min(vktLine1.y(), vktLine2.y());
	// find max -> upper right point of line
	double x2_line = std::max(vktLine1.x(), vktLine2.x());
	double y2_line = std::max(vktLine1.y(), vktLine2.y());

	// define all edges of the bounding box - "pixel"
	//cie::VEKTOR bottomLeft(x1, y1);
	//cie::VEKTOR bottomRight(x2, y1);
	//cie::VEKTOR upperRight(x2, y2);
	//cie::VEKTOR upperLeft(x1, y2);

	// initialize x and y increment, used to discretize the line
	double x_increment = 0.;
	double y_increment = 0.;

	// find x and y increment by dividing x length and y length by desired # of discretizations
	x_increment = (x2_line - x1_line) / numberOfDiscretization;
	y_increment = (y2_line - y1_line) / numberOfDiscretization;

	// loop across discretizations and run "inside" function until an intersection is found
	for (int i = 0; i < numberOfDiscretization; ++i)
	{
		if (insideRectangle(x1_line, y1_line, x1_rect, x2_rect, y1_rect, y2_rect) == 1)
		{
			return true;
		}
		// increment x and y coordinate by the calculated increment
		x1_line += x_increment;
		y1_line += y_increment;
	}
	return false;

	//// check crossing for each line segment on a rectangle
	//if (comparisonFunction(vktLine1, vktLine2, bottomLeft, bottomRight) == true)
	//{
	//	return true;
	//}
	//if (comparisonFunction(vktLine1, vktLine2, bottomRight, upperRight) == true)
	//{
	//	return true;
	//}
	//if (comparisonFunction(vktLine1, vktLine2, upperRight, upperLeft) == true)
	//{
	//	return true;
	//}
	//if (comparisonFunction(vktLine1, vktLine2, upperLeft, bottomLeft) == true)
	//{
	//	return true;
	//}
	//// if they all fail return false for no crossing
	//return false;

}

// function that checks if the points from the list of coordinates
// input list of point to check and two vectors of bounding box
bool doesItCut(std::vector<cie::VEKTOR>& inputPolygon, cie::VEKTOR vkt1, cie::VEKTOR vkt2)
{
	// new variable so input can be passed as a reference
	std::vector<cie::VEKTOR> polygon = inputPolygon;
	// add the fist element to the last place to ease calculations - closed loop
	polygon.push_back(polygon[0]);

	// loop over the whole list of edges
	for (int n = 0; n < polygon.size() - 1; ++n)
	{
		// create a test line segment
		cie::VEKTOR vktLine1(polygon[n].x(), polygon[n].y());
		cie::VEKTOR vktLine2(polygon[n + 1].x(), polygon[n + 1].y());

		// check each point if it is inside
		if (doesItCut(vktLine1, vktLine2, vkt1, vkt2) == true)
		{
			// function returns true when first such point is found
			return true;
		}
	}
	return false;
}


// function that "scans" a square area of points around a specified origin and returns a matrix of 1s and 0s
// points are given in term of a counterclockwise or clockwise defined polygon
cie::DoubleMatrix Scan(std::vector<cie::VEKTOR>& inputPolygon, double size, int resolution, cie::VEKTOR center)
{
	// new variable so input can be passed as a reference
	std::vector<cie::VEKTOR> testPolygon = inputPolygon;
	// add the fist element to the last place to ease calculations - closed loop
	testPolygon.push_back(testPolygon[0]);

	// initialize a square matrix of results - zero intialized
	cie::DoubleMatrix result(resolution, resolution);

	// calculate the width of one bin
	double bin_width = size / resolution;

	// calculate the null vector of the matrix - position of (0,0)
	cie::VEKTOR center_offset(size / 2.0, size / 2.0);
	//cie::VEKTOR center_offset(3.5 , 1.5);

	cie::VEKTOR null = center - center_offset;

	// WARNING - I am not sure about mirroring and orientation of the picture with
	// respect to the global coordiante system

	// scan along x-coordinate -> j
	for (int x = 0; x < resolution; ++x)
	{
		// scan along y-coordinate -> i
		for (int y = 0; y < resolution; ++y)
		{
			// calculate the v1 and v2 vector - min and max of data bin
			// scaning is row-wise along x-coordinate
			cie::VEKTOR vkt1(bin_width * y, bin_width * x);
			cie::VEKTOR vkt2(bin_width * (y + 1), bin_width * (x + 1));
			//cie::VEKTOR vkt1(bin_width * x, bin_width * y);
			//cie::VEKTOR vkt2(bin_width * (x + 1), bin_width * (y + 1));
			vkt1 = vkt1 + null;
			vkt2 = vkt2 + null;

			// call inside funcition that check the entire list if the point is inside
			if (doesItCut(testPolygon, vkt1, vkt2) == true)
			{
				// invert the result as printing occurs in reverse order
				result.at(resolution - 1 - x, y) = 1;
				//result.at(x, y) = 1;
			}
			else 
			{
				result.at(resolution - 1 - x, y) = 0;
				//result.at(x, y) = 0;
			}
			
		}
	}
	return result;
}



 ////------------------------------------------------------------------------------------------------------------





// interpolate points based on linear interpolation
cie::VEKTOR interpolatePoints(cie::VEKTOR& vkt1, cie::VEKTOR& vkt2, double t)
{
	// calculate new coordinates
	double x_new = vkt1.x() * (1.0 - t) + vkt2.x() * t;
	double y_new = vkt1.y() * (1.0 - t) + vkt2.y() * t;
	// return result
	return cie::VEKTOR(x_new, y_new);
}

// linear interpolation function that adds points to a polygon
// user specifies how many points in total there should be with more or less equal division
std::vector<cie::VEKTOR> interpolatePoints(std::vector<cie::VEKTOR>& inputPolygon, int nTotalPoints)
{
	// new variable so input can be passed as a reference
	std::vector<cie::VEKTOR> polygon = inputPolygon;
	// add the fist element to the last place to ease calculations - closed loop
	polygon.push_back(polygon[0]);

	// size of the polygon
	int n = polygon.size();
	// calculate total length of the points
	double totalLength = 0.0;
	for (int i = 0; i < n - 1; ++i)
	{
		// define an edge as a vektor
		cie::VEKTOR vkt1 = polygon[i];
		cie::VEKTOR vkt2 = polygon[i + 1];
		cie::VEKTOR vktEdge = vkt2 - vkt1;
		// add lengths together
		totalLength += vktEdge.length();
	}

	// initialize results
	std::vector<cie::VEKTOR> result;
	// define distance between two points
	double spacing = totalLength / nTotalPoints;

	// create new points on each edge
	for (int i = 0; i < n - 1; ++i)
	{
		// define an edge as a vektor
		cie::VEKTOR vkt1 = polygon[i];
		cie::VEKTOR vkt2 = polygon[i + 1];
		cie::VEKTOR vktEdge = vkt2 - vkt1;
		//current length
		double length = vktEdge.length();
		// calculate increment
		double increment = spacing / length;
		// loop add points on each edge
		for (double t = 0.0; t < 1.0; t += increment)
		{
			result.push_back(interpolatePoints(vkt1, vkt2, t));
		}
	}
	return result;
}








// Monotone chain - Convex hull
// sorting complexity - O(n*log(n))
// convex hull complexity - O(n)
// custom compare function that sorts the points in lexiographical order
// first by X and then by Y
int compareVektors(cie::VEKTOR &vkt1, cie::VEKTOR &vkt2)
{
	// compare the x-coordinates
	if (vkt1.x() < vkt2.x())
	{
		return 1;
	}
	// check X equality and compare the y-coordinates
	if (vkt1.x() == vkt2.x() && vkt1.y() < vkt2.y())
	{
		return 1;
	}
	// if both fail return 0 -> switch the elements
	return 0;
}


// 3D cross product of OA and OB vectors
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
double calculateCrossProduct(cie::VEKTOR &O, cie::VEKTOR &A, cie::VEKTOR &B)
{
	return (A.x() - O.x()) * (B.y() - O.y()) - (A.y() - O.y()) * (B.x() - O.x());
}

// Andrew's monotone chain 
// https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#C++
// Returns a list of points on the convex hull in counter-clockwise order.
std::vector<cie::VEKTOR> calculateConvexHull(std::vector<cie::VEKTOR> points)
{
	int n = points.size();
	int k = 0;

	// if there are only 3 or less points, retun back these points back
	if (n <= 3)
	{
		return points;
	}

	// initialize results 
	//std::vector<cie::VEKTOR> hull(2 * n); // two time the size just to be sure - Optimize
	std::vector<cie::VEKTOR> hull(n);

	// Sort points lexicographically from smallest to biggest x
	// and if two entris have the same value x, sort them by y
	std::sort(points.begin(), points.end(), compareVektors);

	// I have no idead how exactly this work on memory manipulation
	// Build lower hull
	// loop from smallest to biggest X 
	for (int i = 0; i < n; ++i) 
	{
		while (k >= 2 && calculateCrossProduct(hull[k - 2], hull[k - 1], points[i]) <= 0)
		{
			k--;
		}
		hull[k++] = points[i];
	}

	// Build upper hull
	// loop from biggest to smallest X
	for (int i = n - 1, t = k + 1; i > 0; --i)
	{
		while (k >= t && calculateCrossProduct(hull[k - 2], hull[k - 1], points[i - 1]) <= 0)
		{
			k--;
		}
		hull[k++] = points[i - 1];
	}

	//resizes the vector
	hull.resize(k - 1); // last point is not the same
	//hull.resize(k); // the last point is the same as the first one

	return hull;
}



// points must be sequential in counterclockwise or clockwise direction
// calculate area of a polygon defined by coordinates
double calculatePolygonArea(std::vector<cie::VEKTOR>& coordinates)
{
	// intialize area
	double area = 0.0;
	// number of coordinates
	int n = coordinates.size();

	// Calculate value of shoelace formula
	int j = n - 1;
	for (int i = 0; i < n; i++)
	{
		area += (coordinates[j].x() + coordinates[i].x()) * (coordinates[j].y() - coordinates[i].y());
		j = i; // j is previous vertex to i 
	}
	// counterclockwise direction -> positive 
	// clockwise direction -> negative.
	// This is why we must return a absolute value absolute value 
	return std::abs(area * 0.5);
}

// circle defined by three points
// [0] -> x, [1] -> y, [2] -> radius squared,  
std::array<double, 3> calculateCircle(cie::VEKTOR& vkt1, cie::VEKTOR& vkt2, cie::VEKTOR& vkt3)
{
	// get the coordinates
	double x1 = vkt1.x();
	double y1 = vkt1.y();
	double x2 = vkt2.x();
	double y2 = vkt2.y();
	double x3 = vkt3.x();
	double y3 = vkt3.y();

	double a = x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2;
	double b = (x1 * x1 + y1 * y1) * (y3 - y2) + (x2 * x2 + y2 * y2) * (y1 - y3) + (x3 * x3 + y3 * y3) * (y2 - y1);
	double c = (x1 * x1 + y1 * y1) * (x2 - x3) + (x2 * x2 + y2 * y2) * (x3 - x1) + (x3 * x3 + y3 * y3) * (x1 - x2);
	//double d = (x1 * x1 + y1 * y1) * (x3 * y2 - x2 * y3) + (x2 * x2 + y2 * y2) * (x1 * y3 - x3 * y1) + (x3 * x3 + y3 * y3) * (x2 * y1 - x1 * y2);

	double x = -b / (2 * a);
	double y = -c / (2 * a);
	double rSquared = (x - x1) * (x - x1) + (y - y1) * (y - y1);
	//double rSquared = (b * b + c * c - 4 * a * d) / (4 * a * a);

	return { x, y, rSquared };
}

// calculate angle between two vectors
// angle starts from vkt1 in counterclockwise direction
// retun values are positive for 0...pi and negative for pi...2*pi
double calculateAngle(cie::VEKTOR& vkt1, cie::VEKTOR& vkt2)
{
	double dotProduct = vkt1.x() * vkt2.x() + vkt1.y() * vkt2.y();
	double determinant = vkt1.x() * vkt2.y() - vkt1.y() * vkt2.x();
	return std::atan2(determinant, dotProduct);
}



// rotation of coordinates counterclockwise through an angle about the origin in 2D space
std::vector<cie::VEKTOR> rotatePoints(std::vector<cie::VEKTOR> & coordinates, double angle)
{
	// initialize the results
	std::vector<cie::VEKTOR> newCoordinates(coordinates.size());
	// loop over the coordinates
	for (int i = 0; i < coordinates.size(); ++i)
	{
		// get old coordiantes
		double x = coordinates[i].x();
		double y = coordinates[i].y();
		// calculate new coordiantes - simple rotation matrix
		double x_new = x * std::cos(angle) - y * std::sin(angle);
		double y_new = x * std::sin(angle) + y * std::cos(angle);
		// add to the vector of results
		newCoordinates[i] = cie::VEKTOR(x_new, y_new);
	}
	return newCoordinates;
}

// calculate the minimum bounding box, function takes in the convex hull coordinates
// convex hull must be defined in counterclockwise direction
std::vector<cie::VEKTOR> calculateMinimalRectangle(std::vector<cie::VEKTOR>& inputConvexHull)
{
	// new variable so input can be passed as a reference
	std::vector<cie::VEKTOR> convexHull = inputConvexHull;
	// add the fist element to the last place to ease calculations - closed loop
	convexHull.push_back(convexHull[0]);

	// size of the convex hull
	int h = convexHull.size();
	// define area and angle
	double min_area = INFINITY;
	double min_angle = 0.0;
	// define end results
	cie::VEKTOR bottomLeft;
	cie::VEKTOR upperLeft;
	cie::VEKTOR upperRight;
	cie::VEKTOR bottomRight;
	// define x-axis as a vector
	cie::VEKTOR vktX(1.0, 0.0);

	// loop over the whole convex hull and transform coordinates
	// we know the fact that minimum rectangle will be colinear with one edge of the convex hull
	for (int i = 0; i < h - 1; ++i)
	{
		// define an edge as a vektor
		double xEdge = convexHull[i + 1].x() - convexHull[i].x();
		double yEdge = convexHull[i + 1].y() - convexHull[i].y();
		cie::VEKTOR vktEdge(xEdge, yEdge);

		// calclulate angle between edge and x-aixs
		double angle = calculateAngle(vktX, vktEdge);

		// define and calculate new rotated coordiantes, angle must be reversed
		std::vector<cie::VEKTOR> rotatedCoordiantes = rotatePoints(convexHull, -angle);

		// seperate to x and y coordinates to ease find max_function
		std::vector<double> x_coordinates(h);
		std::vector<double> y_coordinates(h);
		for (int j = 0; j < h; ++j)
		{
			x_coordinates[j] = rotatedCoordiantes[j].x();
			y_coordinates[j] = rotatedCoordiantes[j].y();
		}

		// find min and max x
		double min_x = *std::min_element(x_coordinates.begin(), x_coordinates.end());
		double max_x = *std::max_element(x_coordinates.begin(), x_coordinates.end());

		// find min and max y
		double min_y = y_coordinates[i]; // min y is always the active i
		double max_y = *std::max_element(y_coordinates.begin(), y_coordinates.end());
		
		// define witdh and height of the bounding rectangle
		double width = max_x - min_x;
		double height = max_y - min_y;

		// calculate area of the bounding rectangle
		double area = width * height;

		// update the minimal area and coordiantes of actual rectangle
		if (area < min_area)
		{
			min_area = area;
			min_angle = angle;
			// save coordinates of rectangle in rotated CS
			bottomLeft  = cie::VEKTOR(min_x, min_y);
			bottomRight = cie::VEKTOR(max_x, min_y);
			upperRight  = cie::VEKTOR(max_x, max_y);
			upperLeft   = cie::VEKTOR(min_x, max_y);		
		}
	}
	// intialize rotated vector
	std::vector<cie::VEKTOR> rectangle(4);
	rectangle[0] = bottomLeft;
	rectangle[1] = bottomRight;
	rectangle[2] = upperRight;
	rectangle[3] = upperLeft;
	
	// rotated points back in the original plane
	// return results in counterclockwise direction
	return rotatePoints(rectangle, min_angle);
}




int main()
{
	//// input data as a matrix that contains the number of point in each bin
	cie::DoubleMatrix data(8, 8);

	data.at(0, 3) = 1;
	data.at(1, 3) = 1;
	data.at(2, 1) = 1;
	data.at(2, 3) = 1;
	data.at(3, 1) = 1;
	data.at(3, 2) = 1;
	data.at(3, 3) = 1;

	data.at(1, 5) = 1;
	data.at(2, 5) = 1;
	data.at(2, 6) = 1;
	data.at(2, 7) = 1;

	data.at(5, 2) = 1;
	data.at(6, 2) = 1;
	data.at(6, 3) = 1;
	data.at(7, 2) = 1;

	data.at(4, 5) = 1;
	data.at(5, 4) = 1;
	data.at(5, 5) = 1;
	data.at(5, 6) = 1;
	data.at(6, 6) = 1;
	data.at(5, 7) = 1;
	data.at(2, 0) = 1;

	// print the data
	cie::print(data);

	// converts the data into a list of interesting points
	// initialize rows and column indicies and group number (zero for unchanged data)
	std::vector<int> rows_i;
	std::vector<int> columns_j;
	std::vector<int> GroupNumber;

	// create a vector of size 3
	std::vector<std::vector<int>> InterestingPoints(3);

	// scan the matrix to find non-zero entries (find points of interest)
	for (int i = 0; i < data.numberOfRows(); ++i)
	{
		for (int j = 0; j < data.numberOfColumns(); ++j)
		{
			if (data.at(i, j) != 0)
			{
				// store indices
				rows_i.push_back(i); 
				columns_j.push_back(j);
				// store default group number
				GroupNumber.push_back(0);
			}
		}
	}

	InterestingPoints[0] = rows_i;
	InterestingPoints[1] = columns_j;
	InterestingPoints[2] = GroupNumber;

	// find groups of points
	FindGroups(InterestingPoints);

	// print function for InterestingPoints
	for (int i = 0; i < InterestingPoints[0].size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			std::cout << InterestingPoints[j][i] << ",";
		}
		std::cout << std::endl;
	}







	// shape matching in image recognition
	// https://docs.opencv.org/2.4/modules/imgproc/doc/structural_analysis_and_shape_descriptors.html?highlight=cvmatchshapes#humoments
	// https://en.wikipedia.org/wiki/Image_moment
	// https://www.youtube.com/watch?v=O-hCEXi3ymU


	// create a binary picture of 20x20
	cie::DoubleMatrix picture_1(20, 20);
	// create a square 
	for (int i = 4; i < 16; ++i)
	{
		picture_1.at(i, 4) = 1;
		picture_1.at(i, 5) = 1;
		picture_1.at(i, 14) = 1;
		picture_1.at(i, 15) = 1;
	}
	for (int i = 6; i < 14; ++i)
	{
		picture_1.at(4, i) = 1;
		picture_1.at(5, i) = 1;
		picture_1.at(14, i) = 1;
		picture_1.at(15, i) = 1;
	}
	// print the picture_1
	cie::print_dots(picture_1);


	// create a binary picture of 20x20
	cie::DoubleMatrix picture_2(20, 20);
	// create a square 
	for (int i = 0; i < 12; ++i)
	{
		picture_2.at(i, 0) = 1;
		picture_2.at(i, 1) = 1;
		picture_2.at(i, 10) = 1;
		picture_2.at(i, 11) = 1;
	}
	for (int i = 2; i < 10; ++i)
	{
		picture_2.at(0, i) = 1;
		picture_2.at(1, i) = 1;
		picture_2.at(10, i) = 1;
		picture_2.at(11, i) = 1;
	}
	// print the picture_2
	cie::print_dots(picture_2);


	// create a binary picture of 20x20
	cie::DoubleMatrix picture_3(10, 10);
	// create a square 
	for (int i = 4; i < 10; ++i)
	{
		picture_3.at(i, 0) = 1;
		picture_3.at(i, 5) = 1;
	}
	for (int i = 1; i < 5; ++i)
	{
		picture_3.at(4, i) = 1;
		picture_3.at(9, i) = 1;
		
	}
	// print the picture_3
	cie::print_dots(picture_3);


	// create a binary picture of 20x20
	cie::DoubleMatrix picture_4(20, 20);
	// create a square 
	for (int i = 1; i < 9; ++i)
	{
		picture_4.at(i, 10 - i) = 1;
		picture_4.at(i + 1, 10 - i) = 1;
		picture_4.at(i + 2, 10 - i) = 1;
		
		picture_4.at(i + 1, 9 + i) = 1;
		picture_4.at(i + 2, 9 + i) = 1;
		picture_4.at(i + 3, 9 + i) = 1;

		picture_4.at(18 - i, 10 - i) = 1;
		picture_4.at(17 - i, 10 - i) = 1;
		picture_4.at(16 - i, 10 - i) = 1;

		picture_4.at(17 - i, 9 + i) = 1;
		picture_4.at(16 - i, 9 + i) = 1;
		picture_4.at(15 - i, 9 + i) = 1;
	}
	picture_4.at(9, 1) = 1;
	picture_4.at(10, 17) = 0;
	picture_4.at(11, 17) = 0;
	picture_4.at(8, 17) = 0;
	picture_4.at(7, 17) = 0;
	// print the picture_4
	cie::print_dots(picture_4);



	// calculate HU moments for picture 1
	std::array<double, 7> HU_1 = HU_moments(picture_1);
	PrintMoments(HU_1);
	
	// calculate HU moments for picture 2
	std::array<double, 7> HU_2 = HU_moments(picture_2);
	PrintMoments(HU_2);

	// calculate HU moments for picture 3
	std::array<double, 7> HU_3 = HU_moments(picture_3);
	PrintMoments(HU_3);

	// calculate HU moments for picture 4
	std::array<double, 7> HU_4 = HU_moments(picture_4);
	PrintMoments(HU_4);

	// comparison between 1 and 4
	// I don't know exact sensiability of this comparison function
	std::cout << "Normal comparison 3-4: " << Comparison_1(HU_3, HU_4) << std::endl;
	std::cout << "Log comparison 3-4: " << Comparison_2(HU_3, HU_4) << std::endl;

	







	// intialize vectors of coordiantes
	std::vector<double> x_coordinates{ 5.5, 5.5, 5.5, 5.5,  5.5,  6.5,  7.5,  8.5,  9.5,  10.5, 11.5, 12.5, 12.5, 12.5, 12.5, 12.5, 11.5, 10.5, 9.5, 8.5, 7.5, 6.5 };
	std::vector<double> y_coordinates{ 7.5, 8.5, 9.5, 10.5, 11.5, 11.5, 11.5, 11.5, 11.5, 11.5, 11.5, 11.5, 10.5, 9.5,  8.5,  7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5 };

	//std::vector<double> x_coordinates{ 5.1, 5.2, 5.3, 5.4, 5.5, 6.1, 7.2, 8.3, 9.4, 10.5, 11.1, 12.2, 12.3, 12.4, 12.5, 12.1, 11.2, 10.3, 9.4, 8.5, 7.1, 6.2 };
	//std::vector<double> y_coordinates{ 7.1, 8.2, 9.3, 10.4, 11.5, 11.1, 11.2, 11.3, 11.4, 11.5, 11.1, 11.2, 10.3, 9.4, 8.5, 7.1, 7.2, 7.3, 7.4, 7.5, 7.1, 7.2 };

	// center of scan
	cie::VEKTOR center(9.0, 9.0);

	cie::DoubleMatrix picture_5 = Scan(x_coordinates, y_coordinates, 20, 20, center);
	cie::print_dots(picture_5);

	cie::DoubleMatrix picture_6 = Scan(x_coordinates, y_coordinates, 10, 5, center);
	cie::print_dots(picture_6);

	// instantiate four coordinates which form the four corners of a rectangle
	cie::VEKTOR bottomLeftPoint(5.5,7.5);
	cie::VEKTOR topLeftPoint(5.5, 11.5);
	cie::VEKTOR topRightPoint(12.5, 11.5);
	cie::VEKTOR bottomRightPoint(12.5, 7.5);

	// construct an "inputPolygon" with a vector of coordinates. NOTE: must be constructed in a
	// consecutive manner! Can't construct top, bottom, left, right side for instance
	std::vector<cie::VEKTOR> inputPolygon;
	inputPolygon.push_back(topLeftPoint);
	inputPolygon.push_back(bottomLeftPoint);
	inputPolygon.push_back(bottomRightPoint);
	inputPolygon.push_back(topRightPoint);

	std::cout << "The next matrix is the result of the line cutting algorithm: " << std::endl;

	// run the algorithm on the constructed polygon and print the result
	cie::DoubleMatrix picture_7 = Scan(inputPolygon,40,40,center);
	cie::print_dots(picture_7);
	
	// calculate HU moments for picture 5
	std::array<double, 7> HU_5 = HU_moments(picture_5);
	PrintMoments(HU_5);

	// calculate HU moments for picture 6
	std::array<double, 7> HU_6 = HU_moments(picture_6);
	PrintMoments(HU_6);

	// comparison between 5 and 6
	// treshold value must be selected upon experience
	std::cout << "Normal comparison : " << Comparison_1(HU_5, HU_6) << std::endl;
	std::cout << "Log comparison : " << Comparison_2(HU_5, HU_6) << std::endl;

	std::cout << "Log comparison : " << Comparison_2(HU_6) - Comparison_2(HU_5) << std::endl;


	

	// test data that produces axis alligned minimal rectangle
	/*std::vector<double> x_coordinates_1{ 3, 3, 5, 6, 6, 7, 8, 9,10,11,11};
	std::vector<double> y_coordinates_1{ 6, 9, 3, 5, 8, 4, 7, 6, 9, 3, 7};*/

	// testa data that produces rotated minimal rectangle
	std::vector<double> x_coordinates_1{ 2, 3, 3, 4, 6, 7, 8 ,7, 9};
	std::vector<double> y_coordinates_1{ 6, 3, 8, 5, 5, 4, 3, 1, 5 };

	// transform the coordinates into vector of VEKTOR
	std::vector<cie::VEKTOR> arrayOfPoints(x_coordinates_1.size());
	for (int i = 0; i < x_coordinates_1.size(); i++)
	{
		cie::VEKTOR vktPoint(x_coordinates_1[i], y_coordinates_1[i]);

		arrayOfPoints[i] = vktPoint;
		std::cout << "(" << vktPoint.x() << ", " << vktPoint.y() << ") ";
	}
	std::cout << std::endl;








	// calculate the convex hull
	std::vector<cie::VEKTOR> hull = calculateConvexHull(arrayOfPoints);
	// print the results
	for (int i = 0; i < hull.size(); ++i)
	{
		std::cout << "(" << hull[i].x() << ", " << hull[i].y() <<  ") ";
	}
	std::cout << std::endl;




	std::cout << calculatePolygonArea(arrayOfPoints) << std::endl;

	/*
	std::array<double, 3> circle =  calculateCircle(arrayOfPoints[0], arrayOfPoints[1], arrayOfPoints[2]);
	std::cout << "x = " << circle[0] << ", y = " << circle[1] << ", r^2 = " << circle[2] << std::endl;
	*/

	std::cout << "angle = " << calculateAngle(cie::VEKTOR(1, 0), cie::VEKTOR(4, -1)) << std::endl;


	cie::VEKTOR vkt1(1, 4);
	cie::VEKTOR vkt2(7, 1);

	cie::VEKTOR interpolate = interpolatePoints(vkt1, vkt2, 0.66666666);
	std::cout << "x = " << interpolate.x() << ", y = " << interpolate.y() << std::endl;


	// calcualte the minimum rectangle
	std::vector<cie::VEKTOR> minimumRectangle = calculateMinimalRectangle(hull);
	// print the results
	for (int i = 0; i < minimumRectangle.size(); ++i)
	{
		std::cout << "(" << minimumRectangle[i].x() << ", " << minimumRectangle[i].y() << ") ";
	}
	std::cout << std::endl;

	std::cout << "Area of original convex hull: " << calculatePolygonArea(hull) << std::endl;
	std::cout << "Area of minimum rectangle: " << calculatePolygonArea(minimumRectangle) << std::endl;

	// divide points
	std::vector<cie::VEKTOR> morePoints = interpolatePoints(minimumRectangle, 10);
	// print the results
	for (int i = 0; i < morePoints.size(); ++i)
	{
		std::cout << "(" << morePoints[i].x() << ", " << morePoints[i].y() << ") ";
	}
	std::cout << std::endl;
	std::cout << "size = " << morePoints.size() << std::endl;



	// test for faulty vector
	/*std::cout << std::endl;
	cie::VEKTOR a(10.0, 5.0);
	cie::VEKTOR b(5.0, 2.0);
	cie::VEKTOR c = a + b;
	cie::print(a);
	cie::print(b);
	cie::print(c);
*/





	return 0;
}