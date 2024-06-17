#include "UI.h"
#include "Tests.h"
#include "ArrayRectangularMatrix.h"
#include "RectangularMatrix.h"
#include "SquareMatrix.h"
#include "Output.h"
#include "Sequence.h"


int main()
{
	TEST();
	TestarSquareMatrix();
	
	TestarRectangularMatrix();
	std::cout<<std::endl;
	MatrixMenu();
}