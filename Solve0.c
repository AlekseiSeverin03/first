#include <stdio.h>
#include <math.h>
#include <assert.h>

enum OUTPUT_TEST
{
	FALSE = 0,
	TRUE
};

enum number_roots
{
	ZERO_ROOTS = 0,
	ONE_ROOT, 
	TWO_ROOTS, 
	INF_ROOTS
};

const double PRECISION = 0.00001;

enum error_codes 
{
	OK = 0,             	//< OK CODE
	ERROR_WRONG_INPUT,      //< USER EN
	ERROR_IN_TESTS, 
};

//----------------------------------------------------------
//! Checking if two value are equal
//! 
//! @param [in]  value1  first value
//! @param [in]  value2  second value
//!
//! @return result comparison (TRUE or FALSE)
//!
//! @note function can compare NAN (NAN = NAN)  
//-----------------------------------------------------------
int cmp_with_number_on_equality(double value1, double value2)
{
	if (!isnan(value1) && !isnan(value2))
	{
		return (fabs(value1 - value2) <= PRECISION);
	}
	if (isnan(value1) && isnan(value2))
	{
		return TRUE;
	}
	return FALSE;
}

//----------------------------------------------------------
//! Checking if first value less second value
//! 
//! @param [in]  value1  first value
//! @param [in]  value2  second value
//!
//! @return result comparison (TRUE or FALSE) 
//-----------------------------------------------------------

int cmp_with_number_on_sign_less(double value1, double value2)
{
	return value1 - value2 < -PRECISION;
}

//--------------------------------------------------------------------------------------------
//! Solve linear equation b*x + c = 0
//! 
//! @param [in]  b  b-coefficient
//! @param [in]  c  c-coefficient
//! @param [out] x  pointer to root
//!
//! @return number of roots
//!
//! @note if linear equation have one solution then funtion will change value under pointer x
//--------------------------------------------------------------------------------------------

int SolveLinearEquation(double b, double c, double *x) 
{
	assert(x != NULL);
	assert(!isnan(b));
	assert(!isnan(c));

	if (cmp_with_number_on_equality(b, 0) && cmp_with_number_on_equality(c, 0))
	{
		return INF_ROOTS;
	}
	else if (cmp_with_number_on_equality(b, 0))
	{
		return ZERO_ROOTS;
	}
	*x = -c / b;
	return ONE_ROOT;
}

//------------------------------------------------
//! Testing function SolveLinearEquation
//!
//! @param [in]  testNum   test number
//! @param [in]  b         b-coefficient
//! @param [in]  c         c-coefficient
//! @param [in]  nRootsRef true number roots
//! @param [in]  xref      true solution
//! 
//! @return result test (TRUE or FALSE)
//------------------------------------------------

int UnitTestSolveLinear(int testNum, double b, double c, int nRootsRef, double xref) 
{
	double x = NAN;
	int nRoots = SolveLinearEquation(b, c, &x);
	if (nRoots != nRootsRef || !cmp_with_number_on_equality(x, xref))
	{
		printf("Test %d\n"
			   "nRoots = %d, x = %lg\n"
			   "Should be: %d, %lg\n\n", testNum, nRoots, x, nRootsRef, xref);
		return TRUE;
	}
	return FALSE;
}

//----------------------------------------------------------------
//! Solve square equation using discriminant a*x^2 + b*x + c = 0
//!  
//! @param [in]  a  a-coefficient
//! @param [in]  b  b-coefficient
//! @param [in]  c  c-coefficient
//! @param [out] x1 poonter to first root
//! @param [out] x2 pointer to second root
//!
//! @return number of roots
//!
//! @note discriminant solution
//----------------------------------------------------------------

int SolveSquareEquation_Discrim(double a, double b, double c, double *x1, double *x2)
{
	double discriminant = b * b - 4 * a * c;
	if (cmp_with_number_on_sign_less(discriminant, 0))
    {
         return ZERO_ROOTS;
    }
    else if (cmp_with_number_on_equality(discriminant, 0))
    {
         *x1 = -b / (2 * a);
         return ONE_ROOT;
    }
    double sqrt_d = sqrt(discriminant);
    *x1 = (-b - sqrt_d) / (2 * a);
    *x2 = (-b + sqrt_d) / (2 * a);
    return TWO_ROOTS;       
}

//--------------------------------------------------------------------------------------------------------
//! Solve square equation a*x^2 + b*x + c = 0
//!
//! @param [in]  a  a-coefficient
//! @param [in]  b  b-coefficient
//! @param [in]  c  c-coefficient
//! @param [out] x1 poonter to first root
//! @param [out] x2 pointer to second root
//!
//! @return number of roots
//!
//! @note if square equation have one or two solutions then funtion will change value under pointer x1; 
//!       if square equation have two solution then funtion will change value under pointer x2; 
//!       if equation have infinite number of roots then RETURN: INF_ROOTS
//--------------------------------------------------------------------------------------------------------

int SolveSquareEquation(double a, double b, double c, double *x1, double *x2) 
{
	assert(x1 != NULL);
	assert(x2 != NULL);
	assert(x1 != x2);
	
	assert(!isnan(a));
	assert(!isnan(b));
	assert(!isnan(c));

	if (cmp_with_number_on_equality(a, 0))
	{
		return SolveLinearEquation(b, c, x1);
	}
	else if (cmp_with_number_on_equality(b, 0) && cmp_with_number_on_equality(c, 0))
	{
		*x1 = 0;
		return ONE_ROOT;
	}
	else if (cmp_with_number_on_equality(b, 0))
	{
		if (c > 0)
		{
			return ZERO_ROOTS;
		}
		*x1 = -sqrt(fabs(-c / a));
		*x2 = -*x1;
		return TWO_ROOTS;
	}
	else if (cmp_with_number_on_equality(c, 0))
	{
		*x1 = -b / a;
		*x2 = 0;
		return TWO_ROOTS;
	}
	return SolveSquareEquation_Discrim(a, b, c, x1, x2);
}	

//-------------------------------------------------------------------------
//! Testing function SolveSquareEquation and SolveSquareEquation_Discrim
//!
//! @param [in]  testNum   test number
//! @param [in]  a         a-coefficient
//! @param [in]  b         b-coefficient
//! @param [in]  c         c-coefficient
//! @param [in]  nRootsRef true number roots
//! @param [in]  x1ref     true first solution
//! @param [in]  x2ref     true second solution
//! 
//! @return result test (TRUE or FALSE)
//--------------------------------------------------------------------------

int UnitTestSolveSquareEquation(int testNum, double a, double b, double c, int nRootsRef, double x1ref, double x2ref)
{
	double x1 = NAN, x2 = NAN;
	int nRoots = SolveSquareEquation(a, b, c, &x1, &x2);
	if (nRoots != nRootsRef || !cmp_with_number_on_equality(x1, x1ref) || !cmp_with_number_on_equality(x2, x2ref))
	{
		printf("Test %d\n"
			   "nRoots = %d, x1 = %lg, x2 = %lg\n"
			   "Should be: %d, %lg, %lg\n\n", testNum, nRoots, x1, x2, nRootsRef, x1ref, x2ref);
		return TRUE;
	}
	return FALSE;
} 

//--------------------------------
//! Runing tests
//!
//! @return number failed tests
//--------------------------------

int RunUnitTests(void)
{
	int fail = 0;
//                                  testNUM  a  b  c  nRootsRef   x1ref  x2ref
	if (UnitTestSolveSquareEquation(1,       0, 0, 0, INF_ROOTS,  NAN,   NAN)) fail++;
	if (UnitTestSolveSquareEquation(2,       0, 0, 1, ZERO_ROOTS, NAN,   NAN)) fail++;
	if (UnitTestSolveSquareEquation(3,       0, 1, 0, ONE_ROOT,   0,     NAN)) fail++;
	if (UnitTestSolveSquareEquation(4,       0, 1, 1, ONE_ROOT,   -1,    NAN)) fail++;
	if (UnitTestSolveSquareEquation(5,       1, 0, 0, ONE_ROOT,   0,     NAN)) fail++;
	if (UnitTestSolveSquareEquation(6,       1, 0, 1, ZERO_ROOTS, NAN,   NAN)) fail++;	
	if (UnitTestSolveSquareEquation(7,       1, 0, -1,TWO_ROOTS,  -1,    1))   fail++;
	if (UnitTestSolveSquareEquation(8,       1, 1, 0, TWO_ROOTS,  -1,    0))   fail++;
	if (UnitTestSolveSquareEquation(9,       1, 1, 1, ZERO_ROOTS, NAN,   NAN)) fail++;
	if (UnitTestSolveSquareEquation(10,      1, 2, 1, ONE_ROOT,   -1,    NAN)) fail++;
	if (UnitTestSolveSquareEquation(11,      1, 5, 6, TWO_ROOTS,  -3,    -2))  fail++;
	return fail;
}
	
int main()
{
	int count_fail = RunUnitTests();	
	if (count_fail)
	{
		printf("NUMBER ERROR IN UnitTEST: %d\n\n", count_fail);
		return ERROR_IN_TESTS;
	}

	printf("Square equation solver: a * x ^ 2 + b * x + c = 0\n");
	printf("Enter 3 numbers - coefficients square equation: a b c\n\n");

	double a = 0, b = 0, c = 0, x1 = 0, x2 = 0;
	int n = scanf("%lg %lg %lg", &a, &b, &c);
	
	if (n != 3) 
	{
		printf("Error (wrong input)\n");
		return ERROR_WRONG_INPUT;
	}

	/* nroots - количество корней квадратного уравнения. */
        int nroots = SolveSquareEquation(a, b, c, &x1, &x2);

	switch (nroots) 
	{
		case ZERO_ROOTS: {
			printf("No roots\n");
			break;
		}
		case ONE_ROOT: {
			printf("x = %lg\n", x1);
			break;
		}
		case TWO_ROOTS: {
			printf("x1 = %lg , x2 = %lg\n", x1, x2);
			break;
		}
		case INF_ROOTS: {
			printf("Any number\n");
			break;
		}
		default: {
			/* неожиданное количество корней уравнения */
			assert(0);
		}
	}
	return OK;
}
