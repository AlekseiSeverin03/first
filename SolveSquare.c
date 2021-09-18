#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "input.h"

enum OUTPUT_TESTS
{
	FALSE = 0,
	TRUE
};

enum NUMBER_ROOTS
{
	ZERO_ROOTS = 0,
	ONE_ROOT, 
	TWO_ROOTS, 
	INF_ROOTS
};

const double PRECISION = 0.00001;

enum OUTPUT_fMAIN 
{
	OK = 0,             	
	ERROR_WRONG_INPUT,      
	ERROR_IN_TESTS, 
};

//--------------------------------------------------------------------------------------------------------
//! Input three double 
//! 
//! @param [in] 1st double
//! @param [in] 2nd double
//! @param [in] 3rd double
//!
//! @return void
//!
//! @note from function you can terminate main program, for this you need to press 'y' 
//--------------------------------------------------------------------------------------------------------

void InputThreeDouble (double *a, double *b, double *c);

//--------------------------------------------------------------------------------------------------------
//! Output answer task
//!
//! @param [in] nroots number roots of equation
//!
//! @return void
//!
//! @note Output in stdout roots square equation
//--------------------------------------------------------------------------------------------------------

void OutputAnswer (int nroots, double *x1, double *x2)
{
	switch (nroots) 
	{
		case ZERO_ROOTS: 
		{
			printf("No roots\n");
			break;
		}

		case ONE_ROOT: 
		{
			printf("x = %lg\n", *x1);
			break;
		}

		case TWO_ROOTS: 
		{
			printf("x1 = %lg , x2 = %lg\n", *x1, *x2);
			break;
		}

		case INF_ROOTS: 
		{
			printf("Any number\n");
			break;
		}

		default: 
		{
			/* неожиданное количество корней уравнения */
			assert(0);
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//! Checking if two value are equal
//! 
//! @param [in]  value1  first value
//! @param [in]  value2  second value
//!
//! @return result comparison (TRUE or FALSE)
//!
//! @note function can compare NAN (NAN = NAN)  
//--------------------------------------------------------------------------------------------------------

int isEqual (double value1, double value2)
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

//--------------------------------------------------------------------------------------------------------
//! Checking if first value less second value
//! 
//! @param [in]  value1  first value
//! @param [in]  value2  second value
//!
//! @return result comparison (TRUE or FALSE) 
//--------------------------------------------------------------------------------------------------------

int isLess (double value1, double value2)
{
	assert (!isnan(value1));
	assert (!isnan(value2));

	return value1 - value2 < -PRECISION;
}

//--------------------------------------------------------------------------------------------------------
//! Function find max number of 2 double
//!
//! @param [in] x 1st double 
//! @param [in] y 2nd double
//!
//! @return max of 2 double
//--------------------------------------------------------------------------------------------------------

double MaxModule_2_Double (double x, double y)
{
	return fabs (x) > fabs (y) ? x : y;
}

///--------------------------------------------------------------------------------------------------------
//! Function find max number of 3 double
//!
//! @param [in] a 1st double 
//! @param [in] b 2nd double
//! @param [in] c 3rd double
//!
//! @return max of 3 double
//--------------------------------------------------------------------------------------------------------

double MaxModule_3_Double (double a, double b, double c)
{
	return MaxModule_2_Double (MaxModule_2_Double (a, b), MaxModule_2_Double (a, c));
}

//--------------------------------------------------------------------------------------------------------
//! Division of small coefficients by the max them
//!
//! @param [in] a 1st coef
//! @param [in] b 2nd coef
//! @param [in] c 3rd coef
//!
//! @return void
//!
//! @note numbers is "small" if they less const PRECISION
//--------------------------------------------------------------------------------------------------------

void DivSmallCoefByMax (double *a, double *b, double *c)
{
	if (isEqual (*a, 0) && isEqual (*b, 0) && isEqual (*c, 0))
	{
		double max = MaxModule_3_Double (*a, *b, *c);
		*a /= max;
		*b /= max;
		*c /= max;
	}
}

//--------------------------------------------------------------------------------------------------------
//! Solve linear equation b*x + c = 0
//! 
//! @param [in]  b  b-coefficient
//! @param [in]  c  c-coefficient
//! @param [out] x  pointer to root
//!
//! @return number of roots
//!
//! @note if linear equation have one solution then funtion will change value under pointer x
//--------------------------------------------------------------------------------------------------------

int SolveLinearEquation (double b, double c, double *x) 
{
	assert(x != NULL);

	assert(!isnan(b));
	assert(!isnan(c));

	if (isEqual (b, 0) && isEqual (c, 0))
	{
		return INF_ROOTS;
	}
	else if (isEqual (b, 0))
	{
		return ZERO_ROOTS;
	}

	*x = -c / b;
	
	return ONE_ROOT;
}

//--------------------------------------------------------------------------------------------------------
//! Testing function SolveLinearEquation
//!
//! @param [in]  testNum   test number
//! @param [in]  b         b-coefficient
//! @param [in]  c         c-coefficient
//! @param [in]  nRootsRef true number roots
//! @param [in]  xref      true solution
//! 
//! @return result test (TRUE or FALSE)
//--------------------------------------------------------------------------------------------------------

int UnitTestSolveLinear (int testNum, double b, double c, int nRootsRef, double xref) 
{
	double x = NAN;
	int nRoots = SolveLinearEquation (b, c, &x);

	if (nRoots != nRootsRef || !isEqual (x, xref))
	{
		printf("Test %d\n"
			   "nRoots = %d, x = %lg\n"
			   "Should be: %d, %lg\n\n", testNum, nRoots, x, nRootsRef, xref);

		return TRUE;
	}

	return FALSE;
}

//--------------------------------------------------------------------------------------------------------
//! Solve square equation using discriminant a*x^2 + b*x + c = 0
//!  
//! @param [in]  a  a-coefficient
//! @param [in]  b  b-coefficient
//! @param [in]  c  c-coefficient
//! @param [out] x1 pointer to first root
//! @param [out] x2 pointer to second root
//!
//! @return number of roots
//!
//! @note discriminant solution
//--------------------------------------------------------------------------------------------------------

int SolveSquareEquation_Discrim (double a, double b, double c, double *x1, double *x2)
{
	assert (x1 != NULL);
	assert (x2 != NULL);

	assert (!isnan(a));
	assert (!isnan(b));
	assert (!isnan(c));

	double discriminant = b * b - 4 * a * c;
	if (isLess (discriminant, 0))
    {
         return ZERO_ROOTS;
    }
    else if (isEqual (discriminant, 0))
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
//! @param [out] x1 pointer to first root
//! @param [out] x2 pointer to second root
//!
//! @return number of roots
//!
//! @note if square equation have one or two solutions then funtion will change value under pointer x1; 
//!       if square equation have two solution then funtion will change value under pointer x2; 
//!       if equation have infinite number of roots then RETURN: INF_ROOTS
//--------------------------------------------------------------------------------------------------------

int SolveSquareEquation (double a, double b, double c, double *x1, double *x2) 
{
	assert (x1 != NULL);
	assert (x2 != NULL);
	assert (x1 != x2);
	
	assert (!isnan(a));
	assert (!isnan(b));
	assert (!isnan(c));

	if (isEqual (a, 0))
	{
		return SolveLinearEquation (b, c, x1);
	}
	else if (isEqual (b, 0) && isEqual (c, 0))
	{
		*x1 = 0;
		return ONE_ROOT;
	}
	else if (isEqual (b, 0))
	{
		if (c > 0)
		{
			return ZERO_ROOTS;
		}
		*x1 = -sqrt(fabs(-c / a));
		*x2 = -*x1;
		return TWO_ROOTS;
	}
	else if (isEqual (c, 0))
	{
		*x1 = -b / a;
		*x2 = 0;
		return TWO_ROOTS;
	}
	
	return SolveSquareEquation_Discrim (a, b, c, x1, x2);
}	

//--------------------------------------------------------------------------------------------------------
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
//--------------------------------------------------------------------------------------------------------

int UnitTestSolveSquareEquation (int testNum, double a, double b, double c, int nRootsRef, double x1ref, double x2ref)
{
	double x1 = NAN, x2 = NAN;
	int nRoots = SolveSquareEquation (a, b, c, &x1, &x2);
	
	if (nRoots != nRootsRef || !isEqual (x1, x1ref) || !isEqual (x2, x2ref))
	{
		printf("Test %d\n"
			   "nRoots = %d, x1 = %lg, x2 = %lg\n"
			   "Should be: %d, %lg, %lg\n\n", testNum, nRoots, x1, x2, nRootsRef, x1ref, x2ref);
		
		return TRUE;
	}
	
	return FALSE;
} 

//--------------------------------------------------------------------------------------------------------
//! Runing tests
//!
//! @return number failed tests
//--------------------------------------------------------------------------------------------------------

int RunUnitTests (void)
{
	int fail = 0;
//                                   testNUM  a  b  c  nRootsRef   x1ref  x2ref
	if (UnitTestSolveSquareEquation (1,       0, 0, 0, INF_ROOTS,  NAN,   NAN)) fail++;
	if (UnitTestSolveSquareEquation (2,       0, 0, 1, ZERO_ROOTS, NAN,   NAN)) fail++;
	if (UnitTestSolveSquareEquation (3,       0, 1, 0, ONE_ROOT,   0,     NAN)) fail++;
	if (UnitTestSolveSquareEquation (4,       0, 1, 1, ONE_ROOT,   -1,    NAN)) fail++;
	if (UnitTestSolveSquareEquation (5,       1, 0, 0, ONE_ROOT,   0,     NAN)) fail++;
	if (UnitTestSolveSquareEquation (6,       1, 0, 1, ZERO_ROOTS, NAN,   NAN)) fail++;	
	if (UnitTestSolveSquareEquation (7,       1, 0, -1,TWO_ROOTS,  -1,    1))   fail++;
	if (UnitTestSolveSquareEquation (8,       1, 1, 0, TWO_ROOTS,  -1,    0))   fail++;
	if (UnitTestSolveSquareEquation (9,       1, 1, 1, ZERO_ROOTS, NAN,   NAN)) fail++;
	if (UnitTestSolveSquareEquation (10,      1, 2, 1, ONE_ROOT,   -1,    NAN)) fail++;
	if (UnitTestSolveSquareEquation (11,      1, 5, 6, TWO_ROOTS,  -3,    -2))  fail++;
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
		
	printf("\nThis program solves square equation: a * x ^ 2 + b * x + c = 0\n");

	double a = 0, b = 0, c = 0;

	InputThreeDouble (&a, &b, &c);

	DivSmallCoefByMax (&a, &b, &c);
	
	double x1 = 0, x2 = 0;
	
	int nroots = SolveSquareEquation (a, b, c, &x1, &x2);

	OutputAnswer (nroots, &x1, &x2);

	return OK;
}
