#include <stdio.h>
#include <math.h>
#include <assert.h>

#define MANY 3

/*
 *  функция, решающая линейное уравнение (b * x + c = 0, где а = 0),
 *  х - корень этого уравнения; значение самой функции равно количесту корней уравнения;
 *  при работе функция буддет изменять значение под указателем х,
 *  если линейное уравнение имеет одно решение. 
 */
int SolveLinearEquation(double b, double c, double *x) 
{
	assert(x != NULL);
	if (isnan(b) == 0 && isnan(c) == 0) {

	if (b == 0 && c == 0)
	{
		return MANY;
	}
	else if (b == 0)
	{
		return 0;
	}
	*x = -c / b;
	return 1;

	}
	return -1;
}
/* 
 *  функция, решающая квадратное уравнение (a * (x^2) + b * x + c = 0), х1 и х2- его корни; 
 *  значение самой функции равно количеству корней уравнения;
 *  при работе функция будет изменять значение под указателем х1,
 *  если квадратное уравнение будет иметь 1 или 2 корня, а под указателем х2,
 *  если квадратное уравнение имеет ровно 2 корня.
 */
int SolveSquareEquation(double a, double b, double c, double *x1, double *x2) 
{

	assert(x1 != NULL);
	assert(x2 != NULL);
	
	if (isnan(a) == 0 && isnan(b) == 0 && isnan(c) == 0) {	
	
	if (a == 0)
	{
		return SolveLinearEquation(b, c, x1);
	}
	else if (b == 0 && c == 0)
	{
		*x1 = 0;
		return 1;
	}
	else if (b == 0)
	{
		if (c > 0)
		{
			return 0;
		}
		*x1 = sqrt(fabs(-c / a));
		*x2 = -sqrt((-c / a));
		return 2;
	}
	else if (c == 0)
	{
		*x1 = 0;
		*x2 = -b / a;
		return 2;
	}
	
	{
		double discriminant = b * b - 4 * a * c;
		if (discriminant < 0)
		{
			return 0;
		}
		else if (discriminant == 0)
		{
			*x1 = -b / (2 * a);
			return 1;
		}
		else 
		{
			*x1 = (-b - sqrt(discriminant)) / (2 * a);
			*x2 = (-b + sqrt(discriminant)) / (2 * a);
			return 2;
		}
	}
	
	}
	return -1;
}

int main()
{
	printf("Solver square equation\n\n");
	printf("Enter 3 numbers - coefficients square equatio: a b c\n\n");
	
	double a = 0, b = 0, c = 0, x1 = 0, x2 = 0;
	/*
	 * nroots - переменная, в которую будет записано количество корней квадратного уравнения.
	 */
	int n_successfully_processed_values = 0, nroots = -1;  
	
	n_successfully_processed_values = scanf("%lg %lg %lg", &a, &b, &c);   
	if (n_successfully_processed_values != 3) 
	{
		printf("Error (wrong input)\n");
		return 1;
	}

	nroots = SolveSquareEquation(a, b, c, &x1, &x2);

	switch (nroots) 
	{
	case 0: 
		printf("No roots\n");
		break;
	case 1:
		printf("x = %lg\n", x1);
		break;
	case 2:
		printf("x1 = %lg , x2 = %lg\n", x1, x2);
		break;
	case MANY:
		printf("Any number\n");
		break;
	default:
		printf("Error (got unexpected nroots value)\n");
		return 2;
	}
	return 0;
}
