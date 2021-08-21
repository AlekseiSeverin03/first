#include <stdio.h>
#include <math.h>
#include <assert.h>

#define MANY 3

int SolveLinealEquration(double b, double c, double *x1)        /* функция, решающая линейное уравнение (b * x - c = 0, где а = 0), х1 - корень этого уравнения; значение самой функции равно количесту корней уравнения. */
{
	if (b == 0 && c == 0)
		return MANY;
	else if (b == 0)
		return 0;
	else 
	{
		*x1 = -c / b;
		return 1;
	}
}

void SolveSquareEquration(double a, double b, double c, double *x1, double *x2, int *nroots)      /* функция, решеающая квадратное уравнение (a * (x^2) + b * x + c = 0), х1 и х2- его корни. */
{

	assert (x1 != NULL);
	assert (x2 != NULL);
	assert (nroots !=  NULL);

	if (a == 0)
	{
		*nroots = SolveLinealEquration(b, c, x1);
	}
	else
	{
		double discriminant = b * b - 4 * a * c;
		if (discriminant < 0)
		{
			*nroots = 0;
		}
		else if (discriminant == 0)
		{
			*x1 = -b / (2 * a);
			*nroots = 1;
		}
		else 
		{
			*x1 = (-b - sqrt(discriminant)) / (2 * a);
			*x2 = (-b + sqrt(discriminant)) / (2 * a);
			*nroots = 2;
		}
	}
}

int main()
{
	printf("Solver square equation\n\n");
	printf("Enter 3 numbers - coefficients square equatio: a b c\n\n");
	
	double a = 0, b = 0, c = 0, x1, x2;
	int NSuccessfullyProcessedValues = 0, nroots = -1;      /* nroots - переменная, в которую будет записано количество корней квадратного уравнения. */
	
	NSuccessfullyProcessedValues = scanf("%lg %lg %lg", &a, &b, &c);   /* ввод коэффициентов квадратного уравнения. */
	if (NSuccessfullyProcessedValues != 3) 
	{
		printf("Error (wrong input)\n");
		return 1;
	}

	SolveSquareEquration(a, b, c, &x1, &x2, &nroots);

	switch (nroots)   /* вывод ответа (корней квадратного уравнения). */
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
		printf("Error (crash in the program)\n");
		return 2;
	}
	return 0;
}
