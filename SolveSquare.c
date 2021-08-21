#include <stdio.h>
#include <math.h>

#define MANY 3

float discriminant(float a, float b, float c, int *sign)
{
	float discrim = b * b - 4 * a * c;
	if (discrim < 0) 
	{
		*sign = -1;
		return discrim;
	}
	else
	{
		if (discrim == 0)
		{
			*sign = 0;
			return 0;
		}
		else
		{
			*sign = 1;
			return discrim;
		}
	}
}

int main()
{
	float a, b, c, x1, x2;
	int n, sign, nroots;
	printf("Enter 3 numbers: a b c\n");
	n = scanf("%f %f %f", &a, &b, &c);
	if (n != 3) 
	{
		printf("Error (wrong input)\n");
		return 4;
	}
	if (a == 0) 
	{
		if (b == 0)
		{
			if (c == 0)
			{
				nroots = MANY;
			} 
			else
			{
				nroots = 0;
			}
		}
		else
		{
			nroots = 1;
			x1 = -c / b;
		}
	}
	else
	{
		float d = discriminant(a, b, c, &sign);
		switch (sign) 
		{
		case -1:
			nroots = 0;
			break;
		case 0:
			nroots = 1;
			x1 = -b / (2 * a);
			break;
		case 1: 
			nroots = 2;
			x1 = (-b - sqrt(d)) / (2 * a);
			x2 = (-b + sqrt(d)) / (2 * a);
			break;
		default:
			printf("Error (in discriminant)\n");
			return 2;
		}
	}

	switch (nroots) 
	{
	case 0: 
		printf("No roots\n");
		break;
	case 1:
		printf("x = %f\n", x1);
		break;
	case 2:
		printf("x1 = %f , x2 = %f\n", x1, x2);
		break;
	case MANY:
		printf("Any number\n");
		break;
	default:
		printf("Error\n");
		return 1;
	}
	return 0;
}
