#define _USE_MATH_DEFINES
#include <math.h>

#include "../include/libwfs.h"

int bohman(double * out, int n)
{
	int i;
	int n_1;
	double x;
	
	if( n < 1 )
		return 0;
	
	if( !out )
	{
		return WFS_INVALID_OUTPUT;
	}
	
	if( n == 1 )
	{
		out[0] = 1.0;
		return 0;
	}
	else
	{
		n_1 = n - 1;
		for( i = 1; i < n_1; i++ )
		{
			x = fabs( -1.0 + i * 2.0/n_1 );
			out[i] = ( 1 - x ) * cos( M_PI * x ) + M_1_PI * sin( M_PI * x );
		}
		
		out[0] = 0.0;
		out[n_1] = 0.0;
	}
	
	return WFS_NO_ERROR;
}