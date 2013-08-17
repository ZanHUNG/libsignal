#include "../include/libwfs.h"

int bartlett(double * out, int n)
{
	int i;
	int n_1, n_2;
	
	if( n < 1 )
		return WFS_NO_ERROR;
	
	if( !out )
	{
		return WFS_INVALID_OUTPUT;
	}
	
	if( n == 1 )
	{
		out[0] = 1.0;
		return WFS_NO_ERROR;
	}
	else
	{
		n_2 = n - 1;
		n_1 = n_2 / 2;
		for( i = 1; i <= n_1; i++ )
		{
			out[i] = i * 2.0 / n_2;
		}
		for( ; i < n_2; i++ )
		{
			out[i] = 2.0 - i * 2.0 / n_2;
		}
		
		out[0] = 0.0;
		out[n-1] = 0.0;
	}
	
	return WFS_NO_ERROR;
}