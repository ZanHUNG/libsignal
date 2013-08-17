#ifndef _LIBWFS_H_
#define _LIBWFS_H_

#define _USE_MATH_DEFINES
#include <math.h>

#include "libwfs_types.h"

extern int bohman(double * out, int n);
extern int bartlett(double * out, int n);

#endif