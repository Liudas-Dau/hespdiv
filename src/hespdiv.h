#ifndef R_hespdiv_H
#define R_hespdiv_H

#ifdef hespdiv_XPORT
# define hespdiv_PREFIX(name) SP_XPORT(name)
#else
# define hespdiv_PREFIX(name) name
#endif
/* remember to touch local_stubs.c */

#define hespdiv_VERSION "1.1"

#include <R.h>
/* RSB 091203 */
#include <Rdefines.h>
#define R_OFFSET 1
#include <Rinternals.h>
#include <Rmath.h>

/* from insiders.c 

int pipbb(double pt1, double pt2, double *bbs);
int between(double x, double low, double up); 
SEXP insiders(SEXP n1, SEXP bbs); */

/* from pip.c */

#ifndef MIN
# define MIN(a,b) ((a)>(b)?(b):(a))
#endif
#ifndef MAX
# define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#define BUFSIZE 8192

/* polygon structs: */
typedef struct {
	double		x, y;
} PLOT_POINT;

typedef struct {
	PLOT_POINT	min, max;
} MBR;

typedef struct polygon {
	MBR mbr;
	int lines;
	PLOT_POINT	*p;
    int close; /* 1 - is closed polygon */
} POLYGON;

void setup_poly_minmax(POLYGON *pl);
char InPoly(PLOT_POINT q, POLYGON *Poly);
SEXP R_point_in_polygon_sp(SEXP px, SEXP py, SEXP polx, SEXP poly);

/* RSB 091203 */

#define DIM     2               /* Dimension of points */
typedef double  tPointd[DIM];   /* type double point */

#endif
/* remember to touch local_stubs.c */