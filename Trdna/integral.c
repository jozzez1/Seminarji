#include <stdio.h>
#include <math.h>
#include <string.h>

static const double ef = 5;
static const double t  = 1;

typedef struct
{
	int M;
	
	double qx,
	       qy,
	       w,
	       beta,
	       mu,
	       kin;
} Data;

void init_Data (int M,
		double qx,
		double qy,
		double w,
		double beta,
		double mu,
		double kin,
		Data * dat)
{
	dat->M		= M;
	dat->qx		= qx;
	dat->qy		= qy;
	dat->w		= w;
	dat->beta	= beta;
	dat->mu		= mu;
	dat->kin	= kin;
}

void free_Data (Data * dat)
{
	free (dat);
}

typedef struct
{
	double bmin,
	       bmax,
	       dbeta,
	       wmin,
	       wmax,
	       dw;
	
	Data * dat;
} Control;

void init_Control (int M,
		   double bmin,
		   double bmax,
		   double dbeta
		   double wmin,
		   double wmax,
		   double dw,
		   double qx,
		   double qy,
		   double mu,
		   double kin,
		   Control * control)
{
	control->dat = (Data *) malloc (sizeof (Data));
	init_Data (M, qx, qy, wmin, bmin, mu, kin, control->dat);
	
	control->bmin	= bmin;
	control->bmax	= bmax;
	control->dbeta	= dbeta;
	control->wmin	= wmin;
	control->wmax	= wmax;
	control->dw	= dw;
}

void free_Control (Control * control)
{
	free (control->dat);
	free (control);
}

double fI1 (double kx,
           double ky,
	   double qx,
	   double qy,
	   double w,
	   double beta,
	   double mu,
	   double kin)
{
	double pred	= 0.5 / kin,
	       f	= exp ((-2) * beta * (cos (kx) + cos (ky)) - beta * mu) + 1;
	       g	= cos (kx) + cos (ky) +
				sin (kx) * sin (qx) + sin (ky) * sin (qy) -
				cos (kx) * cos (qx) - cos (ky) * cos (qy) -
				w/(2 * kin),
	       re	= pred / (f * g);
	
	return re;
}

double fI2 (double kx,
            double ky,
	    double qx,
	    double qy,
	    double w,
	    double beta,
	    double mu,
	    double kin)
{
	double pred	= (-0.5) / kin,
	       f	= exp ((-2) * kin * beta * (cos (kx)*cos (qx) +
				cos (ky)*cos (qy) - sin (kx)*sin(qx) -
				sin (ky)*sin (qy)) - beta * mu) + 1;
	       g	= cos (kx) + cos (ky) +
				sin (kx) * sin (qx) + sin (ky) * sin (qy) -
				cos (kx) * cos (qx) - cos (ky) * cos (qy) -
				w/(2 * kin),
	       re	= pred / (f * g);
	
	return re;
}

double f (int i, int j, Data * dat)
{
	double U    = M_PI / (dat->M - 1),
	       D    = 0.5 * U,
	       kx   = i * U + D,
	       ky   = j * U + D;
	
	double f1 = fI1 (kx, ky, 
		dat->qx, dat->qy, dat->w, dat->beta, dat->mu, dat->kin),
	       f2 = fI2 (kx, ky,
		dat->qx, dat->qy, dat->w, dat->beta, dat->mu, dat->kin);
	
	return f1 - f2;
}

double integrate (Data * dat)
{
	double U	= M_PI / dat->M,
	       A	= U * U,
	       sum	= 0.0;
	       
	int i, j;
	for (i = 0; i <= dat->M - 1; i++)
	{
		for (j = 0; j <= dat->M - 1; j++)
			sum += f (i, j, dat);
	}
	
	sum *= A;
	
	return sum;
}

void loop-plotter (Control * control, FILE * fout)
{
	double q = sqrt (pow (control->dat->qx,2) + pow (control->dat->qy,2));
	
	do
	{
		do
		{
			double integral = integrate (control->dat);
			control->dat->beta += control->dbeta;
			
			fprintf (fout, "% 9.6lf % 9.6lf % 9.6lf\n",
				control->dat->w,
				control->dat->beta,
				1 - integral/(q*q));
			
		} while (control->dat->beta <= control->bmax);
		
		control->dat->beta	 = control->bmin;
		control->dat->w		+= control->dw;
		fprintf (fout, "\n");
	
	} while (control->dat->w <= control->wmax);
}

void readfile (char * name, double * a)
{
	FILE * fin = fopen (name, "r");
	
	int p	= 0,
	    c	= 0;
	    
	while (1)
	{
		c = fgetc (fin);
		if (c == '#')
		{
			do
			{
				c = fgetc (fin);
			} while (c != '\n' || c == '#');
		}
		
		else if (p < 2)
		{
			fseek (fin, -1, SEEK_CUR);
			fscanf (fin, "%lf %lf %lf %lf %lf %lf %lf %lf",
					&a[0], &a[1], &a[2], &a[3],
					&a[4], &a[5], &a[6], &a[7]);
			p++;
		}
		
		if (p == 2 || c == EOF)
			break;
	}

	fclose (fin);
}

void calculate (int M,
		double qx,
		double qy,
		char * name
		char * dout)
{
	Control * control = (Control *) malloc (sizeof (Control));
	double * a = (double *) malloc (8 * sizeof (double));

	readfile (name, a);

	init_Control (M, a[0], a[1], a[2], a[3], a[4], a[5],
			qx, qy, mu, kin, control);

	free (a);

	FILE * fout = fopen (dout, "w");
	loop-plotter (control, fout);

	fclose (fout);
	free_Control (control);
}

int main (int argc, char ** argv)
{
	int M	= 100;
}










