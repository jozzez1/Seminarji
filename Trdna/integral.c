#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>

static const double ef = 4 / (M_PI * M_PI);
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

	int fix_beta,
	    fix_w;
	
	Data * dat;
} Control;

void init_Control (int M,
		   int fbeta,
		   double bmin,
		   double bmax,
		   double dbeta,
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

	control->fix_beta = fbeta;
	control->fix_w    = 0;
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
	       f	= exp ((-2) * kin * beta * (cos (kx) + cos (ky) - mu)) + 1,
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
				sin (ky)*sin (qy) - mu)) + 1,
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
	
	return f1 + f2;
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
	
	sum *= 4*A;
	
	return sum;
}

void progress_bar (int a, int b, double dsec)
{
	double left;
	if (dsec)
	{
		left = (b - a) * dsec/b;
		printf ("ETA:% 3.2lfs | ", left);
	}

	int percent = 20 * a/b,
	    i;

	printf ("DONE: % d/%d (% 3.0lf%%) [", a, b, 100.0 * a/b);

	for (i = 0; i <= percent - 1; i++)
		printf ("=");
	printf (">");

	for (i = percent; i<= 18; i++)
		printf (" ");
	printf ("]\n\033[F\033[J");
}

void loop_plotter (Control * control, FILE * fout)
{
	double q2 = pow (control->dat->qx,2) + pow (control->dat->qy,2);

	int a	= 0,	// counter as to how far we've come
	    b	= (int) (((control->bmax - control->bmin)/control->dbeta) *
			    ((control->wmax - control->wmin)/control->dw));

	if (control->fix_beta)
		b = (int) ((control->wmax - control->wmin)/control->dw);
	
	time_t start,
	       stop;

	double dtime;

	do
	{
		do
		{
			time (&start);
			double integral = integrate (control->dat);
			
			fprintf (fout, "% 9.6lf % 9.6lf % 9.6lf\n",
				control->dat->w,
				control->dat->beta,
				1 + integral/q2);

			control->dat->w	+= control->dw;

			time (&stop);
			dtime = difftime (stop, start);

			progress_bar (a, b, dtime);
			a++;
			
		} while (control->dat->w <= control->wmax);
		
		if (control->fix_beta)
			break;

		control->dat->w	 = control->wmin;
		control->dat->beta += control->dbeta;
		fprintf (fout, "\n");
	
	} while (control->dat->beta <= control->bmax);
}

void readfile (char * name, double * a)
{
	FILE * fin = fopen (name, "r");

	if (fin == NULL)
	{
		printf ("No input file!\n");
		exit (EXIT_FAILURE);
	}
	
	int p	= 0,
	    gg	= 0,
	    c	= 0;
	    
	while (1)
	{
		c = fgetc (fin);
		if (c == '#')
		{
			do
			{
				c = fgetc (fin);
				gg++;
			} while (c != '\n' || c == '#');
		}
		
		else if (p < 2 && gg != 0)
		{
			fseek (fin, -1, SEEK_CUR);
			fscanf (fin, "%lf %lf %lf %lf %lf %lf",
					&a[0], &a[1], &a[2], &a[3], &a[4], &a[5]);
			p++;
		}
		
		if (p == 2 || c == EOF)
			break;
	}

	fclose (fin);
}

void calculate (int M,
		int fbeta,
		double qx,
		double qy,
		double mu,
		double kin,
		char * name,
		char * dout)
{
	Control * control = (Control *) malloc (sizeof (Control));
	double * a = (double *) malloc (8 * sizeof (double));

	readfile (name, a);

	init_Control (M, fbeta, a[0], a[1], a[2], a[3], a[4], a[5],
			qx, qy, mu, kin, control);

	free (a);

	FILE * fout = fopen (dout, "w");
	loop_plotter (control, fout);

	fclose (fout);
	free_Control (control);
}

int main (int argc, char ** argv)
{
	int M	= 100,
	    fbeta = 0,
	    arg;

	double qx	= M_PI,
	       qy	= M_PI,
	       kin	= t,
	       mu	= ef;


	while ((arg = getopt (argc, argv, "M:x:y:t:e:fh")) != -1)
	{
		switch (arg)
		{
			case 'M':
				M = atoi (optarg);
				break;
			case 'x':
				qx = atof (optarg);
				break;
			case 'y':
				qy = atof (optarg);
				break;
			case 't':
				kin = atof (optarg);
				break;
			case 'e':
				mu = atof (optarg);
				break;
			case 'f':
				fbeta = 1;
				break;
			case 'h':
				printf ("-M <int>\n");
				printf ("-x <double>: qx\n");
				printf ("-y <double>: qy\n");
				printf ("-t <double>: kinetic constant\n");
				printf ("-e <double>: mu -- fermi energy\n");
				printf ("-h -- print this list\n");
				exit (EXIT_SUCCESS);
			default:
				printf ("Unknown command!\n");
				exit (EXIT_FAILURE);
		}
	}

	calculate (M, fbeta, qx, qy, mu, kin, "input.txt", "output.txt");
	exit (EXIT_SUCCESS);
}

