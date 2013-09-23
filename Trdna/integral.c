#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>

static const double t  = 1;

typedef struct
{
	int M;

	double qx,
	       qy,
	       w,
	       beta,
	       mu,
	       kin,
	       alpha;
} Data;

void init_Data (int M,
		double qx,
		double qy,
		double w,
		double beta,
		double mu,
		double kin,
		double alpha,
		Data * dat)
{
	dat->M		= M;
	dat->qx		= qx;
	dat->qy		= qy;
	dat->w		= w;
	dat->beta	= beta;
	dat->mu		= mu;
	dat->kin	= kin;
	dat->alpha	= alpha;
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
	    fix_wmega;
	
	Data * dat;
} Control;

void init_Control (int M,
		   int fbeta,
		   int fwmega,
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
		   double alpha,
		   Control * control)
{
	control->dat = (Data *) malloc (sizeof (Data));
	init_Data (M, qx, qy, wmin, bmin, mu, kin, alpha, control->dat);
	
	control->bmin	= bmin;
	control->bmax	= bmax;
	control->dbeta	= dbeta;
	control->wmin	= wmin;
	control->wmax	= wmax;
	control->dw	= dw;

	control->fix_beta  = fbeta;
	control->fix_wmega = fwmega;
}

void free_Control (Control * control)
{
	free (control->dat);
	free (control);
}

double other_f_for_mi (int i, int j, Data * dat, double mu)
{
	double U	= M_PI / (dat->M - 1),
	       D	= 0.5 * U,
	       kx	= i*U + D,
	       ky	= j*U + D;

	double cx	= cos (kx),
	       cy	= cos (ky),
	       b	= dat->beta,
	       f1	= exp ((-2) * b *(cx + cy) - b*mu),
	       f	= 1.0 / (f1 + 1);

	return f;
}

double crap_integrate_for_mu (Data * dat, double mu)
{
	double U	= M_PI / dat->M,
	       A	= U*U,
	       sum	= 0.0;

	int i, j;
	for (i = 0; i <= dat->M - 1; i++)
	{
		for (j = 0; j <= dat->M - 1; j++)
			sum += other_f_for_mi (i, j, dat, mu);
	}

	sum *= 4*A;
	sum -= 0.5 * (M_PI * M_PI);

	return sum;
}

double find_mu (Data * dat)
{
	double prec	= 0.001,
	       mu1	= 10,
	       mu2	= 5;

	do
	{
		double I2	= crap_integrate_for_mu (dat, mu2),
		       I1	= crap_integrate_for_mu (dat, mu1),
		       mu3	= mu2 - (I2 * ((mu2 - mu1)/(I2 - I1)));

		mu1 = mu2;
		mu2 = mu3;
	} while (fabs (mu2 - mu1) > prec);
/*
	if (isnan (mu1))
	{
		exit (EXIT_FAILURE);
		fprintf (stderr, "Sorry, but mu2 doesn't converge.\n");
	}

	if (isinf (mu1))
	{
		exit (EXIT_FAILURE);
		fprintf (stderr, "Sorry, but mu2 returns infinity.\n");
	}
*/
	return mu1;
}

double f (int i, int j, Data * dat, double mu)
{
	double U    = M_PI / (dat->M - 1),
	       D    = 0.5 * U,
	       kx   = i * U + D,
	       ky   = j * U + D;
	
	double cx	= cos (kx),
	       cy	= cos (ky),
	       w	= dat->w,
	       b	= dat->beta,
	       g1	= 2 * (cx + cy) - w/2,
	       g2	= dat->alpha*0.25,
	       g	= g1 / (g1*g1 + g2*g2),
	       e1	= exp ((-2) * b * (cx + cy) - b * mu),	// ef == 0
	       e2	= exp (2 * b * (cx + cy) - b * mu),	// ef == 0
	       f1	= 1.0/(e1 + 1),
	       f2	= 1.0/(e2 + 1),
	       pred	= 1.0/(4 * M_PI * M_PI),
	       re	= pred * g * (f1 - f2);

	return re;
}

double integrate (Data * dat, double mu)
{
	double U	= M_PI / dat->M,
	       A	= U * U,
	       sum	= 0.0;
	       
	int i, j;
	for (i = 0; i <= dat->M - 1; i++)
	{
		for (j = 0; j <= dat->M - 1; j++)
			sum += f (i, j, dat, mu);
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
	double mu = find_mu (control->dat);

	int a	= 0,	// counter as to how far we've come
	    b	= (int) (((control->bmax - control->bmin)/control->dbeta) *
			    ((control->wmax - control->wmin)/control->dw));

	if (control->fix_beta)
		b = (int) ((control->wmax - control->wmin)/control->dw);

	if (!control->fix_wmega)
	{
		do
		{
			do
			{
				double integral = integrate (control->dat, mu);

				fprintf (fout, "% 9.6lf % 9.6lf % 9.6lf %9.6lf\n",
						control->dat->w,
						control->dat->beta,
						integral,
						mu);

				control->dat->w	+= control->dw;
				progress_bar (a, b, 0);

				a++;

			} while (control->dat->w <= control->wmax);

			if (control->fix_beta)
				break;

			control->dat->w	 = control->wmin;
			control->dat->beta += control->dbeta;
			mu = find_mu (control->dat);
			fprintf (fout, "\n");
		
		} while (control->dat->beta <= control->bmax);
	}

	else if (control->fix_wmega)
	{
		b = (int) ((control->bmax - control->bmin)/control->dbeta);
		
		do
		{
			double integral = integrate (control->dat, mu);

			fprintf (fout, "% 9.6lf % 9.6lf % 9.6lf %9.6lf\n",
					control->dat->beta,
					control->dat->w,
					integral,
					mu);

			control->dat->beta += control->dbeta;
			mu = find_mu (control->dat);
			progress_bar (a, b, 0);

			a++;
		} while (control->dat->beta <= control->bmax);
	}
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
		int fwmega,
		double qx,
		double qy,
		double mu,
		double kin,
		double alpha,
		char * name,
		char * dout)
{
	Control * control = (Control *) malloc (sizeof (Control));
	double * a = (double *) malloc (8 * sizeof (double));

	readfile (name, a);

	init_Control (M, fbeta, fwmega, a[0], a[1], a[2], a[3], a[4], a[5],
			qx, qy, mu, kin, alpha, control);

	free (a);

	FILE * fout = fopen (dout, "w");
	loop_plotter (control, fout);

	fclose (fout);
	free_Control (control);
}

int main (int argc, char ** argv)
{
	int M	   = 100,
	    fbeta  = 0,
	    fwmega = 0,
	    arg;

	double qx	= M_PI,
	       qy	= M_PI,
	       kin	= t,
	       a	= 1.0 / M;


	while ((arg = getopt (argc, argv, "M:x:y:t:a:fwh")) != -1)
	{
		switch (arg)
		{
			case 'M':
				M = atoi (optarg);
				a = 12.5 / M;
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
			case 'f':
				fbeta = 1;
				break;
			case 'w':
				fwmega = 1;
				break;
			case 'a':
				a = atof (optarg);
				break;
			case 'h':
				printf ("-M <int>\n");
				printf ("-x <double>: qx\n");
				printf ("-y <double>: qy\n");
				printf ("-t <double>: kinetic constant\n");
				printf ("-e <double>: mu -- fermi energy\n");
				printf ("-f -- work at fixed beta_min\n");
				printf ("-w -- work at fixed wmega_min\n");
				printf ("-a <double>: alpha -- regularization parameter\n");
				printf ("-h -- print this list\n");
				exit (EXIT_SUCCESS);
			default:
				printf ("Unknown command!\n");
				exit (EXIT_FAILURE);
		}
	}

	calculate (M, fbeta, fwmega, qx, qy, 0.0, kin, a, "input.txt", "output.txt");

	if (!fbeta && !fwmega)
		system ("./slika.sh");

	else if (fbeta)
		system ("./slika1.sh");

	else if (fwmega)
		system ("./slika2.sh");

	exit (EXIT_SUCCESS);
}

