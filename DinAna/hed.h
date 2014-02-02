// hed.h
///////////
#ifndef DINANA_HED
#define DINANA_HED

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int
compare (const void * p1, const void * p2)
{
	double left = *(const double *) p1,
	       right= *(const double *) p2;

	return ((left > right) - (left < right));
}

double mod (double a, double b)
{
	// basically "a mod b", but for doubles
	double n = 0,
	       r = 0;

	if (b == 0)
		r = b;
	else
		r = b * modf (a/b, &n);

	if (r < 0)
		r += b;
	return r;
}

typedef struct
{
	int G;		// number of random trajectories
	size_t T,	// length of the array 't'
	       J;	// length of the Lyapunov exponent array

	double k1,	// k1 of the 1st kicker	
	       k2,	// k2 of the 2nd kicker
	       tau,	// frequency ratio
	       top,	// for lyapunov exponents
	       theta,	// phase-delay between k1 and k2 kicks
	       tmax;	// maximum time of observation

	double * L,	// progress of Lyapunov exponent
	       * C,	// time correlations -- diffusion parameter
	       * D,	// Delta_i -- length of time intervals between kicks
	       * K,	// strength of each kick (at a time t_i)
	       * t;	// array of "times" of kicks
	
	double ** p,	// array for momenta of all the trajectories at all times
	       ** x;	// array for all positions of all trajectories at all times

} ph;

void
getTimes (ph * u)
{
	double tmax = u->tmax,
	       tau  = u->tau,
	       theta= u->theta;

	// if it's uninitialized
	if (!u->t)
		u->t = (double *) malloc (sizeof(double));
	double * cas = (double *) malloc (sizeof(double));

	// now we fill it with all the integers
	int i;
	for (i = 0; i <= (int) u->tmax; i++)
	{
		// first we enlarge it
		cas = (double *) realloc (cas, (i+1)*sizeof(double));
		// now we write in it
		cas [i] = (double) i;
	}

	// we need to fill it with all the multiples of tau
	int T = 0;
	size_t N = (size_t) tmax + 1;
	while (T*tau - theta <= tmax)
	{
		// again we enlarge it
		cas = (double *) realloc (cas, (N + T + 1)*sizeof(double));

		// and again write in it
		cas [N+T] = T*tau + theta;
		T++;
	}

	// now we gotta sort it
	N = N + T;
	qsort (cas, N, sizeof(double), &compare);

	// and now we have to filter out recurring instances and write them out
	// in our output time array
	T = 0;
	for (i = 0; i <= N-1; i++)
	{
		// if we don't have double occurance ...
		if (cas[i+1] != cas[i])
		{
			// we enlarge t to accomodate a new entry ...
			u->t = (double *) realloc (u->t, (T+1)*sizeof(double));

			// and write it in
			u->t[T] = cas[i];
			T++;
			u->T = T;
		}
	}

	// now we remove the auxilliary array and voila!
	if (cas) free (cas);
}

void
Deltas (ph * u)
{
	/* D, K and t don't have to be initialized */
	if (!u->D) u->D = (double *) malloc (sizeof (double));
	if (!u->K) u->K = (double *) malloc (sizeof (double));
	if (!u->t) u->t = (double *) malloc (sizeof (double));

	// now we need to fill up the `t' array
	// first with all the integer 't'-s
	getTimes (u);

	// we enlarge the arrays
	int N = u->T - 1,
	    i;
	
	u->K = (double *) realloc (u->K, N * sizeof (double));
	u->D = (double *) realloc (u->D, N * sizeof (double));

	// and now we put some crap in them
	for (i = 0; i <= N-1; i++)
	{
		// this one is simple
		u->D[i] = u->t[i+1] - u->t[i];

		// this one has some if clauses ...
		if (mod (u->t[i], 1) == 0)
		{
			if (mod(u->t[i], u->tau) - u->theta == 0)
				u->K[i] = u->k1 + u->k2;
			else
				u->K[i] = u->k1;
		}
		else
			u->K[i] = u->k2;
	}
}

void check (ph * u)
{
	size_t i;
	for (i = 0; i <= u->T-2; i++)
		printf("%.2lf\t %.2e\t %.2e\n",
				u->t[i], u->D[i], u->K[i]);
	printf ("%.2lf\n", u->t[u->T-1]);
}

void initPh (ph * u,
		int G,
		double tmax,
		double top,
		double tau,
		double theta,
		double k1,
		double k2)
{
	// first we initialize u
	if (!u) u = (ph *) malloc (sizeof (ph));

	u->G	= G;
	u->tmax	= tmax;
	u->top	= top;
	u->tau 	= tau;
	u->theta= theta;
	u->k1	= k1;
	u->k2	= k2;

	// we initialize the "Deltas"
	Deltas (u);

//	check (u);

	size_t N = u->T;
	u->L = (double *) malloc (sizeof(double));
	u->C = (double *) malloc (N * sizeof(double));

	u->p = (double **) malloc (u->G * sizeof (double *));
	u->x = (double **) malloc (u->G * sizeof (double *));

	int i;
	for (i = 0; i <= G-1; i++)
	{
		u->p[i] = (double *) malloc (N * sizeof (double));
		u->x[i] = (double *) malloc (N * sizeof (double));
	}
}

void freePh (ph * u)
{
	if (u->L) free (u->L);
	if (u->C) free (u->C);
	if (u->D) free (u->D);
	if (u->K) free (u->K);
	if (u->t) free (u->t);

	int i;
	for (i = 0; i <= u->G-1; i++)
	{
		if (u->p[i]) free (u->p[i]);
		if (u->x[i]) free (u->x[i]);
	}

	if (u->p) free (u->p);
	if (u->x) free (u->x);

	if (u) free (u);
}

// propagate
void phaseD (ph * u, int k)
{
	int N = u->T -1,
	    i;

	for (i = 0;i <= N-1; i++)
	{
		u->p[k][i+1] = u->p[k][i] + u->K[i]*sin(u->x[k][i]);
		u->x[k][i+1] = mod(u->x[k][i] + u->D[i]*u->p[k][i+1], 2*M_PI);
	}
}

// all propagate for all vectors
void phases (ph * u)
{
	// initialize the seed for random generator
	srandom (time(NULL) % 23500);

	int k;
	for (k = 0; k <= u->G-1; k++)
	{
		// first we set the initial conditions
		u->x[k][0] = 2*M_PI * ((double) random())/RAND_MAX;
		u->p[k][0] = 12 * ((double) random())/RAND_MAX - 6;

		// now we propagate
		phaseD (u, k);
	}
}

// trajectories
void trajectories (ph * u)
{
	char dat [40];
	sprintf (dat, "t%.2lf-q%.2lf-k%.2lf-K%.2lf.txt",
			u->tau, u->theta, u->k1, u->k2);
	
	FILE * fout = fopen (dat, "w");

	size_t k, i;
	for (k = 0; k <= u->G - 1; k++)
	{
		for (i = 0; i <= u->T - 1; i++)
			fprintf (fout, "% .12e\t% .12e\t 0.01\n",
					u->x[k][i], u->p[k][i]);
	}

	fclose (fout);
}

// in case we need to separate trajectories ...
void sep_trajectories (ph * u)
{
	char dat1 [40],
	     dat2 [40];

	sprintf (dat1, "t%.2lf-q%.2lf-k%.2lf-K%.2lf-1.txt",
			u->tau, u->theta, u->k1, u->k2);
	sprintf (dat2, "t%.2lf-q%.2lf-k%.2lf-K%.2lf-2.txt",
			u->tau, u->theta, u->k1, u->k2);

	FILE * fout1 = fopen (dat1, "w"),
	     * fout2 = fopen (dat2, "w");

	size_t k, i;
	for (k = 0; k <= u->G-1; k++)
	{
		for (i = 0; i <= u->T-1; i++)
		{
			// if it's from the 2nd kick, put it in this plot
			if ((u->K[i] == u->k2 || u->K[i] == u->k1 + u->k2) && mod (u->t[i], u->tau) - u->theta == 0)
				fprintf (fout2, "% .12e\t% .12e\t 0.01\n",
						u->x[k][i], u->p[k][i]);
			// if it's from the 1st kick, put it here ;)
			if ((u->K[i] == u->k1 && mod(u->t[i],1) == 0))
				fprintf (fout1, "% .12e\t% .12e\t 0.01\n",
						u->x[k][i], u->p[k][i]);
		}
	}
}

// plot those trajectories with a shell script
void plot (ph * u)
{
	char call [50];
	if (u->theta != 0 || u->tau != 1)
		sprintf(call, "./plot2.sh %.2lf %.2lf %.2lf %.2lf",
				u->tau, u->theta, u->k1, u->k2);
	else
		sprintf(call, "./plot.sh %.2lf %.2lf %.2lf %.2lf",
				u->tau, u->theta, u->k1, u->k2);
	system (call);
}

// correlations
void corr (ph * u)
{
	char dat [40];
	sprintf (dat, "corr-%.2lf-%.2lf-%.2lf-%.2lf.txt",
			u->tau, u->theta, u->k1, u->k2);

	FILE * fout = fopen (dat, "w");
	int N = u->T,
	    i, k;

	// u->C is actually used ONLY for diffusion constant
	// correlations will not be stored in struct ph * u,
	// but written directly to stdout
	double Cx [N],
	       Cp [N],
	       p0 = 0,
	       pt = 0,
	       x0 = 0,
	       xt = 0,
	       nx = 0,
	       np = 0;

	u->C[0] = 0;

	// now let's get down to business
	for (i = 0; i <= N-1; i++)
	{
		u->C[i] = 0;
		Cx[i]   = 0;
		Cp[i]   = 0;
		pt      = 0;
		xt	= 0;
		for (k = 0; k <= u->G-1; k++)
		{
			u->C[i] += pow(u->p[k][i] - u->p[k][0], 2);
			Cx [i] += u->x[k][i] * u->x[k][0];
			Cp [i] += u->p[k][i] * u->p[k][0];
			pt += u->p[k][i]/u->G;
			xt += u->x[k][i]/u->G;

			if (i == 0)
			{
				p0 += u->p[k][i]/u->G;
				x0 += u->x[k][i]/u->G;
			}
		}
		u->C[i] /= u->G;

		Cx[i] /= u->G;
		Cp[i] /= u->G;
		Cx[i] -= x0*xt;
		Cp[i] -= p0*pt;

		if (i == 0)
		{
			np = Cp[i];
			nx = Cx[i];
		}

		Cx[i] /= nx;
		Cp[i] /= np;

		// now we print them
		fprintf (fout, "%lf\t% .12e\t% .12e\t% .12e\n",
				u->t[i], u->C[i], Cx[i], Cp[i]);
	}

	fclose (fout);
}

// lyapunov exponents
double ljapunov (ph * u)
{
	// starting difference
	int G = u->G,
	    j = 1,
	    k, i;
	
	double lam = 0,
	       d0  = 1e-6,
	       dt  = 0;

	double * cas = (double *) malloc (sizeof (double));

	if (u->G % 2 != 0)
		G = u->G - 1;

	srandom (time(NULL) % 23500);
	for (k = 0; k <= G-2; k += 2)
	{
		u->x[k][0] = 2*M_PI * ((double) random ())/RAND_MAX;
		u->p[k][0] = 12 * ((double) random ())/RAND_MAX - 6;

		u->x[k+1][0] = u->x[k][0] + d0;
		u->p[k+1][0] = u->p[k][0];

		j = 1;
		lam = 0;
		for (i = 0; i <= u->T-2; i++)
		{
			u->p[k][i+1] = u->p[k][i] + u->K[i]*sin(u->x[k][i]);
			u->x[k][i+1] = mod (u->x[k][i] + u->D[i]*u->p[k][i+1], 2*M_PI);

			u->p[k+1][i+1] = u->p[k+1][i] + u->K[i]*sin(u->x[k+1][i]);
			u->x[k+1][i+1] = mod (u->x[k+1][i] + u->D[i]*u->p[k+1][i+1], 2*M_PI);

			dt = fabs (u->x[k][i] - u->x[k+1][i]);
			if (mod(u->t[i], u->top) == 0)
			{
				lam += log (dt/d0);

				if (k == 0)
				{
					u->L = (double *) realloc (u->L, j * sizeof (double));
					cas  = (double *) realloc (cas, j * sizeof (double));
					u->J = j;
					u->L[j-1] = 0;
					cas[j-1] = u->t[i];
				}

				u->x[k+1][i+1] = u->x[k][i+1] + d0;	// we rescale the orbit
				u->p[k+1][i+1] = u->p[k+1][i+1];

				u->L[j-1] += lam/G;
				j++;
			}
		}
	}

	// we fit a straight line to the data
	double g = 0,
	       h = 0,
	       a = 0,
	       b = 0,
	       n = 0,
	       C = 0;
	for (j = 0; j <= u->J-1; j++)
	{
		if (cas[j] >= 400 && cas[j-1] < 400)
			n = u->J - j;
		if (cas[j] > 400)
		{
			g += cas[j];
			h += cas[j] * cas[j];
			a += u->L[j];
			b += cas[j] * u->L[j];
		}
	}

	lam = (n*b - g*a)/(n*h - g*g);	// the ljapunov exponent
	C   = (h*a - g*b)/(n*h - g*g);		// the constant shift

	// now that we have the average Lyapunov stuff, let's print them out in full glory
	FILE * fout = fopen ("ljapun.txt", "w");
	for (j = 0; j <= u->J-1; j++)
		fprintf (fout, "% .12e\t % .12e\n",
				cas[j], u->L[j]);

	fclose (fout);
	free (cas);

	printf ("lam = % .12e\tC = % .12e\n",
			lam, C);

	return lam;
}

void ljpScan (ph * u)
{
	double dk   = 0.025,
	       kmin = 0.1,
	       K    = 1,
	       L    = 0;

	u->k1	= kmin;
	u->k2	= kmin;
	u->top	= 1;

//	Has to be 
//	u->tau = 0.414213562373095
//	in order to make sense ...
	FILE * fout = fopen ("ljp-scan.txt", "w");
	while (u->k1 <= K)
	{
		u->k2 = kmin;
		while (u->k2 <= K)
		{
			L = ljapunov (u);
			fprintf (fout, "%.12e\t %.12e\t %.12e\n",
					u->k1, u->k2, L);

			u->k2 += dk;
		}
		fprintf (fout, "\n");
		u->k1 += dk;
	}
	fclose (fout);
}

#endif
