#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multimin.h>

typedef struct
{
			// dim	|	explanation
			// -----+-------------------------------------
	double * x,	// 3	| the x vectors for phase rotation
	       * y,	// 3	| the y vectors of eta
	       * u,	// 2	| u_1 and u_2 which represent vacuums v_1 and v_2
	       * d,	// 4 	| these are the yukawas
	       b,	// 1	| the angle beta
	       vd_bar,	// 1	| v_d = vd_bar \cos\beta, parametrized with gamma_d
	       vu_bar;	// 1	| v_u = vu_bar \sin\beta, parametrized with gamma_u

	// eigenvalues of each matrix -- 3 for each
	double * M,	// so this has dimension 9: first 3 or slect=1, and so on
	       * m,	// the correct masses + CKM matrix orthogonal angles
	       * e;	// errors of the masses + CKM matrix errors of the angles

	// CKM matrix parameters
	double * V,	// the CKM matrix, or Lu^T
	       * L;	// the L_d matrix
} hod;

typedef struct
{
	// the mass matrix and vectors
	double * a,
	       * b,
	       * M,
	       * d;
} mass;

// p is the parameter vector, containing calculated RG MSSM masses
// see arXiv:0712.1419v3 [hep-ph], 2008
// Table V, on page 20
void hod_init (hod * u, gsl_vector * v, double * p)
{
	u->x = (double *) malloc (3 * sizeof (double));
	u->y = (double *) malloc (3 * sizeof (double));
	u->u = (double *) malloc (2 * sizeof (double));
	u->d = (double *) malloc (4 * sizeof (double));
	u->M = (double *) malloc (12 * sizeof (double));
	u->m = (double *) malloc (12 * sizeof (double));
	u->e = (double *) malloc (12 * sizeof (double));
	u->V = (double *) malloc (9 * sizeof (double));
	u->L = (double *) malloc (9 * sizeof (double));

	int i;
	for (i = 0; i <= 2; i++)
		u->x[i] = gsl_vector_get (v, i);

	for (i = 0; i <= 2; i++)
		u->y[i] = gsl_vector_get (v, i+3);

	for (i = 0; i <= 1; i++)
		u->u[i] = gsl_vector_get (v, i+6);

	for (i = 0; i <= 3; i++)
		u->d[i] = gsl_vector_get (v, i+8);

	u->b		= gsl_vector_get (v, 12);

	// the two vacuums are both bounded between
	// 0 and 246/sqrt{2} GeV. To bound them to
	// that interval, we will use the funtion
	// 264/(pi*sqrt(2)) * (arctan gamma_{u,d} + pi/2)
	// which maps real axis to an open interval
	// (0, 246/sqrt(2))

	double gamma_d	= gsl_vector_get (v, 13),
	       gamma_u	= gsl_vector_get (v, 14);

//	these equations use -pi/2 < alpha_{u,d} < pi/2
//	u->vd_bar = (246.0/(M_PI * sqrt(2))) * (atan (gamma_d) + 0.5*M_PI);
//	u->vu_bar = (246.0/(M_PI * sqrt(2))) * (atan (gamma_u) + 0.5*M_PI);
//	we will use a full angle for alpha_{u,d}
	u->vd_bar = ((246.0 * sqrt(2))/M_PI) * atan (gamma_d);
	u->vu_bar = ((246.0 * sqrt(2))/M_PI) * atan (gamma_u);

	// here we have squares of masses
	for (i = 0; i <= 8; i++)
		u->m[i] = pow(p[i],2);

	// now we add the CKM parameters
	for (i = 9; i <= 11; i++)
		u->m[i] = p[i];

	// here we have errors for mass squares
	for (i = 0; i <= 8; i++)
		u->e[i] = 2 * p[i] * p[i+12];

	// and errors for the orthogonal angles
	for (i = 9; i <= 11; i++)
		u->e[i] = p[i+12];
}

void mass_init (mass * m)
{
	m->a = (double *) malloc (3 * sizeof (double));
	m->b = (double *) malloc (3 * sizeof (double));
	m->d = (double *) malloc (4 * sizeof (double));
	m->M = (double *) malloc (9 * sizeof (double));
}

void destroy_hod (hod * u)
{
	if (u->x) free (u->x);
	if (u->y) free (u->y);
	if (u->u) free (u->u);
	if (u->d) free (u->d);
	if (u->M) free (u->M);
	if (u->m) free (u->m);
	if (u->V) free (u->V);
	if (u->e) free (u->e);
	if (u->L) free (u->L);
	if (u) free (u);
}

void destroy_mass (mass * m)
{
	if (m->a) free (m->a);
	if (m->b) free (m->b);
	if (m->M) free (m->M);
	if (m->d) free (m->d);
	if (m) free (m);
}

// we assign vectors their values
void vec_gen (mass * m, hod * u, int select)
{
	// - EXPLANATION -
	// select = 0 => M_D: a = x_{d^c}, b = x_Q, d -> v_d * d = vd_bar * cos(beta) * d
	// select = 1 => M_U: a = x_{u^c}, b = x_Q, d -> v_u * d = vu_bar * sin(beta) * d
	// select = 2 => M_E: a = x_{e^c}, b = x_L, d -> (-3) * vd_bar * d = (-3) * v * cos(beta) * d

	int i;
	switch (select)
	{
		case 0:
			for (i = 0; i <= 2; i++)
			{
				m->a[i] = (u->x[i] + (u->u[0] - u->u[1])*u->y[i])/(1 + u->u[0] - u->u[1]);
				m->b[i] = (u->x[i] - u->u[0]*u->y[i])/(1 - u->u[0]);
				m->d[i] = u->d[i] * u->vd_bar * cos(u->b);
			}
			m->d[3] = u->d[3] * u->vd_bar * cos(u->b);
			break;
		case 1:
			for (i = 0; i <= 2; i++)
			{
				m->a[i] = (u->x[i] + (u->u[0] + u->u[1])*u->y[i])/(1 + u->u[0] + u->u[1]);
				m->b[i] = (u->x[i] - u->u[0]*u->y[i])/(1 - u->u[0]);
				m->d[i] = u->d[i] * u->vu_bar * sin(u->b);
			}
			m->d[3] = u->d[3] * u->vu_bar * sin(u->b);
			break;
		case 2:
			for (i = 0; i <= 2; i++)
			{
				m->a[i] = (u->x[i] - (3*u->u[0] + u->u[1])*u->y[i])/(1 - 3*u->u[0] - u->u[1]);
				m->b[i] = (u->x[i] + 3*u->u[0]*u->y[i])/(1 + 3*u->u[0]);
				m->d[i] = (-3)*u->d[i] * u->vd_bar * cos(u->b);
			}
			m->d[3] = u->d[3] * (-3)*u->vd_bar * cos(u->b);
			break;
		default:
			printf ("Bad selection option!\n");
			printf ("Aborting program.\n");
			exit (EXIT_FAILURE);
	}
}

// we construct the mass matrix -- non-diagonal
void matrix_gen (mass * m, hod * u, int select)
{
	double anorm	= m->a[0]*m->a[0] + m->a[1]*m->a[1] + m->a[2]*m->a[2],
	       bnorm	= m->b[0]*m->b[0] + m->b[1]*m->b[1] + m->b[2]*m->b[2],
	       dnorm	= m->a[0]*m->b[0]*m->d[0] + m->a[1]*m->b[1]*m->d[1] + m->a[2]*m->b[2]*m->d[2],
	       lambdaA	= 1.0/sqrt(1 + anorm),
	       lambdaB	= 1.0/sqrt(1 + bnorm),
	       alpha	= pow(lambdaA, 2)/(1 + lambdaA),
	       beta	= pow(lambdaB, 2)/(1 + lambdaB),
	       gamma	= beta*bnorm + alpha*anorm - m->d[3]*(1 + alpha*beta*anorm*bnorm) - alpha*beta*dnorm;

	int i,j;
	for (i = 0; i <= 2; i++)
	{
		// the general matrix without the diagonal d_{ij} part
		for (j = 0; j <= 2; j++)
		{
			m->M[j + 3*i] = (-1)* gamma * m->a[i] * m->b[j];
			m->M[j + 3*i] -= beta * m->b[j] * m->a[i] * m->d[i];
			m->M[j + 3*i] -= alpha* m->a[j] * m->d[j] * m->a[i];
		}

		// we also fix the diagonal parts
		m->M[0] += m->d[0];
		m->M[4] += m->d[1];
		m->M[8] += m->d[2];
	}
}

// thise function calculates M -> M^T M
void mass2 (mass * m)
{
	double * A = (double *) malloc (9 * sizeof(double));

	A[0] = m->M[0]*m->M[0] + m->M[3]*m->M[3] + m->M[6]*m->M[6];
	A[1] = m->M[0]*m->M[1] + m->M[3]*m->M[4] + m->M[6]*m->M[7];
	A[2] = m->M[0]*m->M[2] + m->M[3]*m->M[5] + m->M[6]*m->M[8];

	A[3] = m->M[1]*m->M[0] + m->M[4]*m->M[3] + m->M[7]*m->M[6];
	A[4] = m->M[1]*m->M[1] + m->M[4]*m->M[4] + m->M[7]*m->M[7];
	A[5] = m->M[1]*m->M[2] + m->M[4]*m->M[5] + m->M[7]*m->M[8];

	A[6] = m->M[2]*m->M[0] + m->M[5]*m->M[3] + m->M[8]*m->M[6];
	A[7] = m->M[2]*m->M[1] + m->M[5]*m->M[4] + m->M[8]*m->M[7];
	A[8] = m->M[2]*m->M[2] + m->M[5]*m->M[5] + m->M[8]*m->M[8];

	int i;
	for (i = 0; i <= 8; i++)
		m->M[i] = A[i];

	free (A);
}

// now we have to diagonalize the matrix
// can we assume that the matrix is positive definite?
void diagonalize (mass * m, hod * u, int select)
{
	int i, j;

	gsl_matrix_view mass_matrix = gsl_matrix_view_array (m->M, 3, 3);

	gsl_vector * eval = gsl_vector_alloc (3);
	gsl_matrix * evec = gsl_matrix_alloc (3, 3);
	
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
	gsl_eigen_symmv (&mass_matrix.matrix, eval, evec, w);
	gsl_eigen_symmv_free (w);
	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);

	// we use evec for our CKM matrix
	if (select == 0)
	{
		for (i = 0; i <= 2; i++)
		{
			for (j = 0; j <= 2; j++)
				u->L[j + i*3] = gsl_matrix_get (evec, j, i);
		}
	}

	if (select == 1)
	{
		double M [9];
		for (i = 0; i <= 2; i++)
		{
			for (j = 0; j <= 2; j++)
				M[j + i*3] = gsl_matrix_get (evec, i, j);
		}

		// now we can hard-code the CKM matrix
		u->V[0] = M[0]*u->L[0] + M[1]*u->L[3] + M[2]*u->L[6];
		u->V[1] = M[0]*u->L[1] + M[1]*u->L[4] + M[2]*u->L[7];
		u->V[2] = M[0]*u->L[2] + M[1]*u->L[5] + M[2]*u->L[8];

		u->V[3] = M[3]*u->L[0] + M[4]*u->L[3] + M[5]*u->L[6];
		u->V[4] = M[3]*u->L[1] + M[4]*u->L[4] + M[5]*u->L[7];
		u->V[5] = M[3]*u->L[2] + M[4]*u->L[5] + M[5]*u->L[8];

		u->V[6] = M[6]*u->L[0] + M[7]*u->L[3] + M[8]*u->L[6];
		u->V[7] = M[6]*u->L[1] + M[7]*u->L[4] + M[8]*u->L[7];
		u->V[8] = M[6]*u->L[2] + M[7]*u->L[5] + M[8]*u->L[8];

		// and calculate the tangeses orthogonal Euler angles
		u->M[9]	 = u->V[6]/u->V[3];		// \tan\theta_1
		u->M[10] = fabs(u->V[1]/u->V[0]);	// \tan\theta_2
		u->M[11] = (-1.0)*u->V[2]/u->V[1];	// \tan\theta_3
	}

	gsl_matrix_free (evec);

	for (i = 0; i <= 2; i++)
		u->M[i + 3*select] = gsl_vector_get (eval, i);

	gsl_vector_free (eval);
}

// we do this in all three cases for fermion masses
// and return the chi2
double diagonalize_and_chi (hod * u)
{

	mass * m = (mass *) malloc (sizeof (mass));
	mass_init (m);

	int i;
	for (i = 0; i <= 2; i++)
	{
		vec_gen (m, u, i);
		matrix_gen (m, u, i);
		mass2 (m);

		diagonalize (m, u, i);
	}
	destroy_mass (m);

	// the chi2 that comes from masses and other things
	double chi2 = 0;
	for (i = 0; i <= 11; i++)
		chi2 += pow (((u->M[i] - u->m[i])/u->e[i]), 2);

	chi2 /= 12.0;

	return chi2;
}

// we will use the minimizer from GSL without derivatives
double chi2func (const gsl_vector * v, void * p)
{
	hod * u = (hod *) malloc (sizeof (hod));
	hod_init (u, (gsl_vector *) v, (double *) p);

	double chi2 = diagonalize_and_chi (u);
	
	destroy_hod (u);
	
	return chi2;
}
/*
 *              *** printing functions ***
 *
 * These functions serve no other purpose, than to make results
 * more presentable to the eyes of the reader. More or less it's
 * just some switch statements.
 *
 *
 */

void matrix_print (double * M)
{
	int i, j;
	printf ("M\n");
	for (i = 0; i <= 2; i++)
	{
		for (j = 0; j <= 2; j++)
			printf("% 4.2lf", M[j + 3*i]);
		
		printf ("\n");
	}
}

void vectors_print (mass * m)
{
	int i;
	printf ("a\n");
	for (i = 0; i <= 2; i++)
		printf ("% 4.2lf ", m->a[i]);
	
	printf ("\nb\n");
	for (i = 0; i <= 2; i++)
		printf ("% 4.2lf ", m->b[i]);

	printf ("\nd\n");
	for (i = 0; i <= 3; i++)
		printf ("% 4.2lf ", m->d[i]);
	
	printf ("\n");
}

void print_fermion (FILE * fout, int select)
{
	switch (select)
	{
		case 0:
			fprintf (fout, "\td\n");
			break;
		case 1:
			fprintf (fout, "\ts\n");
			break;
		case 2:
			fprintf (fout, "\tb\n");
			break;
		case 3:
			fprintf (fout, "\tu\n");
			break;
		case 4:
			fprintf (fout, "\tc\n");
			break;
		case 5:
			fprintf (fout, "\tt\n");
			break;
		case 6:
			fprintf (fout, "\te\n");
			break;
		case 7:
			fprintf (fout, "\tmu\n");
			break;
		case 8:
			fprintf (fout, "\ttau\n");
			break;
		case 9:
			fprintf (fout, "\tt1\n");
			break;
		case 10:
			fprintf (fout, "\tt2\n");
			break;
		case 11:
			fprintf (fout, "\tt3\n");
			break;
		default:
			printf ("That's no charged fermion!\n");
			exit (EXIT_FAILURE);
	}
}

void parameter_print (FILE * fout, int select, gsl_vector * v)
{
	double gam = 0,
	       vac = 0;

	switch (select)
	{
		case 0:
			fprintf (fout, "% .15e", gsl_vector_get (v, 0));
			fprintf (fout, "\tx1\n");
			break;
		case 1:
			fprintf (fout, "% .15e", gsl_vector_get (v, 1));
			fprintf (fout, "\tx2\n");
			break;
		case 2:
			fprintf (fout, "% .15e", gsl_vector_get (v, 2));
			fprintf (fout, "\tx3\n");
			break;
		case 3:
			fprintf (fout, "% .15e", gsl_vector_get (v, 3));
			fprintf (fout, "\ty1\n");
			break;
		case 4:
			fprintf (fout, "% .15e", gsl_vector_get (v, 4));
			fprintf (fout, "\ty2\n");
			break;
		case 5:
			fprintf (fout, "% .15e", gsl_vector_get (v, 5));
			fprintf (fout, "\ty3\n");
			break;
		case 6:
			fprintf (fout, "% .15e", gsl_vector_get (v, 6));
			fprintf (fout, "\tu1\n");
			break;
		case 7:
			fprintf (fout, "% .15e", gsl_vector_get (v, 7));
			fprintf (fout, "\tu2\n");
			break;
		case 8:
			fprintf (fout, "% .15e", gsl_vector_get (v, 8));
			fprintf (fout, "\td1\n");
			break;
		case 9:
			fprintf (fout, "% .15e", gsl_vector_get (v, 9));
			fprintf (fout, "\td2\n");
			break;
		case 10:
			fprintf (fout, "% .15e", gsl_vector_get (v, 10));
			fprintf (fout, "\td3\n");
			break;
		case 11:
			fprintf (fout, "% .15e", gsl_vector_get (v, 11));
			fprintf (fout, "\td4\n");
			break;
		case 12:
			fprintf (fout, "% .15e", gsl_vector_get (v, 12));
			fprintf (fout, "\tbeta\n");
			break;
		case 13:
			gam = gsl_vector_get (v, 13);
			vac = ((246.0 * sqrt(2))/M_PI) * atan (gam);

			fprintf (fout, "% .15e\tgd", gam);
			fprintf (fout, "\t(% 4.4e\tvd_bar)\n", vac);
			break;
		case 14:
			gam = gsl_vector_get (v, 14);
			vac = ((246.0 * sqrt(2))/M_PI) * atan (gam);

			fprintf (fout, "% .15e\tgu", gam);
			fprintf (fout, "\t(% 4.4e\tvu_bar)\n", vac);
			break;
		default:
			printf ("There is no such parameter to print!\n");
			printf ("Error in parameter_print ()\n");
			exit (EXIT_FAILURE);
	}
}

void print_all (FILE * fout, const gsl_vector * v, double * p, double * a, int full)
{
	hod * u = (hod *) malloc (sizeof (hod));
	hod_init (u, (gsl_vector *) v, p);
	double chi2 = diagonalize_and_chi (u);

	// if we're outputting to a file, we only want
	// relavant solutions, that are small enough
	if (fout != stdout)
	{
		if (chi2 > 10)
		{
			fclose (fout);
			fout = fopen ("/dev/null", "a");
		}

		else
			printf ("\tHuzzah!\n");
	}

	// we're going to give full report
	if (full)
	{
		fprintf (fout, "==============================\n");
		fprintf (fout, "Chi^2 = % 4.4e\n", chi2);
		fprintf (fout, "------------------------------\n");
		fprintf (fout, "1.) Input parameters:\n\n");

		int i;
		for (i = 0; i <= 14; i++)
			fprintf (fout, "% lf ", a[i]);
		fprintf (fout, "\n");

		fprintf (fout, "------------------------------\n");
		fprintf (fout, "2.) Final parameters:\n\n");

		for (i = 0; i <= 14; i++)
			parameter_print (fout, i, (gsl_vector *) v);

		fprintf (fout, "\n");
		fprintf (fout, "------------------------------\n");
		fprintf (fout, "3.) Masses + CKM angles\n\n");
		
		for (i = 0; i <= 11; i++)
		{
			fprintf (fout, "M[%d] = % 4.4e \t m[%d] = % 4.4e",
					i, u->M[i], i, u->m[i]);
			print_fermion (fout, i);
		}

		fprintf (fout, "\n");
		fprintf (fout, "------------------------------\n");
		fprintf (fout, "END chi2 = % 4.4e\n", chi2);
	}

	// we only include the bare necessities to reproduce
	// the result
	else
	{
		int i;
		for (i = 0; i <= 14; i++)
			fprintf (fout, "% lf ", a[i]);
		fprintf (fout, "% 4.4e\n", chi2);
	}

	destroy_hod (u);
}

/*
 *              *** main function ***
 * Initial guess is accepted from a .txt file, which
 * can be specified manually, and passed as an argument
 * to the program, using the `-i some_file' or
 * `--input=some_file.txt' flag. It has to be formatted
 * in one row, which has parameters in this order:
 *
 * x_1 x_2 x_3 y_1 y_2 y_3 u_1 u_2 d_1 d_2 d_3 d_4 beta gu gd
 *
 * if there is no such file, it will be created with
 * default values, which can later on be modified, to
 * suit one's needs.
 *
 * Results can either be ouptuted to `stdout', or to
 * a file. In the latter case, cuts will be applied so
 * that only solutions with chi2 < 100 will be shown.
 *
 * For the rest, se and see the `-h', or `--help' options.
 *
 */

int main (int argc, char ** argv)
{

	//
	// start of the function is merely
	// bureaucracy and variable definitions/declarations
	//

	char * input	= NULL,
	     * output 	= NULL;

	int arg, c,
	    maxI	= 20000,
	    sloppy	= 0,
	    w		= 0,
	    f		= 0;

	double E_length	= 1e-5,		// final symplex step length
	       I_length	= 100.0,	// initial symplex step length
	       sl_err	= 0.05;		// sloppy relative error to use

	struct option longopts[] =
	{
		{"input-file",      required_argument,   NULL,   'i'},
		{"min-step",        required_argument,   NULL,   'M'},
		{"initial-step",    required_argument,   NULL,   'I'},
		{"max-iter",        required_argument,   NULL,   'n'},
		{"write",           optional_argument,   NULL,   'w'},
		{"write-full",      optional_argument,   NULL,   'f'},
		{"sloppy-errors",   optional_argument,   NULL,   's'},
		{"help",            no_argument,         NULL,   'h'},
		{"version",         no_argument,         NULL,   'v'},
		{NULL,              0,                   NULL,     0}
	};

	while ((arg = getopt_long (argc, argv, "i:M:I:n:w::f::s::hv", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 'i':
				input = (char *) malloc (30 * sizeof (char));
				strcpy (input, optarg);
				break;
			case 'M':
				E_length = atof (optarg);
				break;
			case 'I':
				I_length = atof (optarg);
				break;
			case 'n':
				maxI = atoi (optarg);
				break;
			case 'v':
				printf ("Charged fermion mass fitter v0.1\n");
				printf ("By Jozze Zobec, 2013\n");
				printf ("This program is licenced under WTFPLv2.0\n");
				exit (EXIT_SUCCESS);
			case 'w':
				w = 1;

				if (optarg)
				{
					output = (char *) malloc (30 * sizeof (char));
					strcpy (output, optarg);
				}

				break;
			case 'f':
				f = 1;
				w = 1;

				if (optarg)
				{
					output = (char *) malloc (30 * sizeof (char));
					strcpy (output, optarg);
				}

				break;
			case 's':
				sloppy = 1;

				if (optarg)
					sl_err = atof (optarg);

				break;
			case 'h':
				printf ("List of commands\n");
				printf ("------------------------------\n");
				printf ("-i <file>\n");
				printf ("--input-file=<file>\n");
				printf ("    Program reads its initial parameters from a configurable file.\n");
				printf ("    If it doesn't exist, it  will be  automatically created by the\n");
				printf ("    program, using default values. The hash `#' symbol can be used\n");
				printf ("    as a comment sign. File must contain one line with 15 entries,\n");
				printf ("    which mean the following: x_1, x_2, x_3,  y_1, y_2, y_3,  u_1,\n");
				printf ("    u_2, d_1, d_2, d_3, d_4, beta, gu, gd.\n\n");
				printf ("-M <real number>\n");
				printf ("--min-step=<real number>\n");
				printf ("    This number represents which is the minimal length for our sym-\n");
				printf ("    plex algorithm to stop searching for a minimum. In other words,\n");
				printf ("    it's a stopping condition for our program.\n\n");
				printf ("-I <real number>\n");
				printf ("--initial-step=<real number>\n");
				printf ("    Initial step length -- if our initial guess was  too close to a \n");
				printf ("    local minimum, we can remedy this by using a larger initial\n");
				printf ("    step.\n\n");
				printf ("-n <integer>\n");
				printf ("--max-iter=<integer>\n");
				printf ("    Another stopping condition is maximum number of iterations. It's\n");
				printf ("    possible  that minimizer  falls into a  region where the mininum\n");
				printf ("    step length condition will never be met, so we bound the  number\n");
				printf ("    of iterations.\n\n");
				printf ("-w\n");
				printf ("--write\n");
				printf ("    Don't output files to `stdout', but to a file, called `runner-o-\n");
				printf ("    utput.txt'. The program will output only the bare neccesities.\n\n");
				printf ("-f [output-file.txt]\n");
				printf ("--write-full[=output-file.txt]\n");
				printf ("    This flag  also sets the output  to a file, but  will output the\n");
				printf ("    results in a full format, which will print the final  results as\n");
				printf ("    if they were taken from directly from the `stdout'.\n\n");
				printf ("-s [integer]\n");
				printf ("--sloppy-errors[=integer]\n");
				printf ("    Don't use the errors from the arXiv:0712.1419v3 [hep-ph] report,\n");
				printf ("    but use 5%% relative error for each value, whose  error is lower\n");
				printf ("    than 5%%.\n\n");
				printf ("-h\n");
				printf ("--help\n");
				printf ("    Prints this list.\n\n");
				printf ("-v\n");
				printf ("--version\n");
				printf ("    Prints program version information.\n");
				exit (EXIT_SUCCESS);
			default:
				printf ("Unknown command!\n");
				printf ("Wrong program usage!\n");
				printf ("See `%s -h' for instructions!\n", argv[0]);
				exit (EXIT_FAILURE);
		}
	}

	if (input == NULL)
	{
		input = (char *) malloc (9 * sizeof (char));
		sprintf (input, "input.txt");
	}

	FILE * fin = fopen (input, "r");

	// we will create one
	if (fin == NULL)
	{
		printf ("Creating the default config file...\t");
		FILE * fdefault = fopen ("input.txt", "w");
		fprintf (fdefault, "# This file has been generated automatically by the %s program and can be\n",
				argv[0]);
		fprintf (fdefault, "# modified at any time, to change the initial parameters. The  hash sign,\n");
		fprintf (fdefault, "# `#', can be used  to comment out lines, or  portions to the end of  the\n");
		fprintf (fdefault, "# line. In this particular case, use `-I2.1' options to see the best results\n");
		fprintf (fdefault, "# These values should give chi^2 = 1.0004e+05, which is pretty good.\n");
		fprintf (fdefault, "######################################################################################\n");
		fprintf (fdefault, "#x1  |  x2 |  x3 |  y1 |  y2 |  y3 |  u1 |  u2 |  d1 |  d2 |  d3 |  d4 |  b |  gd | gu\n");
		fprintf (fdefault, "100    2     3.6    1     2      0   180     2   0.1     1    10     1    2  1000 -0.2\n");
		fclose (fdefault);
		printf ("Done.\n");

		fin = fopen ("input.txt", "r");
	}

	// now we read the input parameters from the file
	double * a = (double *) malloc (15 * sizeof (double));
	do
	{
		c = fgetc (fin);
		if (c == '#')
		{
			do
			{
				c = fgetc (fin);
			} while (c != '\n' || c == '#');
		}

		else
		{
			fseek (fin, -1, SEEK_CUR);
			fscanf (fin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
					&a[0], &a[1], &a[2],		// parameters of the vector `x'
					&a[3], &a[4], &a[5],		// parameters of the vector `y'
					&a[6], &a[7],			// parameters of the vector `u'
					&a[8], &a[9], &a[10], &a[11],	// parameters of the vector `d'
					&a[12],				// the angle beta
					&a[13], &a[14]);		// vacuum angles gd and gu
			fclose (fin);
		}
	} while (c != EOF);

	gsl_vector_const_view v = gsl_vector_const_view_array (a, 15);

	// Known (masses & co.) + (CKM angles) go here. Relevant units are GeV
	double Masses [24] = { 0.70e-3,		// d
			      13.0e-3,		// s
			      0.79,		// b
		              0.49e-3,		// u
			      0.236,		// c
			      92.0,		// t
			      0.283755495e-3,	// e
			      59.9033617e-3,	// mu
			      1021.95e-3,	// tau
			      -3.587e-3,	// \theta_1
			      0.2316,		// \theta_2
			      5.5979,		// \theta_3
			      0.30e-3,		// d (error)
			      4e-3,		// s (error)
			      0.04,		// b (error)
			      0.17e-3,		// u (error)
			      0.036,		// c (error)
			      7.8,		// t (error)
			      24e-12,		// e (error)
			      54e-9,		// mu (error)
			      0.11e-3,		// tau (error)
			      0.928e-3,		// \theta_1 (error)
			      0.0001849,	// \theta_2 (error)
			      1.5899		// \theta_3 (error)
		};

	// we use the sloppy error values
	if (sloppy == 1)
	{
		for (c = 0; c <= 11; c++)
		{
			double new_error = sl_err * Masses [c];

			if (Masses[c + 12] < new_error)
				Masses [c + 12] = new_error;
		}
	}

	
	//
	// the actual minimizer code
	//

	const gsl_multimin_fminimizer_type * T =
		gsl_multimin_fminimizer_nmsimplex2;

	gsl_vector * step_s = gsl_vector_alloc (15);
	gsl_vector_set_all (step_s, I_length);

	gsl_multimin_function minex_func;
	minex_func.f		= &chi2func;
	minex_func.n 		= 15;
	minex_func.params	= Masses; 

	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 15);
	gsl_multimin_fminimizer_set (s, &minex_func, &v.vector, step_s);

	int iter = 0,
	    status;

	double size;

	do
	{
//		printf ("Iter =% 3d\n", iter);
		iter++;
		status = gsl_multimin_fminimizer_iterate (s);

		if (status)
			break;

		size	= gsl_multimin_fminimizer_size (s);
		status	= gsl_multimin_test_size (size, E_length);

//		if (status == GSL_SUCCESS)
//			printf ("Converged to minimum in %d iterations!\n", iter);
	} while (status == GSL_CONTINUE && iter < maxI);

	if (w)
	{
		if (output == NULL)
		{
			output = (char *) malloc (30 * sizeof(char));
			sprintf (output, "output.txt");
		}

		FILE * fout = fopen (output, "a");
		print_all (fout, s->x, Masses, a, f);

		fclose (fout);
	}

	else
		print_all (stdout, s->x, Masses, a, 1);

	if (s) gsl_multimin_fminimizer_free (s);
	if (output) free (output);
	gsl_vector_free (step_s);
	free (input);
	free (a);

	return 0;
}

