/*
 *                  A FOREWORD ON THE RUNNER PROGRAM
 *
 * This program is given a config file, much like the one fitter.c requires.
 * The difference here is, that it takes two lines. First line is for lower
 * boundary of the parameter, and the 2nd line is for upper boundary. The
 * program then takes 50 000 random vectors within that portion of phase spa-
 * ce, although that number can be changed. Program will then create a file
 * with results for each of these vectors: displaying initial parameters and
 * chi^2 after minimization. It will only include those vectors, whose chi^2
 * is lower than a certain value, which can be changed dynamically.
 *
 */


#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>

// we read the date from the file
void reader (FILE * fin, double * bound)
{
	int pass = 0,
	    c;

	// we want to use `#' as comment markers
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

		else if (pass < 2)
		{
			fseek (fin, -1, SEEK_CUR);
			fscanf (fin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
					&bound[0 + pass*15], &bound[1 + pass*15], &bound[2 + pass*15],
					&bound[3 + pass*15], &bound[4 + pass*15], &bound[5 + pass*15],
					&bound[6 + pass*15], &bound[7 + pass*15], &bound[8 + pass*15],
					&bound[9 + pass*15], &bound[10+ pass*15], &bound[11+ pass*15],
					&bound[12+ pass*15], &bound[13+ pass*15], &bound[14+ pass*15]);
			pass++;
		}

		// we only read two lines
		if (pass == 2 || c == EOF)
			break;
	}
}

void randvec (gsl_rng * rand, double * param, double * bound)
{
	int i;
	for (i = 0; i <= 14; i++)
	{
		double uniform = gsl_rng_uniform (rand);
		param[i] = bound[i] + uniform*(bound[i+15] - bound[i]);
	}
}

void fitter_minimize (double * param, FILE * finit, char * command)
{
	int i;

	// we create input for fitter.c
	for (i = 0; i <= 14; i++)
		fprintf (finit, "%lf ", param[i]);

	rewind (finit);

	system (command);
}

int main (int argc, char ** argv)
{
	int N		= 50000000,	// number of random points we generate
	    maxI	= 2000000,	// maximum iterations for fitter.c
	    seed	= 0,		// Marsenne Twister seed
	    sloppy	= 0,		// sloppyness settings
	    full	= 0,		// output format
	    arg, i;

	char * input	= NULL,
	     * output	= NULL;
	
	double sl_error		= 0.05,
	       min_length	= 0.0001,
	       init_length	= 3.0;

	// we declare our options
	struct option longopts[] =
	{
		{"input-file",        required_argument,     NULL,    'i'},
		{"no-random-points",  required_argument,     NULL,    'N'},
		{"min-step",          required_argument,     NULL,    'M'},
		{"initial-step",      required_argument,     NULL,    'I'},
		{"max-iter",          required_argument,     NULL,    'n'},
		{"output-file",       required_argument,     NULL,    'o'},
		{"full",              no_argument,           NULL,    'f'},
		{"sloppy-errors",     optional_argument,     NULL,    's'},
		{"help",              no_argument,           NULL,    'h'},
		{"version",           no_argument,           NULL,    'v'},
		{NULL,                0,                     NULL,      0}
	};

	while ((arg = getopt_long (argc, argv, "i:N:M:I:n:o:fs::hv", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 'i':
				input = (char *) malloc (30 * sizeof (char));
				strcpy (input, optarg);
				break;
			case 'N':
				N = atoi (optarg);
				break;
			case 'M':
				min_length = atof (optarg);
				break;
			case 'I':
				init_length = atof (optarg);
				break;
			case 'n':
				maxI = atoi (optarg);
				break;
			case 'o':
				output = (char *) malloc (30 * sizeof (char));
				strcpy (output, optarg);
				break;
			case 'f':
				full = 1;
				break;
			case 's':
				sloppy = 1;
				if (optarg)
					sl_error = atof (optarg);

				break;
			case 'h':
				printf ("List of commands:\n");
				printf ("---------------------------\n");
				printf ("-i <filename>\n");
				printf ("--input-file=<filename>\n");
				printf ("    Program reads boundaries for each parameter from a configuration\n");
				printf ("    file. It needs to be  formatted, like the file for  the fitter.c\n");
				printf ("    with the exception, that this time there are two lines. The  fi-\n");
				printf ("    rst line is lower boundary,  and the 2nd is  the upper  boundary\n");
				printf ("    for the specific variable in the column\n\n");
				printf ("-N <integer>\n");
				printf ("--num-random-points=<integer>\n");
				printf ("    Number of random points (vectors)  we want our program to select\n");
				printf ("    for the initial condition. Each  of those points is selected  so\n");
				printf ("    that is within the aforementioned boundaries.\n\n");
				printf ("-M <real number>\n");
				printf ("--min-step=<real number>\n");
				printf ("    This is the stopping parameter to our fitter.c program.\n\n");
				printf ("-I <real number>\n");
				printf ("--initial-step=<real number>\n");
				printf ("    The initial step length to our fitter.c program\n\n");
				printf ("-n <integer>\n");
				printf ("--max-iter=<integer>\n");
				printf ("    Another stepping  parameter for the fitter.c program, is the nu-\n");
				printf ("    mber of maximum allowed interations.\n\n");
				printf ("-o <filename>\n");
				printf ("--output-file=<filename>\n");
				printf ("    The name of the output file for the fitter.c program.\n\n");
				printf ("-f\n");
				printf ("--full\n");
				printf ("    Enable the more detailed format for the fitter.c program output.\n");
				printf ("\n-s[=real number]\n");
				printf ("--sloppy-errors[=real number]\n");
				printf ("    The  minimum allowed relative error for the fitter.c parameters.\n");
				printf ("\n-h\n");
				printf ("--help\n");
				printf ("    Prints this list.\n");
				printf ("-v\n");
				printf ("--version\n");
				printf ("    Prints version details.\n");
				exit (EXIT_SUCCESS);
			case 'v':
				printf ("Yeah, this is my program. Got a problem with that?\n");
				exit (EXIT_SUCCESS);
			default:
				fprintf (stderr, "Sorry, but this is an invalid option.\nFsck you too!\n");
				exit (EXIT_FAILURE);
		}
	}

	if (!input)
	{
		input = (char *) malloc (30 * sizeof (char));
		sprintf (input, "runner-in.txt");
	}

	if (!output)
	{
		output = (char *) malloc (30 * sizeof (char));
		sprintf (output, "output-runner.txt");
	}

	FILE * fin = fopen (input, "r"),
	     * ftmp = fopen ("runner4fitter.txt", "w"),
	     * urand = fopen ("/dev/urandom", "r");
	
	if (!fin)
	{
		printf ("Error! No input file!\nUsing the default one.\n");
		
		FILE * fdefault = fopen (input, "w");
		fprintf (fdefault, "# This file has been automatically generated by the %s program, and can\n", argv[0]);
		fprintf (fdefault, "# modified at any time. First line is the lower boundary, lower line is\n");
		fprintf (fdefault, "# the upper boundary. The hash sign, `#', can be used as a comment mar-\n");
		fprintf (fdefault, "# ker. For ease of use, the meaning of each line has been marked with a\n");
		fprintf (fdefault, "# correspinding variable symbol.\n");
		fprintf (fdefault, "######################################################################################\n");
		fprintf (fdefault, "#x1  |  x2 |  x3 |  y1 |  y2 |  y3 |  u1 |  u2 |  d1 |  d2 |  d3 |  d4 |  b |  gd | gu\n");
		fprintf (fdefault, "-32    +50    0   -177   -220  -0.1  -237   -75  -12   -12   +50  -1010  10   -660 +65\n");
		fprintf (fdefault, "-5     +70   +40  -150   -160  +0.1  -190  +105   +8    +8    +75   -990  32  -610 +87\n");
		fprintf (fdefault, "\n");

		fclose (fdefault);
		fin = fopen (input, "r");
	}

	if (!output)
	{
		printf ("Error! No output designated!\n");
		printf ("Yes, not even the default one ...\nExiting.\n");
		exit (EXIT_FAILURE);
	}
	
	double * p = (double *) malloc (30 * sizeof (double)),
	       * v = (double *) malloc (15 * sizeof (double));
	
	char * command = (char *) malloc (100 * sizeof (char));
	sprintf (command, "./fermion -i runner4fitter.txt -M%lf -n%d -f%s -I%lf",
			min_length, maxI, output, init_length);

	if (sloppy)
		sprintf (command, "%s -s%lf", command, sl_error);

	gsl_rng * rand = gsl_rng_alloc (gsl_rng_mt19937);
	fread (&seed, sizeof (int), 1, urand);
	gsl_rng_set (rand, seed);
	fclose (urand);

	reader (fin, p);
	fclose (fin);

	for (i = 0; i <= N-1; i++)
	{
		randvec (rand, v, p);
		fitter_minimize (v, ftmp, command);
		printf ("n = % d\n", i);
	}

	gsl_rng_free (rand);
	fclose (ftmp);
	free (p);
	free (v);

	exit(EXIT_SUCCESS);
}

