// prg.c
///////////
#include "hed.h"
#include <getopt.h>
#include <string.h>

int main (int argc, char ** argv)
{
	// default values
	int G = 600,
	    arg;

	double tmax  = 2000,
	       tau   = 1,
	       top   = 30,
	       theta = 0,
	       k1    = 0.3,
	       k2    = 0.6;

	// some interface so we don't have to recompile
	struct option longopts [] =
	{
		{ "tmax",	required_argument,	NULL,	'T' },
		{ "tau",	required_argument,	NULL,	't' },
		{ "k1",		required_argument,	NULL,	'k' },
		{ "k2",		required_argument,	NULL,	'K' },
		{ "theta",	required_argument,	NULL,	'h' },
		{ "top",	required_argument,	NULL,	'D' },
		{ "Gr",		required_argument,	NULL,	'G' },
	};

	while ((arg = getopt_long (argc, argv, "T:t:k:K:h:G:", longopts, NULL)) != -1)
	{
		switch (arg)
		{
			case 'T':
				tmax = atof (optarg);
				break;
			case 't':
				tau = atof (optarg);
				break;
			case 'k':
				k1 = atof (optarg);
				break;
			case 'K':
				k2 = atof (optarg);
				break;
			case 'h':
				theta = atof (optarg);
				break;
			case 'D':
				top = atof (optarg);
				break;
			case 'G':
				G = atoi (optarg);
				break;
			default:
				printf ("Error!\nUnknown command!\n");
				exit (EXIT_FAILURE);
		}
	}

	// we initialize our phases struct
	ph * u = (ph *) malloc (sizeof(ph));
	initPh (u, G, tmax, top, tau, theta, k1, k2);
//	phases (u);

	// we put them in a file
//	trajectories (u);
//	if (theta != 0 || tau != 1)
//		sep_trajectories (u);
//	// which we now plot
//	plot (u);
//	corr (u);
//	ljapunov (u);
	ljpScan (u);

	freePh (u);

	exit (EXIT_SUCCESS);
}
