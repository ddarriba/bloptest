/*
 ============================================================================
 Name        : pll-jctest.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include "bloptest.h"

#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>
#include <nlopt.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>

unsigned num_branches = 0;
double time_limit = 1.0;
double epsilon = 0.001;

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data) {
	tree_data * trd = (tree_data *) my_func_data;

#ifdef PRINT_BRANCHES
	printf("Br: ");
#endif
	int i;
	for (i = 0; i < num_branches; i++) {
		trd->branches[i]->z[0] = x[i];
		trd->branches[i]->back->z[0] = x[i];
#ifdef PRINT_BRANCHES
		printf("%8.4lf ", pllGetBranchLength(trd->tr, trd->branches[i], 0));
#endif
	}

	pllEvaluateLikelihood(trd->tr, trd->pr, trd->tr->start, true, false);
	double lk = trd->tr->likelihood;

#ifdef PRINT_BRANCHES
	printf(" (%lf)\n", lk);
#endif

	return -lk;
}

static void nl_opt_brlen_optimize(pllInstance * tree,
		partitionList * partitions, node ** branches, nlopt_algorithm algo) {
	int i;
	double minf; /* the minimum objective value, upon return */

	/* create algorithm */
	nlopt_opt opt;
	opt = nlopt_create(algo, num_branches);

	tree_data trd;
	trd.tr = tree;
	trd.pr = partitions;
	trd.branches = branches;

	/* set bounds */
	nlopt_set_lower_bounds1(opt, 0.0);
	nlopt_set_upper_bounds1(opt, 1.0);

	double * x = (double *) alloca(num_branches * sizeof(double));

	for (i = 0; i < num_branches; i++) {
		x[i] = branches[i]->z[0];
	}

	nlopt_set_min_objective(opt, myfunc, (void *) &trd);
	nlopt_set_xtol_rel(opt, epsilon);
	if (time_limit > 1e-5) {
		nlopt_set_maxtime(opt, time_limit);
	}
	//nlopt_set_maxeval(opt, 500);

	if (nlopt_optimize(opt, x, &minf) < 0) {
		printf("nlopt failed!\n");
	}
	else {
		printf("found minimum with %0.10g\n", minf);
	}

	nlopt_destroy(opt);
}

void print_model(partitionList * partitions) {
	int i;
	pInfo * part = partitions->partitionData[0];
	printf(" --- Model --- \n");
	printf("  F: ");
	for (i = 0; i < 4; i++) {
		printf("%4.2lf ", part->frequencies[i]);
	}
	printf("\n  R: ");
	for (i = 0; i < 6; i++) {
		printf("%4.2lf ", part->substRates[i]);
	}
	printf("\n  a: %4.2lf", part->alpha);
	printf("\n  l: %6.2lf\n", part->partitionLH);
	printf(" --- ----- --- \n");
}

void print_tree(pllInstance * tree, partitionList * partitions,
		FILE * foutput) {
	pllTreeToNewick(tree->tree_string, tree, partitions, tree->start->back,
	true,
	true,
	false, false, false, PLL_SUMMARIZE_LH, false,
	false);
	tree->tree_string[tree->treeStringLength - 1] = '\0';
	pllEvaluateLikelihood(tree, partitions, tree->start, true, false);
	if (foutput)
		fprintf(foutput, "%lf %s\n", tree->likelihood, tree->tree_string);
	else
		printf("%lf %s\n", tree->likelihood, tree->tree_string);
}

int main(int argc, char **argv) {
	int i, opt = 0, long_index = 0;
	int args_validate = true;

	int use_nlopt = false;
	nlopt_algorithm algo;

	char * input_file = 0;
	char * input_tree = 0;

	pllInstance * tree;
	pllQueue * parts;
	pllAlignmentData * phylip;
	partitionList * partitions;
	int number_of_taxa, seq_len, number_of_partitions;
	int partitionId = 0;

	static struct option long_options[] = {
			{ "algorithm", required_argument, 0, 'a' },
			{ "epsilon", required_argument, 0, 'e' },
			{ "help", no_argument, 0, 'h' },
			{ "input", required_argument, 0, 'i' },
			{ "limit", required_argument, 0, 'l' },
			{ "tree", required_argument, 0, 't' },
			{ 0, 0, 0, 0 } };

	while ((opt = getopt_long(argc, argv, "a:e:hi:l:t:", long_options, &long_index))
			!= -1) {
		switch (opt) {
		case 'a': {
			long algo_id;
			char * endptr;
			use_nlopt = true;
			errno = 0;
			algo_id = strtol(optarg, &endptr, 10);
			if (errno != 0) {
				printf("Algorithm %s is invalid\n", optarg);
				args_validate = false;
			}
			if (algo_id < 0 || algo_id >= NLOPT_NUM_ALGORITHMS) {
				printf("Algorithm %s is invalid\n", optarg);
				args_validate = false;

			} else {
				algo = (nlopt_algorithm) algo_id;
			}
			break;
		}
		case 'e':
			epsilon = atof(optarg);
			break;
		case 'h':
			exit_with_usage(argv[0]);
			break;
		case 'i':
			input_file = (char *) malloc(strlen(optarg) + 1);
			strcpy(input_file, optarg);
			break;
		case 'l':
			time_limit = atof(optarg);
			break;
		case 't':
			input_tree = (char *) malloc(strlen(optarg) + 1);
			strcpy(input_tree, optarg);
			break;
		}
	}

	if (!input_file) {
		printf("Input file is mandatory (-i)\n");
		args_validate = false;
	}
	if (!input_tree) {
		printf("Input tree is mandatory (-t)\n");
		args_validate = false;
	}
	if (!args_validate) {
		exit(EXIT_FAILURE);
	}

	if (use_nlopt) {
		printf("Algorithm: %s\n", nlopt_algorithm_name(algo));
	} else {
		printf("Algorithm: PLL\n");
	}

	phylip = pllParseAlignmentFile(PLL_FORMAT_PHYLIP, input_file);
	if (!phylip) {
		printf("ERROR! Cannot parse file %s\n", input_file);
		exit(EXIT_FAILURE);
	}

	char * partitionstring = (char *) malloc((size_t) 500);
	sprintf(partitionstring, "DNA, p1=1-%d", phylip->sequenceLength);
	parts = pllPartitionParseString(partitionstring);
	free(partitionstring);
	if (!parts) {
		printf("ERROR! Something failed parsing partitions for %s\n",
				input_file);
		exit(EXIT_FAILURE);
	}

	if (!pllPartitionsValidate(parts, phylip)) {
		printf("ERROR! Partitions are not valid.\n");
		exit(EXIT_FAILURE);
	}

	partitions = pllPartitionsCommit(parts, phylip);
	pllQueuePartitionsDestroy(&parts);

	number_of_taxa = phylip->sequenceCount;
	seq_len = phylip->sequenceLength;
	number_of_partitions = partitions->numberOfPartitions;

	pllAlignmentRemoveDups(phylip, partitions);

	pllInstanceAttr attr;
	attr.fastScaling = false;
	attr.randomNumberSeed = 12345;
	attr.rateHetModel = PLL_GAMMA;
	attr.saveMemory = false;
	attr.useRecom = false;
	attr.numberOfThreads = 1;

	tree = pllCreateInstance(&attr);
	tree->perGeneBranchLengths = true;

	pllNewickTree * newick = pllNewickParseFile(input_tree);
	if (!pllValidateNewick(newick)) {
		printf("ERROR Newick tree is not valid\n");
		exit(EXIT_FAILURE);
	}
	pllTreeInitTopologyNewick(tree, newick, DEFAULT_BRLEN);
	pllNewickParseDestroy(&newick);
	if (!pllLoadAlignment(tree, phylip, partitions)) {
		printf("ERROR loading alignment\n");
		exit(EXIT_FAILURE);
	}
	pllAlignmentDataDestroy(phylip);

	printf("\nSettings:\n\n");
	printf("  Input file:          %s\n", input_file);
	printf("     N: %d\n", number_of_taxa);
	printf("     L: %d\n", seq_len);
	printf("     K: %d\n", number_of_partitions);
	printf("  Input tree:          %s\n", input_tree);


//	partitions->partitionData[partitionId]->optimize
	tree->likelihoodEpsilon = 0.001;
	/** MODEL **/
	double frequencies[4] = { 0.25, 0.25, 0.25, 0.25 };
	pllInitModel(tree, partitions);

	pllSetSubstitutionRateMatrixSymmetries("0,0,0,0,0,0", partitions,
			partitionId);
	pllSetFixedBaseFrequencies(frequencies, 4, partitionId, partitions, tree);
	pllSetFixedAlpha(PLL_ALPHA_MAX, partitionId, partitions, tree);
//	pllSetFixedSubstitutionMatrix(rates, 6, partitionId, partitions, tree);
//	pllInitReversibleGTR(tree, partitions, partitionId);
	pllEvaluateLikelihood(tree, partitions, tree->start, true, false);

	num_branches = 2 * tree->mxtips - 3;
	/* extract branches */
	node ** branches = (node **) malloc(num_branches * sizeof(node *));
	branches[0] = tree->nodep[1];
	for (i = 1; i <= (2 * tree->mxtips - 2); i++) {
		tree->nodep[i]->support = 0;
		tree->nodep[i]->next->support = 0;
		tree->nodep[i]->next->next->support = 0;
	}
	node * startnode = tree->nodep[1];
	node * curnode = startnode->next->back;
	i = 0;
	while (curnode != startnode) {
		if (!curnode->back->support) {
			branches[i] = curnode;
			i++;
			curnode->back->support = 1;
		}
		curnode->support = 1;
		curnode = curnode->next->back;
	}
	assert(i == num_branches);

	/** OPTIMIZATION **/
	//print_model(partitions);
	//printf("\n");
//	pllOptimizeModelParameters(tree, partitions, 100.1);
//	pllOptimizeBranchLengths(tree, partitions, 32);
	struct timeval stop, start;
	gettimeofday(&start, NULL);
	if (use_nlopt) {
		nl_opt_brlen_optimize(tree, partitions, branches, algo);
	} else {
		pllOptimizeModelParameters(tree, partitions, 100.1);
		pllOptimizeBranchLengths(tree, partitions, 32);
	}
	gettimeofday(&stop, NULL);
	pllEvaluateLikelihood(tree, partitions, tree->start, true, false);
	__suseconds_t tt = stop.tv_sec * 1000 + stop.tv_usec / 1000
			- start.tv_sec * 1000 - start.tv_usec / 1000;
	printf("   Result: %lf (%lu ms) \n\n", tree->likelihood, tt);
	print_model(partitions);
	print_tree(tree, partitions, stdout);

	pllPartitionsDestroy(tree, &partitions);
	pllDestroyInstance(tree);

	free(branches);
	free(input_file);
	free(input_tree);

	return EXIT_SUCCESS;
}

void exit_with_usage(const char * command) {
	int i;
	printf("Usage:\n");
	printf("  %s -i input-filename -t tree-filename [-a algorithm]\n", command);
	printf("\n");
	printf("  -a, --algorithm algorithm-id       Select the optimization algorithm\n\n");
	for (i=0; i<NLOPT_NUM_ALGORITHMS; i++) {
	printf("      --algorithm %2d    %s\n", i, nlopt_algorithm_name((nlopt_algorithm) i));
	}
	printf("\n");
	printf("  -e, --epsilon xtol                 Set the relative tolerance (default: 0.001)\n");
	printf("\n");
	printf("  -h, --help                         Shows this help message\n");
	printf("\n");
	printf("  -i, --input input-filename         Set the input alignment (mandatory)\n");
	printf("\n");
	printf("  -l, --limit time-limit             Set the time limit in seconds (default: 1.0)\n");
	printf("                                     0 = unlimited\n");
	printf("\n");
	printf("  -t, --tree input-tree              Set the tree (mandatory)\n");
	printf("\n");
	printf("Examples:\n");
		printf("  (1)  %s -i input -t tree -l 5.0\n", command);
		printf("  (2)  %s -t 5 -m 3 -n myTest2\n", command);
	exit(EXIT_SUCCESS);
}
