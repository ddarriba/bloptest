/*
 * bloptest.h
 *
 *  Created on: Dec 9, 2014
 *      Author: diego
 */

#ifndef BLOPTEST_H_
#define BLOPTEST_H_

#include <pll.h>

#define true 1
#define false 0

//#define PRINT_BRANCHES 1

#define DEFAULT_BRLEN false

typedef struct {
	pllInstance * tr;
	partitionList * pr;
	node ** branches;
} tree_data;

extern unsigned num_branches;
extern double time_limit;
extern double epsilon;

void exit_with_usage(const char * command);

void print_model(partitionList * partitions);
void print_tree(pllInstance * tree, partitionList * partitions,	FILE * foutput);

#endif /* BLOPTEST_H_ */
