/*
 * typedefs.h
 *
 *  Created on: 5.5.2017
 *      Author: heine
 */

#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_

struct Process {
  double *y;
  long n;
  int dim;
};

struct Estimates {
  double nzer;
  double *mu;
  double w;
  //double minx;
  //double maxx;
};

void delete_process(struct Process p);

#endif /* TYPEDEFS_H_ */
