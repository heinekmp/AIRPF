/*
 * typedefs.c
 *
 *  Created on: 5.5.2017
 *      Author: heine
 */

#include <stdlib.h>
#include "typedefs.h"

void delete_process( struct Process p ) {
  free(p.y);
}
