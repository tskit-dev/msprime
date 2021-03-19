/*
** Copyright (C) 2019 University of Oxford
**
** This file is part of msprime.
**
** msprime is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** msprime is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with msprime.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __LIKELIHOOD_H__
#define __LIKELIHOOD_H__

#include <tskit.h>

int msp_unnormalised_log_likelihood_mut(tsk_treeseq_t *ts, double mu, double *lik);
int msp_log_likelihood_arg(tsk_treeseq_t *ts, double r, double Ne, double *lik);

#endif /*__LIKELIHOOD_H__*/
