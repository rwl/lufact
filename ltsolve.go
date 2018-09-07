// Copyright 1988 John Gilbert and Tim Peierls
// All rights reserved.

package gp

import "fmt"

// ltsolve: Modified from lsolve to solve with L transpose.
// Sivan: removed error checking marked by cs comments.

// ltsolve solves lower triangular systems.
//
// This routine takes an LU factorization from lufact (i.e. P, L, U with
// PA = LU) and solves Lx = Pb for x.  There is nothing clever at all
// about sparse right-hand sides here; we always look at every nonzero
// of L.  We do make some checks for consistency of the LU data
// structure.
//
// Input parameters:
//   n    Dimension of the system.
//   lu, lurow, lcolst, ucolst, rperm, cperm  LU factorization
//   b    Right-hand side, as a dense n-vector.
//
// Output parameter:
//   x    Solution, as a dense n-vector.
//   error 0 if successful, 1 otherwise
func ltsolve(n int, lu []float64, lurow, lcolst, ucolst, rperm, cperm []int, b, x []float64) error {
	if n <= 0 {
		return fmt.Errorf("ltsolve called with nonpositive n=%v", n)
	}

	// Check that rperm is really a permutation.
	//
	//      do 10 i = 1, n
	//          x(i) = 0.0
	//10        continue
	//      do 20 i = 1, n
	//          if (rperm(i) .lt. 1  .or.  rperm(i) .gt. n) { "ltsolve, rpermutation is illegal in position i =" }
	//          if (x(rperm(i)) .ne. 0.0) goto 803
	//          x(rperm(i)) = 1.0
	//20        continue
	//
	// Check that cperm is really a permutation.
	//
	//      do 110 i = 1, n
	//          x(i) = 0.0
	//110       continue
	//      do 120 i = 1, n
	//          if (cperm(i) .lt. 1  .or.  cperm(i) .gt. n) { "lsolve, cpermutation is illegal in position i =" }
	//          if (x(cperm(i)) .ne. 0.0) goto 804
	//          x(cperm(i)) = 1.0
	//120        continue

	// Solve the system.
	for i := 1; i <= n; i++ {
		x[i-off] = b[i-off]
	}

	for j := n; j >= 1; j-- {
		nzst := lcolst[j-off]
		nzend := ucolst[j+1-off] - 1
		if nzst < 1 || nzst > nzend+1 {
			return fmt.Errorf("ltsolve, inconsistent column of L: j=%v, nzst=%v, nzend=%v", j, nzst, nzend)
		}
		if nzst <= nzend {
			for nzptr := nzst; nzptr <= nzend; nzptr++ {
				i := lurow[nzptr-off]
				if i <= j || i > n {
					return fmt.Errorf("ltsolve, illegal row i in column j of L: i=%v, j=%v, nzptr=%v", i, j, nzptr)
				}
				x[j-off] = x[j-off] - lu[nzptr-off]*x[i-off]
			}
		}
	}

	for i := 1; i <= n; i++ {
		b[i-off] = x[i-off]
	}

	for i := 1; i <= n; i++ {
		//x[rperm[i-off]-off] = b[i-off]
		x[i-off] = b[rperm[i-off]-off]
	}
	return nil
}
