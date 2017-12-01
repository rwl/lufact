// Copyright 1988 John Gilbert and Tim Peierls
// All rights reserved.

package lufact

import "fmt"

// lsolve solves lower triangular system.
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
func lsolve(n int, lu []float64, lurow, lcolst, ucolst, rperm, cperm []int, b, x []float64) error {
	if n <= 0 {
		return fmt.Errorf("lsolve called with nonpositive n = %v", n)
	}

	// Check that rperm is really a permutation.
	for i := 1; i <= n; i++ {
		x[i] = 0
	}
	for i := 1; i <= n; i++ {
		if rperm[i] < 1 || rperm[i] > n {
			return fmt.Errorf("lsolve, rpermutation is illegal in position i = %v", rperm[i])
		}
		if x[rperm[i]] != 0.0 {
			return fmt.Errorf("lsolve, rpermutation is illegal in position i = %v", rperm[i])
		}
		x[rperm[i]] = 1.0
	}

	// Check that cperm is really a permutation.
	/*for i := 1; n; i++ {
		x[i] = 0
	}
	for i := 1; n; i++ {
		if cperm[i] < 1 || cperm[i] > n {
			return fmt.Errorf("lsolve, cpermutation is illegal in position i = %v", cperm[i])
		}
		if x[cperm[i]] != 0.0 {
			return fmt.Errorf("lsolve, cpermutation is illegal in position i = %v", cperm[i])
		}
		x[cperm[i]] = 1.0
	}*/

	// Solve the system.
	for i := 1; i <= n; i++ {
		x[rperm[i]] = b[i]
	}

	for j := 1; j <= n; j++ {
		nzst := lcolst[j]
		nzend := ucolst[j+1] - 1
		if nzst < 1 || nzst > nzend+1 {
			return fmt.Errorf("lsolve, inconsistent column of L: j=%v nzst=%v, nzend=%v", j, nzst, nzend)
		}
		if nzst <= nzend {
			for nzptr := nzst; nzptr <= nzend; nzptr++ {
				i := lurow[nzptr]
				if i <= j || i > n {
					return fmt.Errorf("lsolve, illegal row i in column j of L: i=%v, j=%v, nzptr=%v", i, j, nzptr)
				}
				x[i] = x[i] - lu[nzptr]*x[j]
			}
		}
	}

	return nil
}
