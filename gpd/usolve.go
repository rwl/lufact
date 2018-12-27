// Code generated with gpgen. DO NOT EDIT.

// Copyright 1988 John Gilbert and Tim Peierls
// All rights reserved.

package gpd

import "fmt"

// usolve solves the upper triangular system.
//
// This routine takes an LU factorization from lufact (i.e. L, U
// with PA = LU) and solves Ux = b for x.  Note that P is not used
// and is not a parameter.  There is nothing clever at all about
// sparse right-hand sides here; we always look at every nonzero of U.
// We do make some checks for consistency of the LU data structure.
//
// Input parameters:
//   n    Dimension of the system.
//   lu, lurow, lcolst, ucolst  LU factorization; see lufact for format.
//   b    Right-hand side, as a dense n-vector.
//
// Output parameter:
//   x    Solution, as a dense n-vector.
//   error 0 if successful, 1 otherwise
func usolve(n int, lu []float64, lurow, lcolst, ucolst, rperm, cperm []int, b, x []float64) error {
	if n <= 0 {
		return fmt.Errorf("usolve called with nonpositive n=%v", n)
	}
	for i := 1; i <= n; i++ {
		x[i-off] = b[i-off]
	}

	for jj := 1; jj <= n; jj++ {
		j := n + 1 - jj
		nzst := ucolst[j-off]
		nzend := lcolst[j-off] - 1
		if nzst < 1 || nzst > nzend {
			return fmt.Errorf("usolve, inconsistent column of U: j=%v, nzst=%v, nzend=%v", j, nzst, nzend)
		}
		if lurow[nzend-off] != j {
			return fmt.Errorf("usolve, diagonal elt of col j is not in last place: j=%v, nzend=%v, lurow[nzend]=%v", j, nzend, lurow[nzend-off])
		}
		if lu[nzend-off] == 0 {
			return fmt.Errorf("usolve, zero diagonal element in column j=%v", j)
		}
		x[j-off] = x[j-off] / lu[nzend-off]
		nzend = nzend - 1
		if nzst > nzend {
			goto l150
		}
		for nzptr := nzst; nzptr <= nzend; nzptr++ {
			i := lurow[nzptr-off]
			if i <= 0 || i >= j {
				return fmt.Errorf("usolve, illegal row i in column j of U: i=%v, j=%v, nzptr=%v", i, j, nzptr)
			}
			x[i-off] -= lu[nzptr-off] * x[j-off]
		}
	l150:
	}

	for i := 1; i <= n; i++ {
		b[i-off] = x[i-off]
	}
	for i := 1; i <= n; i++ {
		x[cperm[i-off]-off] = b[i-off]
	}

	return nil
}

// utsolve solves the upper triangular system.
//
// This routine takes an LU factorization from lufact (i.e. L, U
// with PA = LU) and solves Ux = b for x.  Note that P is not used
// and is not a parameter.  There is nothing clever at all about
// sparse right-hand sides here; we always look at every nonzero of U.
// We do make some checks for consistency of the LU data structure.
//
// Input parameters:
//   n    Dimension of the system.
//   lu, lurow, lcolst, ucolst  LU factorization; see lufact for format.
//   b    Right-hand side, as a dense n-vector.
//
// Output parameter:
//   x    Solution, as a dense n-vector.
//   error 0 if successful, 1 otherwise
func utsolve(n int, lu []float64, lurow, lcolst, ucolst, rperm, cperm []int, b, x []float64) error {
	if n <= 0 {
		return fmt.Errorf("utsolve called with nonpositive n=%v", n)
	}

	//     do 60 i = 1, n
	//         x(rperm(i)) = b(i)
	//60        continue
	//
	//     do 50 i = 1, n
	//         x(i) = b(i)
	//50        continue

	for i := 1; i <= n; i++ {
		x[i-off] = b[cperm[i-off]-off]
	}

	for j := 1; j <= n; j++ {
		nzst := ucolst[j-off]
		nzend := lcolst[j-off] - 1
		if nzst < 1 || nzst > nzend {
			return fmt.Errorf("utsolve, inconsistent column of U: j=%v, nzst=%v, nzend=%v", j, nzst, nzend)
		}
		if lurow[nzend-off] != j {
			return fmt.Errorf("utsolve, diagonal elt of col j is not in last place: j=%v, nzend=%v, lurow[nzend]=%v", j, nzend, lurow[nzend-off])
		}
		if lu[nzend-off] == 0 {
			return fmt.Errorf("utsolve, zero diagonal element in column j=%v", j)
		}
		nzend = nzend - 1
		if nzst > nzend {
			goto l150
		}
		for nzptr := nzst; nzptr <= nzend; nzptr++ {
			i := lurow[nzptr-off]
			if i <= 0 || i >= j {
				return fmt.Errorf("utsolve, illegal row i in column j of U: i=%v, j=%v, nzptr=%v", i, j, nzptr)
			}
			x[j-off] -= lu[nzptr-off] * x[i-off]
		}
	l150:
		x[j-off] = x[j-off] / lu[nzend+1-off]
	}
	//l200:

	return nil
}
