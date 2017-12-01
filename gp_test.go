// Sparse Pivoting in Time Proportional to Arithmetic Operations
package lufact

import (
	"fmt"
	"math"
	"testing"
	"time"
)

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// This is a sample driver for lufact and lusolv, to illustrate how
// they are called.  It reads a matrix from a file in a format
// described below into the arrays a, arow, and acolst, then factors
// this matrix with lufact, overwriting the original values in a and
// arow.  For a description of the matrix data structure, see lufact.
//
// If the matrix is square and no error occurs in factorization, the
// vector composed of the sums of the rows of the original matrix is
// solved for this factorization using lusolv, with which solution
// the relative error of the factorization is then computed.
//
// If the problem is rectangular, the matrices PtL and U resulting
// from the factorization are written out to a file.
//
// The variables pivot and thresh may be set to do partial or
// threshold pivoting.
func TestLUFact(t *testing.T) {
	t.Skip("lufact is alpha")

	// Storage required by this routine (with overwriting):

	//a := make([]float64, lastlu)
	//dense := make([]float64, n)
	//arow := make([]int, lastlu)
	//perm := make([]int, n)
	////iw := make([]int, 3*n partial / 4*n threshold)
	//acolst := make([]int, n+1)
	//lcolst := make([]int, n)
	//ucolst := make([]int, n+1)

	var a /*, dense*/ []float64
	var arow, rperm, cperm, iw, acolst, lcolst, ucolst []int

	var maxn, maxlu, iwdim, maxnp1 int
	maxn = 12000
	maxlu = 3000000
	iwdim = 4 * maxn
	maxnp1 = maxn + 1

	a = make([]float64, maxlu)
	var thresh float64
	x := make([]float64, maxn)
	var tstart time.Time
	var tfact, tsolv, ttotl time.Duration
	var norm, tnorm /*, mflops*/ float64
	//var opcnt float64
	var pivot, nrow, ncol, lastlu int
	arow = make([]int, maxlu)
	acolst = make([]int, maxnp1)
	lcolst = make([]int, maxn)
	ucolst = make([]int, maxnp1)
	rperm = make([]int, maxn)
	iw = make([]int, iwdim)
	var lasta, nprobs /*i,*/, k, nnz, n /*, nzptr*/ int
	var nreals, nints, nrpi, n2rpi int
	var title string

	//var second float64
	//real*8    flops

	pivot = 1
	thresh = 0.1
	nprobs = 1

	// Read in the matrix.

	//     up to 60 characters           title
	//     integer nrow                  number of rows
	//     integer ncol                  number of cols
	//     (for each column:)
	//           integer nnz             number of nzs in col
	//                   integer index   row index of nonzero
	//                   real value      value of nonzero

	lasta = 1
	//read (5, 901, end=900) title
	//901       format (60a1)
	//          write (6, 902) title
	//902       format (60a1)
	//          read (5,*) nrow, ncol

	n, err := fmt.Scanf("%s\n", &title)
	if err != nil || n != 1 {
		panic(err)
	}
	n, err = fmt.Scanf("%d %d\n", &nrow, &ncol)
	if err != nil || n != 1 {
		panic(err)
	}

	if nrow > maxn || ncol > maxn {
		goto l851
	}
	for i := 1; i <= ncol; i++ {
		n, err = fmt.Scanf("%d\n", &nnz)
		if err != nil || n != 1 {
			panic(err)
		}

		//read (5,*) nnz
		acolst[i] = lasta
		if nnz+lasta-1 > maxlu {
			goto l852
		}
		for k := 1; k <= nnz; k++ {
			n, err = fmt.Scanf("%d %f\n", &arow[lasta], &a[lasta])
			if err != nil || n != 1 {
				panic(err)
			}
			//read (5,*) arow(lasta), a(lasta)
			lasta = lasta + 1
		}
	}
	acolst[ncol+1] = lasta
	lasta = lasta - 1

	//     Put the row sums of A in the vector X.

	rfill(x, nrow, 0.0)
	for k := 1; k <= ncol; k++ {
		for nzptr := acolst[k]; nzptr <= acolst[k+1]-1; nzptr++ {
			i := arow[nzptr]
			x[i] = x[i] + a[nzptr]
		}
	}

	// Perform numerical factorization (with overwriting).
	tstart = time.Now()
	err = lufact(pivot, thresh, 0, nrow, ncol, a, arow, acolst, maxlu, lastlu, a, arow, lcolst, ucolst, rperm, cperm)
	if err != nil {
		panic(err)
	}
	tfact = time.Since(tstart)

	// Set RELERR and TSOLV to appropriate values in case
	// problem is rectangular or there is a zero pivot.

	norm = 0.0
	tsolv = 0.0
	if err != nil || nrow != ncol {
		goto l500
	}

	// Find solution and compute infinity norm of X - (1,1,...,1)'.

	tstart = time.Now()
	err = lusolv(ncol, a, arow, lcolst, ucolst, rperm, cperm, x)
	if err != nil {
		panic(err)
	}
	tsolv = time.Since(tstart)

	for k := 1; k <= nrow; k++ {
		tnorm = math.Abs(x[k] - 1.0)
		if tnorm > norm {
			norm = tnorm
		}
	}

	// Print statistics.
l500: /*continue*/

	//opcnt = flops(nrow, ncol, arow, lcolst, ucolst, iw)
	//mflops = 1.e-6 * opcnt / tfact
	ttotl = tfact + tsolv
	nreals = lastlu + nrow
	nints = lastlu + 7*nrow + 2
	if pivot == 2 {
		nints = nints + nrow
	}
	nrpi = nreals + nints
	n2rpi = 2*nreals + nints
	//write(nprobs, err, nrow, ncol, lasta, lastlu, opcnt, tfact, tsolv, ttotl, mflops, norm, nreals, nints, nrpi, n2rpi)
	//l599:
	fmt.Printf("problem       %v\n", nprobs)
	fmt.Printf("err           %v\n", err)
	fmt.Printf("nrow  %v\n", nrow)
	fmt.Printf("ncol  %v\n", ncol)
	fmt.Printf("m             %v\n", lasta)
	fmt.Printf("m*            %v\n", lastlu)
	//fmt.Printf("Ops(factor)   %v\n", opcnt)
	fmt.Printf("t(factor)     %v\n", tfact)
	fmt.Printf("t(solve)      %v\n", tsolv)
	fmt.Printf("t(total)      %v\n", ttotl)
	//fmt.Printf("Mflops(factor)%v\n", )
	fmt.Printf("rel. error    %v\n", norm)
	fmt.Printf("#reals        %v\n", nreals)
	fmt.Printf("#integers     %v\n", nints)
	fmt.Printf("#r + #i       %v\n", nrpi)
	fmt.Printf("2#r + #i      %v\n", n2rpi)

	// End of loop for square problems.
	if nrow == ncol {
		goto l900
	}

	// Invert the permutation for results.

	for i := 1; i <= nrow; i++ {
		iw[rperm[i]] = i
	}

	//     Output the results (problem is rectangular).
	n = minInt(nrow, ncol)
	fmt.Println(nrow, n)
	for i := 1; i <= n-1; i++ {
		fmt.Println(ucolst[i+1] - lcolst[i] + 1)
		fmt.Println(iw[i], 1.0)
		for k = lcolst[i]; k <= ucolst[i+1]-1; k++ {
			fmt.Println(iw[arow[k]], a[k])
		}
	}
	fmt.Println(1)
	fmt.Println(iw[n], 1.0)

	fmt.Println(n, ncol)
	for i := 1; i <= ncol; i++ {
		fmt.Println(lcolst[i] - ucolst[i])
		for k := ucolst[i]; k <= lcolst[i]-1; k++ {
			fmt.Println(arow[k], a[k])
		}
	}

	goto l900

l851:
	fmt.Printf("lufact: only room for %v rows or columns", maxn)
	goto l900

l852:
	fmt.Printf("lufact: only room for %v nonzeros", maxlu)

l900: /*continue*/
}
