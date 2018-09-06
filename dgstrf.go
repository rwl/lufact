// Copyright 1988 John Gilbert and Tim Peierls
// All rights reserved.

package lufact

import (
	"fmt"
	"math"
	"strings"
)

type Desc struct {
	m, n, nnz      int
	base           int
	rowind, colptr []int
}

type GP struct {
	PivotPolicy    int // 0=none, 1=partial, 2=threshold
	PivotThreshold float64
	DropThreshold  float64
	ColFillRatio   float64
	FillRatio      float64
	ExpandRatio    float64
	ColPerm        []int
	ColPermBase    int

	UserColPerm     []int
	UserColPermBase int
}

func NewGP() *GP {
	return &GP{
		PivotPolicy:    1,
		PivotThreshold: 1,
		DropThreshold:  0,  // do not drop
		ColFillRatio:   -1, // do not limit column fill ratio
		FillRatio:      4,
		ExpandRatio:    1.2,
		ColPerm:        nil,
		ColPermBase:    0,
	}
}

type LU struct {
	luSize   int
	luNZ     []float64
	luRowInd []int
	lColPtr  []int
	uColPtr  []int

	rowPerm []int
	colPerm []int
}

func DGSTRF(gp *GP, nrow, ncol int, nzA []float64, descA *Desc) (*LU, error) {
	//var lu *LU

	var pivtRow, origRow, thisCol, othrCol int

	// Extract data from gp object.
	if gp == nil {
		return nil, fmt.Errorf("gp must not be nil")
	}

	pivotPolicy := gp.PivotPolicy
	pivotThreshold := gp.PivotThreshold
	dropThreshold := gp.DropThreshold
	colFillRatio := gp.ColFillRatio
	fillRatio := gp.FillRatio
	expandRatio := gp.ExpandRatio
	userColPerm := gp.UserColPerm
	userColPermBase := gp.UserColPermBase

	fmt.Printf("piv pol=%d piv_thr=%v drop_thr=%v col_fill_rt=%v\n", pivotPolicy, pivotThreshold, dropThreshold, colFillRatio)

	//PivotThreshold = 0.001
	//if PivotThreshold == 0.0 { PivotPolicy = 0 } // no pivoting
	//PivotPolicy = 0 // no pivoting

	// If a column permutation is specified, it must be a length ncol permutation.
	if gp.UserColPerm != nil && len(gp.UserColPerm) != ncol {
		//*info = -1
		//goto free_and_exit
		return nil, fmt.Errorf("column permutation (%v) must be a length ncol %v", len(gp.UserColPerm), ncol)
	}

	// Extract data from a's array descriptor.
	//a_m := descA.m
	nA := descA.n
	nnzA := descA.nnz
	baseA := descA.base
	colptrA := descA.colptr
	rowindA := descA.rowind

	// Convert the descriptor to 1-base if necessary.
	if baseA == 0 {
		for jcol := 0; jcol < nA+1; jcol++ {
			colptrA[jcol]++
		}
		for jcol := 0; jcol < nnzA; jcol++ {
			rowindA[jcol]++
		}
		descA.base = 1
		baseA = 1
	}

	// Allocate work arrays.
	rwork := make([]float64, nrow)
	twork := make([]float64, nrow)
	found := make([]int, nrow)
	child := make([]int, nrow)
	parent := make([]int, nrow)
	pattern := make([]int, nrow)

	// Create lu structure.
	luSize := int(float64(nnzA) * fillRatio)
	lu := &LU{
		luSize:   luSize,
		luNZ:     make([]float64, luSize),
		luRowInd: make([]int, luSize),
		uColPtr:  make([]int, ncol+1),
		lColPtr:  make([]int, ncol),
		rowPerm:  make([]int, nrow),
		colPerm:  make([]int, ncol),
	}

	// Compute max matching. We use elements of the lu structure
	// for all the temporary arrays needed.
	cmatch := make([]int, ncol)
	rmatch := make([]int, nrow)

	err := maxmatch(nrow, ncol, colptrA, rowindA, lu.lColPtr, lu.uColPtr,
		lu.rowPerm, lu.colPerm, lu.luRowInd, rmatch, cmatch)
	if err != nil {
		return nil, err
	}

	for jcol := 0; jcol < ncol; jcol++ {
		if cmatch[jcol] == 0 {
			fmt.Printf("Warning: Perfect matching not found\n")
			break
		}
	}

	//for jcol := 0; jcol < ncol; jcol++ {
	//	cmatch[jcol] = jcol + 1
	//	rmatch[jcol] = jcol + 1
	//}

	// Initialize useful values and zero out the dense vectors.
	// If we are threshold pivoting, get row counts.
	var lastlu = 0

	//lasta := colptrA[ncol] - 1
	lu.uColPtr[0] = 1

	ifill(pattern, nrow, 0)
	ifill(found, nrow, 0)
	rfill(rwork, nrow, 0)
	ifill(lu.rowPerm, nrow, 0)

	if userColPerm == nil {
		for jcol := 0; jcol < ncol; jcol++ {
			lu.colPerm[jcol] = jcol + 1
		}
	} else {
		fmt.Printf("UserColPermBase = %d\n", userColPermBase)
		for jcol := 0; jcol < ncol; jcol++ {
			lu.colPerm[jcol] = userColPerm[jcol] + (1 - userColPermBase)
		}
	}

	// Compute one column at a time.
	for jcol := 1; jcol <= ncol; jcol++ {

		// Mark pointer to new column, ensure it is large enough.
		if lastlu+nrow >= lu.luSize {
			newSize := int(float64(lu.luSize) * expandRatio)

			//fmt.Fprintf(os.Stderr, "expanding to %d nonzeros...\n", newSize)

			lu.luNZ = make([]float64, newSize) // FIXME: realloc
			lu.luRowInd = make([]int, newSize)

			lu.luSize = newSize
		}

		// Set up nonzero pattern.
		{
			jjj := lu.colPerm[jcol-1]
			for i := colptrA[jjj-1]; i < colptrA[jjj]; i++ {
				pattern[rowindA[i-1]-1] = 1
			}

			thisCol = lu.colPerm[jcol-1]
			origRow = cmatch[thisCol-1]

			pattern[origRow-1] = 2

			if lu.rowPerm[origRow-1] != 0 {
				return nil, fmt.Errorf("pivot row from max-matching already used")
			}
			// pattern[ thisCol - 1 ] = 2
		}

		// Depth-first search from each above-diagonal nonzero of column
		// jcol of A, allocating storage for column jcol of U in
		// topological order and also for the non-fill part of column
		// jcol of L.
		err := ludfs(jcol, nzA, rowindA, colptrA, &lastlu,
			lu.luRowInd, lu.lColPtr, lu.uColPtr,
			lu.rowPerm, lu.colPerm, rwork, found, parent, child)
		if err != nil {
			return nil, err
		}

		// Compute the values of column jcol of L and U in the dense
		// vector, allocating storage for fill in L as necessary.

		lucomp(jcol, &lastlu, lu.luNZ, lu.luRowInd, lu.lColPtr, lu.uColPtr,
			lu.rowPerm, lu.colPerm, rwork, found)

		//if rwork[origRow-1] == 0.0 {
		//	fmt.Printf("Warning: Matching to a zero\n")
		//
		//	for i := colptrA[jcol-1]; i < colptrA[jcol]; i++ {
		//		fmt.Printf("(%d,%v) ", rowindA[i-1], nzA[i-1])
		//		fmt.Printf(". origRow=%d\n", origRow)
		//	}
		//}

		// Copy the dense vector into the sparse data structure, find the
		// diagonal element (pivoting if specified), and divide the
		// column of L by it.
		nzCountLimit := int(colFillRatio * (float64(colptrA[thisCol] - colptrA[thisCol-1] + 1)))

		zpivot, err := lucopy(pivotPolicy, pivotThreshold, dropThreshold, nzCountLimit,
			jcol, ncol, &lastlu, lu.luNZ, lu.luRowInd, lu.lColPtr, lu.uColPtr,
			lu.rowPerm, lu.colPerm, rwork, pattern, twork)
		if err != nil {
			return nil, err
		}
		if zpivot == -1 {
			return nil, fmt.Errorf("jcol=%v", jcol)
		}

		{
			jjj := lu.colPerm[jcol-1]
			for i := colptrA[jjj-1]; i < colptrA[jjj]; i++ {
				pattern[rowindA[i-1]-1] = 0
			}

			pattern[origRow-1] = 0

			pivtRow = zpivot
			othrCol = rmatch[pivtRow-1]

			cmatch[thisCol-1] = pivtRow
			cmatch[othrCol-1] = origRow
			rmatch[origRow-1] = othrCol
			rmatch[pivtRow-1] = thisCol

			//pattern[thisCol - 1] = 0
		}

		// If there are no diagonal elements after this column, change the pivot mode.
		if jcol == nrow {
			pivotPolicy = -1
		}
	}

	// Fill in the zero entries of the permutation vector, and renumber the
	// rows so the data structure represents L and U, not PtL and PtU.
	jcol := ncol + 1
	for i := 0; i < nrow; i++ {
		if (lu.rowPerm)[i] == 0 {
			lu.rowPerm[i] = jcol
			jcol = jcol + 1
		}
	}

	for i := 0; i < lastlu; i++ {
		lu.luRowInd[i] = lu.rowPerm[lu.luRowInd[i]-1]
	}

	//fmt.Printf("rperm:\n[")
	//for i := 0; i < ncol; i++ {
	//	fmt.Printf("%d ", lu.rowPerm[i])
	//}
	//fmt.Printf("]\n")
	//
	//fmt.Printf("cperm:\n[")
	//for i := 0; i < ncol; i++ {
	//	fmt.Printf("%d ", lu.ColPerm[i])
	//}
	//fmt.Printf("]\n")

	{
		var ujj float64
		minujj := math.Inf(0)

		for jcol := 1; jcol <= ncol; jcol++ {
			ujj = math.Abs(lu.luNZ[lu.lColPtr[jcol-1]-2])
			if ujj < minujj {
				minujj = ujj
			}
		}

		//fmt.Printf(">>> last = %v, min = %v\n", ujj, minujj)
	}

	return lu, nil
}

func DGSTRS(gp *GP, trans string, n, nrhs int, lu *LU /*ia, ja int,*/, b []float64 /*, ib, jb int*/) error {
	var rwork []float64

	if gp == nil {
		return fmt.Errorf("gp must not be nil")
	}
	if nrhs != 1 {
		return fmt.Errorf("nrhs must be 1")
	}

	rwork = make([]float64, n)

	if strings.ToUpper(trans) == "N" {
		lsolve(n, lu.luNZ, lu.luRowInd, lu.lColPtr, lu.uColPtr, lu.rowPerm, lu.colPerm, b, rwork)
		usolve(n, lu.luNZ, lu.luRowInd, lu.lColPtr, lu.uColPtr, lu.rowPerm, lu.colPerm, rwork, b)
	} else if strings.ToUpper(trans) == "T" {
		utsolve(n, lu.luNZ, lu.luRowInd, lu.lColPtr, lu.uColPtr, lu.rowPerm, lu.colPerm, b, rwork)
		ltsolve(n, lu.luNZ, lu.luRowInd, lu.lColPtr, lu.uColPtr, lu.rowPerm, lu.colPerm, rwork, b)
	} else {
		return fmt.Errorf("trans %q must be N or T", trans)
	}

	return nil
}
