package lufact

import (
	"fmt"
	"math"
	"strings"
)

type desc struct {
	m, n, nnz      int
	base           int
	rowind, colptr []int
}

type gp_t struct {
	pivot_policy    int // 0=none, 1=partial, 2=threshold
	pivot_threshold float64
	drop_threshold  float64
	col_fill_ratio  float64
	fill_ratio      float64
	expand_ratio    float64
	col_perm        []int
	col_perm_base   int

	user_col_perm      []int
	user_col_perm_base int
}

func NewGP() *gp_t {
	return &gp_t{
		pivot_policy:    1,
		pivot_threshold: 1,
		drop_threshold:  0,
		col_fill_ratio:  -1,
		fill_ratio:      4,
		expand_ratio:    1.2,
		col_perm:        nil,
		col_perm_base:   0,
	}
}

type lu_t struct {
	lu_size   int
	lu_nz     []float64
	lu_rowind []int
	l_colptr  []int
	u_colptr  []int

	row_perm []int
	col_perm []int
}

func DGSTRF(gp *gp_t, nrow, ncol int, a_nz []float64, desc_a *desc) (*lu_t, error) {
	//var lu *lu_t

	var pivt_row, orig_row, this_col, othr_col int

	// Extract data from gp object.
	if gp == nil {
		return nil, fmt.Errorf("gp must not be nil")
	}

	pivot_policy := gp.pivot_policy
	pivot_threshold := gp.pivot_threshold
	drop_threshold := gp.drop_threshold
	col_fill_ratio := gp.col_fill_ratio
	fill_ratio := gp.fill_ratio
	expand_ratio := gp.expand_ratio
	user_col_perm := gp.user_col_perm
	user_col_perm_base := gp.user_col_perm_base

	fmt.Printf("piv pol=%d piv_thr=%v drop_thr=%v col_fill_rt=%v\n", pivot_policy, pivot_threshold, drop_threshold, col_fill_ratio)

	//pivot_threshold = 0.001
	//if pivot_threshold == 0.0 { pivot_policy = 0 } // no pivoting
	//pivot_policy = 0 // no pivoting

	// If a column permutation is specified, it must be a length ncol permutation.
	if gp.user_col_perm != nil && len(gp.user_col_perm) != ncol {
		//*info = -1
		//goto free_and_exit
		return nil, fmt.Errorf("column permutation (%v) must be a length ncol %v", len(gp.user_col_perm), ncol)
	}

	// Extract data from a's array descriptor.
	//a_m := desc_a.m
	a_n := desc_a.n
	a_nnz := desc_a.nnz
	a_base := desc_a.base
	a_colptr := desc_a.colptr
	a_rowind := desc_a.rowind

	// Convert the descriptor to 1-base if necessary.
	if a_base == 0 {
		for jcol := 0; jcol < a_n+1; jcol++ {
			a_colptr[jcol]++
		}
		for jcol := 0; jcol < a_nnz; jcol++ {
			a_rowind[jcol]++
		}
		desc_a.base = 1
		a_base = 1
	}

	// Allocate work arrays.
	rwork := make([]float64, nrow)
	twork := make([]float64, nrow)
	found := make([]int, nrow)
	child := make([]int, nrow)
	parent := make([]int, nrow)
	pattern := make([]int, nrow)

	// Create lu structure.
	lu_size := int(float64(a_nnz) * fill_ratio)
	lu := &lu_t{
		lu_size:   lu_size,
		lu_nz:     make([]float64, lu_size),
		lu_rowind: make([]int, lu_size),
		u_colptr:  make([]int, ncol+1),
		l_colptr:  make([]int, ncol),
		row_perm:  make([]int, nrow),
		col_perm:  make([]int, ncol),
	}

	// Compute max matching. We use elements of the lu structure
	// for all the temporary arrays needed.
	cmatch := make([]int, ncol)
	rmatch := make([]int, nrow)

	err := maxmatch(nrow, ncol, a_colptr, a_rowind, lu.l_colptr, lu.u_colptr,
		lu.row_perm, lu.col_perm, lu.lu_rowind, rmatch, cmatch)
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
	lastlu := 0

	//lasta := a_colptr[ncol] - 1
	lu.u_colptr[0] = 1

	ifill(pattern, nrow, 0)
	ifill(found, nrow, 0)
	rfill(rwork, nrow, 0)
	ifill(lu.row_perm, nrow, 0)

	if user_col_perm == nil {
		for jcol := 0; jcol < ncol; jcol++ {
			lu.col_perm[jcol] = jcol + 1
		}
	} else {
		fmt.Printf("user_col_perm_base = %d\n", user_col_perm_base)
		for jcol := 0; jcol < ncol; jcol++ {
			lu.col_perm[jcol] = user_col_perm[jcol] + (1 - user_col_perm_base)
		}
	}

	// Compute one column at a time.
	for jcol := 1; jcol <= ncol; jcol++ {

		// Mark pointer to new column, ensure it is large enough.
		if lastlu+nrow >= lu.lu_size {
			var new_size int = int(float64(lu.lu_size) * expand_ratio)

			//fmt.Fprintf(os.Stderr, "expanding to %d nonzeros...\n", new_size)

			lu.lu_nz = make([]float64, new_size)
			lu.lu_rowind = make([]int, new_size)

			lu.lu_size = new_size
		}

		// Set up nonzero pattern.
		{
			jjj := lu.col_perm[jcol-1]
			for i := a_colptr[jjj-1]; i < a_colptr[jjj]; i++ {
				pattern[a_rowind[i-1]-1] = 1
			}

			this_col = lu.col_perm[jcol-1]
			orig_row = cmatch[this_col-1]

			pattern[orig_row-1] = 2

			if lu.row_perm[orig_row-1] != 0 {
				return nil, fmt.Errorf("pivot row from max-matching already used")
			}
			// pattern[ this_col - 1 ] = 2
		}

		// Depth-first search from each above-diagonal nonzero of column
		// jcol of A, allocating storage for column jcol of U in
		// topological order and also for the non-fill part of column
		// jcol of L.
		err := ludfs(jcol, a_nz, a_rowind, a_colptr, lastlu,
			lu.lu_rowind, lu.l_colptr, lu.u_colptr,
			lu.row_perm, lu.col_perm, rwork, found, parent, child)
		if err != nil {
			return nil, err
		}

		// Compute the values of column jcol of L and U in the dense
		// vector, allocating storage for fill in L as necessary.

		lucomp(jcol, lastlu, lu.lu_nz, lu.lu_rowind, lu.l_colptr, lu.u_colptr,
			lu.row_perm, lu.col_perm, rwork, found)

		//if rwork[orig_row-1] == 0.0 {
		//	fmt.Printf("Warning: Matching to a zero\n")
		//
		//	for i := a_colptr[jcol-1]; i < a_colptr[jcol]; i++ {
		//		fmt.Printf("(%d,%v) ", a_rowind[i-1], a_nz[i-1])
		//		fmt.Printf(". orig_row=%d\n", orig_row)
		//	}
		//}

		// Copy the dense vector into the sparse data structure, find the
		// diagonal element (pivoting if specified), and divide the
		// column of L by it.
		nz_count_limit := int(col_fill_ratio * (float64(a_colptr[this_col] - a_colptr[this_col-1] + 1)))

		zpivot, err := lucopy(pivot_policy, pivot_threshold, drop_threshold, nz_count_limit,
			jcol, ncol, &lastlu, lu.lu_nz, lu.lu_rowind, lu.l_colptr, lu.u_colptr,
			lu.row_perm, lu.col_perm, rwork, pattern, twork)
		if err != nil {
			return nil, err
		}
		if zpivot == -1 {
			return nil, fmt.Errorf("jcol=%v", jcol)
		}

		{
			jjj := lu.col_perm[jcol-1]
			for i := a_colptr[jjj-1]; i < a_colptr[jjj]; i++ {
				pattern[a_rowind[i-1]-1] = 0
			}

			pattern[orig_row-1] = 0

			pivt_row = zpivot
			othr_col = rmatch[pivt_row-1]

			cmatch[this_col-1] = pivt_row
			cmatch[othr_col-1] = orig_row
			rmatch[orig_row-1] = othr_col
			rmatch[pivt_row-1] = this_col

			//pattern[this_col - 1] = 0
		}

		// If there are no diagonal elements after this column, change the pivot mode.
		if jcol == nrow {
			pivot_policy = -1
		}
	}

	// Fill in the zero entries of the permutation vector, and renumber the
	// rows so the data structure represents L and U, not PtL and PtU.
	jcol := ncol + 1
	for i := 0; i < nrow; i++ {
		if (lu.row_perm)[i] == 0 {
			lu.row_perm[i] = jcol
			jcol = jcol + 1
		}
	}

	for i := 0; i < lastlu; i++ {
		lu.lu_rowind[i] = lu.row_perm[lu.lu_rowind[i]-1]
	}

	//fmt.Printf("rperm:\n[")
	//for i := 0; i < ncol; i++ {
	//	fmt.Printf("%d ", lu.row_perm[i])
	//}
	//fmt.Printf("]\n")
	//
	//fmt.Printf("cperm:\n[")
	//for i := 0; i < ncol; i++ {
	//	fmt.Printf("%d ", lu.col_perm[i])
	//}
	//fmt.Printf("]\n")

	{
		var ujj float64
		minujj := math.Inf(0)

		for jcol := 1; jcol <= ncol; jcol++ {
			ujj = math.Abs(lu.lu_nz[lu.l_colptr[jcol-1]-2])
			if ujj < minujj {
				minujj = ujj
			}
		}

		//fmt.Printf(">>> last = %v, min = %v\n", ujj, minujj)
	}

	return lu, nil
}

func DGSTRS(gp *gp_t, trans string, n, nrhs int, lu *lu_t /*ia, ja int,*/, b []float64 /*, ib, jb int*/) error {
	var rwork []float64

	if gp == nil {
		return fmt.Errorf("gp must not be nil")
	}
	if nrhs != 1 {
		return fmt.Errorf("nrhs must be 1")
	}

	rwork = make([]float64, n)

	if strings.ToUpper(trans) == "N" {
		lsolve(n, lu.lu_nz, lu.lu_rowind, lu.l_colptr, lu.u_colptr, lu.row_perm, lu.col_perm, b, rwork)
		usolve(n, lu.lu_nz, lu.lu_rowind, lu.l_colptr, lu.u_colptr, lu.row_perm, lu.col_perm, rwork, b)
	} else if strings.ToUpper(trans) == "T" {
		utsolve(n, lu.lu_nz, lu.lu_rowind, lu.l_colptr, lu.u_colptr, lu.row_perm, lu.col_perm, b, rwork)
		ltsolve(n, lu.lu_nz, lu.lu_rowind, lu.l_colptr, lu.u_colptr, lu.row_perm, lu.col_perm, rwork, b)
	} else {
		return fmt.Errorf("trans %q must be N or T", trans)
	}

	return nil
}
