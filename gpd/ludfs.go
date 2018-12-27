// Code generated with gpgen. DO NOT EDIT.

// Copyright 1988 John Gilbert and Tim Peierls
// All rights reserved.

package gpd

import "fmt"

// ludfs :  Depth-first search to allocate storage for U
//
// Input parameters:
//   jcol             current column number.
//   a, arow, acolst  the matrix A; see lufact for format.
//   rperm            row permutation P.
//                    perm(r) = s > 0 means row r of A is row s < jcol of PA.
//                    perm(r) = 0 means row r of A has not yet been used as a
//                    pivot and is therefore still below the diagonal.
//   cperm            column permutation.
//
// Modified parameters (see below for exit values):
//   lastlu           last used position in lurow array.
//   lurow, lcolst, ucolst  nonzero structure of Pt(L-I+U);
//                          see lufact for format.
//   dense            current column as a dense vector.
//   found            integer array for marking nonzeros in this column of
//                    Pt(L-I+U) that have been allocated space in lurow.
//                    Also, marks reached columns in depth-first search.
//                    found(i)=jcol if i was found in this column.
//   parent           parent(i) is the parent of vertex i in the dfs,
//                    or 0 if i is a root of the search.
//   child            child(i) is the index in lurow of the next unexplored
//                    child of vertex i.
//                    Note that parent and child are also indexed according to
//                    the vertex numbering of A, not PA; thus child(i) is
//                    the position of a nonzero in column rperm(i),
//                    not column i.
//
// Output parameters:
//   error            0 if successful, 1 otherwise
//
// On entry:
//   found(*)<jcol
//   dense(*)=0.0
//   ucolst(jcol)=lastlu+1 is the first free index in lurow.
//
// On exit:
//   found(i)=jcol iff i is a nonzero of column jcol of PtU or
//     a non-fill nonzero of column jcol of PtL.
//   dense(*)=column jcol of A.
//     Note that found and dense are kept according to the row
//     numbering of A, not PA.
//   lurow has the rows of the above-diagonal nonzeros of col jcol of U in
//     reverse topological order, followed by the non-fill nonzeros of col
//     jcol of L and the diagonal elt of U, in no particular order.
//     These rows also are numbered according to A, not PA.
//   lcolst(jcol) is the index of the first nonzero in col j of L.
//   lastlu is the index of the last non-fill nonzero in col j of L.
func ludfs(jcol int, a []float64, arow, acolst []int, lastlu *int, lurow, lcolst, ucolst, rperm, cperm []int, dense []float64, found, parent, child []int) error {
	// Depth-first search through columns of L from each nonzero of
	// column jcol of A that is above the diagonal in PA.

	// For each krow such that A(krow,jcol) is nonzero do...

	// Range of indices in arow for column jcol of A.
	nzast := acolst[cperm[jcol-off]-off]
	nzaend := acolst[cperm[jcol-off]] //+1-off]

	if nzaend < nzast {
		return fmt.Errorf("ludfs, negative length for column %v of A. nzast=%v nzend=%v", jcol, nzast, nzaend)
	}
	nzaend = nzaend - 1
	for nzaptr := nzast - 1; nzaptr < nzaend; nzaptr++ { // pointer to current position in arow (zero based)
		// Current vertex in depth-first search (numbered according to A, not PA) (zero based).
		krow := arow[nzaptr] - 1

		// Copy A(krow,jcol) into the dense vector. If above diagonal in
		// PA, start a depth-first search in column rperm(krow) of L.

		dense[krow] = a[nzaptr]
		if rperm[krow] == 0 || found[krow] == jcol || dense[krow] == 0 {
			continue
		}
		parent[krow] = 0
		found[krow] = jcol
		chdptr := lcolst[rperm[krow]-off] // Index of current child of current vertex.

		// The main depth-first search loop starts here.
		// repeat
		//   if krow has a child that is not yet found
		//   then step forward
		//   else step back
		// until a step back leads to 0
	l100:
		// Look for an unfound child of krow.
		chdend := ucolst[rperm[krow]] // Next index after last child of current vertex.

	l200:
		if chdptr < chdend {
			// Possible next vertex in depth-first search (zero based).
			nextk := lurow[chdptr-off] - 1
			chdptr = chdptr + 1
			if rperm[nextk] == 0 {
				goto l200
			}
			if found[nextk] == jcol {
				goto l200
			}
			// Take a step forward.

			//l300:
			child[krow] = chdptr
			parent[nextk] = krow + 1
			krow = nextk
			found[krow] = jcol
			chdptr = lcolst[rperm[krow]-off]
			goto l100
		}
		// Take a step back.

		// Allocate space for U(rperm(k),jcol) = PtU(krow,jcol) in the sparse data structure.
		*lastlu = *lastlu + 1
		lurow[*lastlu-off] = krow + 1
		krow = parent[krow] - 1
		if krow >= 0 {
			chdptr = child[krow]
			goto l100
		}
		// The main depth-first search loop ends here.
	}
	// Close off column jcol of U and allocate space for the non-fill
	// entries of column jcol of L.
	// The diagonal element goes in L, not U, until we do the column
	// division at the end of the major step.

	lcolst[jcol-off] = *lastlu + 1
	for nzaptr := nzast - 1; nzaptr < nzaend; nzaptr++ {
		krow := arow[nzaptr]
		if rperm[krow-off] == 0 {
			found[krow-off] = jcol
			*lastlu += 1
			lurow[*lastlu-off] = krow
		}
	}

	return nil
}
