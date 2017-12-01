// Copyright 1988 John Gilbert and Tim Peierls
// All rights reserved.

package lufact

// lucomp computes one column of L and U in the dense vector.
//
// This routine computes column jcol of L and U (except for dividing
// through by U(jcol,jcol)) in the dense vector, which is equal to
// column jcol of A on entry.  It also allocates space in the sparse
// data structure for the fill entries in column jcol of L.
//
// Input parameters:
//   jcol    current column number.
//   rperm   row permutation P.
//           rperm(r) = s > 0 means row r of A is row s < jcol of PA.
//           rperm(r) = 0 means row r of A has not yet been used as a
//           pivot and is therefore still below the diagonal.
//   cperm   column permutation.
//
// Modified parameters:
//   lastlu  number of positions used in lurow array.
//   lu, lurow, lcolst, ucolst  nonzeros in Pt(L-I+U); see lufact for format.
//           On entry, columns 1 through jcol-1 are complete,
//           ucolst(jcol) and lcolst(jcol) point to the storage for
//           column jcol, and lurow has entries for column jcol of U and
//           the non-fill entries in column jcol of L.  The diagonal
//           element is allocated in L, not U.  On exit, storage has
//           been allocated for all of column jcol of L, and
//           ucolst(jcol+1) has been set up.  No values of column jcol
//           are in lu yet.
//   dense   column jcol as a dense vector.  On entry, column jcol of A;
//           on exit, column jcol of Pt*(U(jcol,jcol)*(L-I) + U).
//   found   found(i)=jcol if storage for position (i,jcol) has been
//           allocated in the sparse data structure.
//   flops   flop count
//
//           Both dense and found are indexed according to the row
//           numbering of A, not PA.
func lucomp(jcol, lastlu int, lu []float64, lurow, lcolst, ucolst, rperm, cperm []int, dense []float64, found /*, pattern*/ []int /*, flops float64*/) {
	// Local variables:
	//   nzuptr                pointer to current nonzero PtU(krow,jcol).
	//   nzuend, nnzu, nzuind  used to compute nzuptr.
	//   krow                  row index of current nonzero (according to A,
	//                         not PA), which is PtU(krow,jcol) = U(kcol,jcol).
	//   kcol                  rperm(krow), that is, row index of current
	//                         nonzero according to PA.
	//   ukj                   value of PtU(krow,jcol).
	//   nzlptr                pointer to PtL(irow,kcol), being used for update.
	//   nzlst, nzlend         used to compute nzlptr.
	//   irow                  row index in which update is taking place,
	//                         according to A (not PA).

	//    For each krow with PtU(krow,jcol) != 0, in reverse postorder, use
	//    column kcol = rperm(krow) of L to update the current column.
	nzuend := lcolst[jcol]
	nnzu := nzuend - ucolst[jcol]
	if nnzu != 0 {
		for nzuind := 1; nzuind <= nnzu; nzuind++ {
			nzuptr := nzuend - nzuind
			krow := lurow[nzuptr]
			kcol := rperm[krow]
			ukj := dense[krow]
			//if pattern[rperm[krow]] == 0 {
			//	ukj = 0
			//}

			// For each irow with PtL(irow,kcol) != 0, update PtL(irow,jcol) or PtU(irow,jcol)

			nzlst := lcolst[kcol]
			nzlend := ucolst[kcol+1] - 1
			if nzlend < nzlst {
				ucolst[jcol+1] = lastlu + 1
				return
			}
			for nzlptr := nzlst; nzlptr <= nzlend; nzlptr++ {
				irow := lurow[nzlptr]
				dense[irow] = dense[irow] - ukj*lu[nzlptr]

				// If this is a new nonzero in L, allocate storage for it.
				if found[irow] != jcol {
					found[irow] = jcol
					lastlu = lastlu + 1
					lurow[lastlu] = irow
				}
			}
		}
	}
	ucolst[jcol+1] = lastlu + 1
	return
}
