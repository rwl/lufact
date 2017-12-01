// Copyright 1988 John Gilbert and Tim Peierls
// All rights reserved.

package lufact

import "fmt"

// maxmatch does maximum matching
//
// maxmatch uses depth-first search to find an augmenting path from
// each column node to get the maximum matching.
//
// Alex Pothen and Chin-Ju Fan, Penn State University, 1988
// last modifed: Alex Pothen July 1990
// last bcs modifications:  John Lewis, Sept. 1990
//
// input variables :
//
//    nrows -- number of row nodes in the graph.
//    ncols -- number of column nodes in the graph.
//    colstr, rowind -- adjacency structure of graph, stored by
//                      columns
//
// output variables :
//
//    rowset -- describe the matching.
//              rowset (row) = col > 0 means column "col" is matched
//                                     to row "row"
//                           = 0       means "row" is an unmatched
//                                     node.
//    colset -- describe the matching.
//              colset (col) = row > 0 means row "row" is matched to
//                             column "col"
//                           = 0       means "col" is an unmatched
//                                     node.
func maxmatch(nrows, ncols int, colstr, rowind, prevcl, prevrw, marker, tryrow, nxtchp, rowset, colset []int) error {

	// Working variables :
	//
	//     prevrw (ncols) -- pointer toward the root of the depth-first
	//                       search from a column to a row.
	//     prevcl (ncols) -- pointer toward the root of the depth-first
	//                       search from a column to a column.
	//                       the pair (prevrw,prevcl) represent a
	//                       matched pair.
	//     marker (nrows) -- marker (row) <= the index of the root of the
	//                       current depth-first search.  row has been
	//                       visited in current pass when equality holds.
	//     tryrow (ncols) -- tryrow (col) is a pointer into rowind to
	//                       the next row to be explored from column col
	//                       in the depth-first search.
	//     nxtchp (ncols) -- nxtchp (col) is a pointer into rowind to the
	//                       next row to be explored from column col for
	//                       the cheap assignment.  set to -1 when
	//                       all rows have been considered for
	//                       cheap assignment

	var row int
	for nodec := 1; nodec <= ncols; nodec++ {
		// Initialize node 'col' as the root of the path.
		col := nodec
		prevrw[col] = 0
		prevcl[col] = 0
		nxtchp[col] = colstr[col]
		// Main loop begins here. Each time through, try to find a
		// cheap assignment from node col.
	l100:
		nextrw := nxtchp[col]
		lastrw := colstr[col+1] - 1

		if nextrw > 0 {
			for xrow := nextrw; xrow <= lastrw; xrow++ {
				row = rowind[xrow]
				if rowset[row] == 0 {
					goto l400
				}
			}

			// Mark column when all adjacent rows have been
			// considered for cheap assignment.
			nxtchp[col] = -1

		}

		// Each time through, take a step forward if possible, or
		// backtrack if not .  Quit when backtracking takes us back
		// to the beginning of the search.

		tryrow[col] = colstr[col]
		nextrw = tryrow[col]
		//lastrw = colstr [col+1] - 1

		if lastrw >= nextrw {
			for xrow := nextrw; xrow <= lastrw; xrow++ {
				// next line inserted by Alex Pothen, July 1990
				// ii  = xrow
				row = rowind[xrow]
				if marker[row] < nodec {

					// Row is unvisited yet for this pass.
					// Take a forward step.

					tryrow[col] = xrow + 1
					marker[row] = nodec
					nxtcol := rowset[row]

					if nxtcol < 0 {
						return fmt.Errorf("maxmatch : search reached a forbidden column")
					} else if nxtcol == col {
						return fmt.Errorf("maxmatch : search followed a matching edge")
					} else if nxtcol > 0 {

						// The forward step led to a matched row
						// try to extend augmenting path from
						// the column matched by this row.

						prevcl[nxtcol] = col
						prevrw[nxtcol] = row
						tryrow[nxtcol] = colstr[nxtcol]
						col = nxtcol
						goto l100

					} else {
						// Unmatched row
						goto l400
					}

				}
				//l300:
			}
		}

		// No forward step -- backtrack.
		// If we backtrack all the way, the search is done

		col = prevcl[col]
		if col > 0 {
			goto l100
		} else {
			goto l600
		}

		// Update the matching by alternating the matching
		// edge backward toward the root.
	l400:
		rowset[row] = col
		prow := prevrw[col]
		pcol := prevcl[col]

	l500:
		if pcol > 0 {
			if rowset[prow] != col {
				return fmt.Errorf("maxmatch : pointer toward root disagrees with matching. prevcl[%v]=%v but colset[%v]=%v", col, row, row, rowset[row])
			}
			rowset[prow] = pcol
			col = pcol
			prow = prevrw[pcol]
			pcol = prevcl[pcol]
			goto l500
		}
	}
l600:
	// Compute the matching from the view of column nodes.
	for row := 1; row <= nrows; row++ {
		col := rowset[row]
		if col > 0 {
			colset[col] = row
		}
	}
	return nil
}
