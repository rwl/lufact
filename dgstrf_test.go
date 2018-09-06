// Copyright 1988 John Gilbert and Tim Peierls
// All rights reserved.

package lufact

import "testing"

func TestDGSTRFS(t *testing.T) {
	var (
		rows   = 10
		cols   = 10
		arow   = []int{0, 7, 8, 1, 4, 9, 2, 9, 3, 6, 7, 8, 9, 1, 4, 5, 3, 6, 9, 0, 3, 7, 8, 0, 3, 7, 8, 1, 2, 3, 6, 9}
		acolst = []int{0, 3, 6, 8, 13, 15, 16, 19, 23, 27, 32}
		a      = []float64{2.1, 0.14, 0.09, 1.1, 0.06, 0.03, 1.7, 0.04, 1, 0.32, 0.19, 0.32, 0.44, 0.06, 1.6, 2.2, 0.32, 1.9, 0.43, 0.14, 0.19, 1.1, 0.22, 0.09, 0.32, 0.22, 2.4, 0.03, 0.04, 0.44, 0.43, 3.2}

		b = []float64{0.403, 0.28, 0.55, 1.504, 0.812, 1.32, 1.888, 1.168, 2.473, 3.695}
	)

	permC := make([]int, cols)
	for i := range permC {
		permC[i] = i
	}

	descA := &Desc{
		//Type: CSC,
		m:      cols, // transpose
		n:      rows,
		nnz:    len(a),
		base:   0,
		rowind: arow,
		colptr: acolst,
	}

	gp := NewGP()
	gp.ColPerm = permC
	gp.UserColPerm = []int{6, 5, 2, 4, 1, 9, 7, 8, 0, 3}

	lu, err := DGSTRF(gp, rows, cols, a, descA)
	if err != nil {
		t.Fatal(err)
	}

	err = DGSTRS(gp, "T", len(b), 1, lu, b)
	if err != nil {
		t.Fatal(err)
	}

	t.Logf("%v", b)
}
