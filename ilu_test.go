// Copyright 1988 John Gilbert and Tim Peierls
// All rights reserved.

package gp_test

import (
	"os"
	"testing"

	"math"

	"github.com/rwlincoln/lufact"
)

func TestMain(m *testing.M) {
	gp.Logger = os.Stdout
	m.Run()
}

func TestFactor(t *testing.T) {
	// Compute the matrix and right-hand-side vector that define
	// the linear system, Ax = b.
	n, rowind, colst, nzA := lhr01()

	// Set exact solution; then compute right-hand-side vector.
	x0 := make([]float64, n)
	for i := range x0 {
		x0[i] = 1
	}

	for i, opts := range [][]gp.OptFunc{
		{
			//gp.DropThreshold(0),
			//gp.PartialPivoting(1),
			//gp.ColFillRatio(-1),
			//gp.ColPerm(nil), // natural
			gp.ExpandRatio(2),
		},
		//{
		//	gp.DropThreshold(0.01),
		//	gp.PartialPivoting(1),
		//	gp.ColFillRatio(-1),
		//	gp.ColPerm(nil), // natural
		//},
		//{
		//	gp.DropThreshold(0.1),
		//	gp.PartialPivoting(1),
		//	gp.ColFillRatio(-1),
		//	gp.ColPerm(nil), // natural
		//},
		//{
		//	gp.DropThreshold(0),
		//	gp.PartialPivoting(1),
		//	gp.ColFillRatio(1),
		//	gp.ColPerm(nil), // natural
		//},
		//{
		//	gp.DropThreshold(0),
		//	gp.PartialPivoting(1),
		//	gp.ColFillRatio(12),
		//	gp.ColPerm(nil), // natural
		//},
		//
		//{
		//	gp.DropThreshold(0.01),
		//	gp.PartialPivoting(0),
		//	gp.ColFillRatio(-1),
		//	gp.ColPerm(nil), // natural
		//},
		//{
		//	gp.DropThreshold(0),
		//	gp.PartialPivoting(0),
		//	gp.ColFillRatio(-1),
		//	gp.ColPerm(nil), // natural
		//},
	} {
		b := matVec(n, rowind, colst, nzA, x0)

		lu, err := gp.Factor(n, rowind, colst, nzA, opts...)
		if err != nil {
			t.Fatalf("factor[%d]: %v", i, err)
		}

		err = gp.Solve(lu, [][]float64{b}, false)
		if err != nil {
			t.Fatalf("solve[%d]: %v", i, err)
		}

		const eps = 1e-10

		resid := residual(b)
		if resid > eps {
			t.Fatalf("resid[%d], expected < %v actual %v", i, eps, resid)
		}
	}

}

func matVec(n int, rowind, colst []int, nzA, x []float64) []float64 {
	y := make([]float64, n)
	for j := 0; j < n; j++ {
		start := colst[j]
		end := colst[j+1]

		for ii := int(start); ii < end; ii++ {
			i := rowind[ii]
			y[i] += nzA[ii] * x[j]
		}
	}
	return y
}

func residual(x []float64) float64 {
	norm := math.Inf(-1)
	for i := range x {
		//abs := math.Abs(x[i] - b[i])
		abs := math.Abs(x[i] - 1)
		if abs > norm {
			norm = abs
		}
	}
	return norm
}
