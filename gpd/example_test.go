// Copyright 1988 John Gilbert and Tim Peierls
// All rights reserved.

package gpd_test

import (
	"fmt"

	"github.com/rwl/lufact/gpd"
)

// ExampleFactor factorizes a 10x10 matrix, given a column permutation
// vector, and solves for a single right-hand-side.
func ExampleFactor() {
	// A = [
	//	[2.10                               0.14 0.09     ]
	//	[     1.10           0.06                     0.03]
	//	[          1.70                               0.04]
	//	[               1.00           0.32 0.19 0.32 0.44]
	//	[     0.06           1.60                         ]
	//	[                         2.20                    ]
	//	[               0.32           1.90           0.43]
	//	[0.14           0.19                1.10 0.22     ]
	//	[0.09           0.32                0.22 2.40     ]
	//	[     0.03 0.04 0.44           0.43           3.20]
	// ]
	var (
		n      = 10
		arow   = []int{0, 7, 8, 1, 4, 9, 2, 9, 3, 6, 7, 8, 9, 1, 4, 5, 3, 6, 9, 0, 3, 7, 8, 0, 3, 7, 8, 1, 2, 3, 6, 9}
		acolst = []int{0, 3, 6, 8, 13, 15, 16, 19, 23, 27, 32}
		a      = []float64{2.1, 0.14, 0.09, 1.1, 0.06, 0.03, 1.7, 0.04, 1, 0.32, 0.19, 0.32, 0.44, 0.06, 1.6, 2.2, 0.32, 1.9, 0.43, 0.14, 0.19, 1.1, 0.22, 0.09, 0.32, 0.22, 2.4, 0.03, 0.04, 0.44, 0.43, 3.2}

		b = []float64{0.403, 0.28, 0.55, 1.504, 0.812, 1.32, 1.888, 1.168, 2.473, 3.695}

		colPerm []int
	)
	colPerm = []int{6, 5, 2, 4, 1, 9, 7, 8, 0, 3}

	lu, err := gpd.Factor(n, arow, acolst, a, gpd.ColPerm(colPerm))
	if err != nil {
		panic(err)
	}

	err = gpd.Solve(lu, [][]float64{b}, true)
	if err != nil {
		panic(err)
	}

	fmt.Print("x =")
	for i := 0; i < n; i++ {
		fmt.Printf(" %.1f", b[i])
	}

	// Output:
	//  x = 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
}
