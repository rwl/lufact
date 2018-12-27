// Code generated with gpgen. DO NOT EDIT.

// Copyright 1988 John Gilbert and Tim Peierls
// All rights reserved.

// Package gp provides sparse LU factorization with partial pivoting.
//
// The algorithm is described in "Sparse Partial Pivoting in Time Proportional
// to Arithmetic Operations" by John R. Gilbert and Tim Peierls.
//
//  @article{Gilbert1988,
//    doi = {10.1137/0909058},
//    url = {https://doi.org/10.1137/0909058},
//    year  = {1988},
//    month = {sep},
//    publisher = {Society for Industrial {\&} Applied Mathematics ({SIAM})},
//    volume = {9},
//    number = {5},
//    pages = {862--874},
//    author = {John R. Gilbert and Tim Peierls},
//    title = {Sparse Partial Pivoting in Time Proportional to Arithmetic Operations},
//    journal = {SIAM Journal on Scientific and Statistical Computing}
//  }
//
// This package is translated from the gp FORTRAN code distributed in
// Sivan Toledo's work on incomplete-factorization, from PARC in the
// early 1990s, as published in the ILU package on Netlib:
//
// http://www.netlib.org/linalg/ilu.tgz
//
// This source code is distributed, with the kind permission of John Gilbert
// and Tim Peierls, under a 3-clause BSD license.
package gpz
