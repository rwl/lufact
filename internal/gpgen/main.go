// Copyright 2018 Richard Lincoln. All rights reserved.

package main

import (
	"bytes"
	"flag"
	"fmt"
	"go/format"
	"io/ioutil"
	"os"
	"path/filepath"
	"text/template"
)

var (
	tmplDir = filepath.Join(os.Getenv("GOPATH"), "src",
		"github.com", "rwlincoln", "lufact", "internal", "gpgen", "templates")
	outDir = filepath.Join(os.Getenv("GOPATH"), "src", "github.com", "rwlincoln", "lufact")

	formatOutput = flag.Bool("fmt", true, "format generated files")

	files = []string{
		"doc",
		"factor",
		"gp",
		"lsolve",
		"lucomp",
		"lucopy",
		"ludfs",
		//"lufact",
		"maxmatch",
		"usolve",
	}
)

type GPData struct {
	Package    string
	ScalarType string
}

func (GPData) Header() string {
	return `// Code generated with gpgen. DO NOT EDIT.

// Copyright 1988 John Gilbert and Tim Peierls
// All rights reserved.`
}

func execute() error {
	for _, filename := range files {
		tpath := fmt.Sprintf("%s/%s.tmpl", tmplDir, filename)
		tmpl, err := template.ParseFiles(tpath)
		if err != nil {
			return err
		}
		for _, t := range []GPData{
			{Package: "gpd", ScalarType: "float64"},
			{Package: "gpz", ScalarType: "complex128"},
		} {
			buf := new(bytes.Buffer)

			if err := tmpl.Execute(buf, t); err != nil {
				return err
			}

			var code []byte
			if *formatOutput {
				code, err = format.Source(buf.Bytes())
				if err != nil {
					return err
				}
			} else {
				code = buf.Bytes()
			}
			out := filepath.Join(outDir, t.Package, filename+".go")
			err := ioutil.WriteFile(out, code, 0644)
			if err != nil {
				return err
			}
		}
	}
	return nil
}

func main() {
	err := execute()
	if err != nil {
		fmt.Fprintln(os.Stderr, err.Error())
		os.Exit(2)
	}
}
