package circuitcompiler

import (
	"errors"
	"math/big"
	"strconv"

	"github.com/paulgoleary/go-snark/r1csqap"
)

// Circuit is the data structure of the compiled circuit
type Circuit struct {
	NVars         int
	NPublic       int
	NSignals      int
	Inputs        []string
	Signals       []string
	PublicSignals []string
	Witness       []*big.Int
	Constraints   []Constraint
	R1CS          struct {
		A [][]*big.Int
		B [][]*big.Int
		C [][]*big.Int
	}
}

// Constraint is the data structure of a flat code operation
type Constraint struct {
	// v1 op v2 = out
	Op      string
	V1      string
	V2      string
	Out     string
	Literal string

	Inputs []string // in func declaration case
}

func indexInArray(arr []string, e string) int {
	for i, a := range arr {
		if a == e {
			return i
		}
	}
	return -1
}
func isValue(a string) (bool, int) {
	v, err := strconv.Atoi(a)
	if err != nil {
		return false, 0
	}
	return true, v
}
func insertVar(arr []*big.Int, signals []string, v string, used map[string]bool) ([]*big.Int, map[string]bool) {
	isVal, value := isValue(v)
	valueBigInt := big.NewInt(int64(value))
	if isVal {
		arr[0] = new(big.Int).Add(arr[0], valueBigInt)
	} else {
		if !used[v] {
			panic(errors.New("using variable before it's set"))
		}
		arr[indexInArray(signals, v)] = new(big.Int).Add(arr[indexInArray(signals, v)], big.NewInt(int64(1)))
	}
	return arr, used
}
func insertVarNeg(arr []*big.Int, signals []string, v string, used map[string]bool) ([]*big.Int, map[string]bool) {
	isVal, value := isValue(v)
	valueBigInt := big.NewInt(int64(value))
	if isVal {
		arr[0] = new(big.Int).Add(arr[0], valueBigInt)
	} else {
		if !used[v] {
			panic(errors.New("using variable before it's set"))
		}
		arr[indexInArray(signals, v)] = new(big.Int).Add(arr[indexInArray(signals, v)], big.NewInt(int64(-1)))
	}
	return arr, used
}

// GenerateR1CS generates the R1CS polynomials from the Circuit
func (circ *Circuit) GenerateR1CS() ([][]*big.Int, [][]*big.Int, [][]*big.Int) {
	// from flat code to R1CS

	var a [][]*big.Int
	var b [][]*big.Int
	var c [][]*big.Int

	used := make(map[string]bool)
	for _, constraint := range circ.Constraints {
		aConstraint := r1csqap.ArrayOfBigZeros(len(circ.Signals))
		bConstraint := r1csqap.ArrayOfBigZeros(len(circ.Signals))
		cConstraint := r1csqap.ArrayOfBigZeros(len(circ.Signals))

		// if existInArray(constraint.Out) {
		if used[constraint.Out] {
			panic(errors.New("out variable already used: " + constraint.Out))
		}
		used[constraint.Out] = true
		if constraint.Op == "in" {
			for i := 0; i < len(constraint.Inputs); i++ {
				aConstraint[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Add(aConstraint[indexInArray(circ.Signals, constraint.Out)], big.NewInt(int64(1)))
				aConstraint, used = insertVar(aConstraint, circ.Signals, constraint.Out, used)
				bConstraint[0] = big.NewInt(int64(1))
			}
			continue

		} else if constraint.Op == "+" {
			cConstraint[indexInArray(circ.Signals, constraint.Out)] = big.NewInt(int64(1))
			aConstraint, used = insertVar(aConstraint, circ.Signals, constraint.V1, used)
			aConstraint, used = insertVar(aConstraint, circ.Signals, constraint.V2, used)
			bConstraint[0] = big.NewInt(int64(1))
		} else if constraint.Op == "-" {
			cConstraint[indexInArray(circ.Signals, constraint.Out)] = big.NewInt(int64(1))
			aConstraint, used = insertVarNeg(aConstraint, circ.Signals, constraint.V1, used)
			aConstraint, used = insertVarNeg(aConstraint, circ.Signals, constraint.V2, used)
			bConstraint[0] = big.NewInt(int64(1))
		} else if constraint.Op == "*" {
			cConstraint[indexInArray(circ.Signals, constraint.Out)] = big.NewInt(int64(1))
			aConstraint, used = insertVar(aConstraint, circ.Signals, constraint.V1, used)
			bConstraint, used = insertVar(bConstraint, circ.Signals, constraint.V2, used)
		} else if constraint.Op == "/" {
			cConstraint, used = insertVar(cConstraint, circ.Signals, constraint.V1, used)
			cConstraint[indexInArray(circ.Signals, constraint.Out)] = big.NewInt(int64(1))
			bConstraint, used = insertVar(bConstraint, circ.Signals, constraint.V2, used)
		}

		a = append(a, aConstraint)
		b = append(b, bConstraint)
		c = append(c, cConstraint)

	}
	return a, b, c
}

func grabVar(signals []string, w []*big.Int, vStr string) *big.Int {
	isVal, v := isValue(vStr)
	vBig := big.NewInt(int64(v))
	if isVal {
		return vBig
	} else {
		return w[indexInArray(signals, vStr)]
	}
}

// CalculateWitness calculates the Witness of a Circuit based on the given inputs
func (circ *Circuit) CalculateWitness(inputs []*big.Int) ([]*big.Int, error) {
	if len(inputs) != len(circ.Inputs) {
		return []*big.Int{}, errors.New("given inputs != circuit.Inputs")
	}
	w := r1csqap.ArrayOfBigZeros(len(circ.Signals))
	w[0] = big.NewInt(int64(1))
	for i, input := range inputs {
		w[i+2] = input
	}
	for _, constraint := range circ.Constraints {
		if constraint.Op == "in" {
		} else if constraint.Op == "+" {
			w[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Add(grabVar(circ.Signals, w, constraint.V1), grabVar(circ.Signals, w, constraint.V2))
		} else if constraint.Op == "-" {
			w[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Sub(grabVar(circ.Signals, w, constraint.V1), grabVar(circ.Signals, w, constraint.V2))
		} else if constraint.Op == "*" {
			w[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Mul(grabVar(circ.Signals, w, constraint.V1), grabVar(circ.Signals, w, constraint.V2))
		} else if constraint.Op == "/" {
			w[indexInArray(circ.Signals, constraint.Out)] = new(big.Int).Div(grabVar(circ.Signals, w, constraint.V1), grabVar(circ.Signals, w, constraint.V2))
		}
	}
	return w, nil
}
