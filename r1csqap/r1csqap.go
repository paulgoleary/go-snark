package r1csqap

import (
	"fmt"
	"github.com/paulgoleary/go-snark/amcl/BLS381"
	"github.com/paulgoleary/go-snark/fields"
)

// Transpose transposes the *big.Int matrix
func Transpose(matrix [][]*BLS381.FP) [][]*BLS381.FP {
	var r [][]*BLS381.FP
	for i := 0; i < len(matrix[0]); i++ {
		var row []*BLS381.FP
		for j := 0; j < len(matrix); j++ {
			row = append(row, matrix[j][i])
		}
		r = append(r, row)
	}
	return r
}

// ArrayOfBigZeros creates a *BLS381.FP array with n elements to zero
func ArrayOfBigZeros(num int) []*BLS381.FP {
	blsZero := BLS381.NewFPint(0)
	var r []*BLS381.FP
	for i := 0; i < num; i++ {
		r = append(r, blsZero)
	}
	return r
}

func BigArraysEqual(a, b []*BLS381.FP) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		if !a[i].Equals(b[i]) {
			return false
		}
	}
	return true
}

// PolynomialField is the Polynomial over a Finite Field where the polynomial operations are performed
type PolynomialField struct {
	F fields.Fq
}

// NewPolynomialField creates a new PolynomialField with the given FiniteField
func NewPolynomialField(f fields.Fq) PolynomialField {
	return PolynomialField{
		f,
	}
}

// Mul multiplies two polinomials over the Finite Field
func (pf PolynomialField) Mul(a, b []*BLS381.FP) []*BLS381.FP {
	r := ArrayOfBigZeros(len(a) + len(b) - 1)
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(b); j++ {
			r[i+j] = pf.F.Add(
				r[i+j],
				pf.F.Mul(a[i], b[j]))
		}
	}
	return r
}

// Div divides two polinomials over the Finite Field, returning the result and the remainder
func (pf PolynomialField) Div(a, b []*BLS381.FP) ([]*BLS381.FP, []*BLS381.FP) {
	// https://en.wikipedia.org/wiki/Division_algorithm
	r := ArrayOfBigZeros(len(a) - len(b) + 1)
	rem := a
	for len(rem) >= len(b) {
		l := pf.F.Div(rem[len(rem)-1], b[len(b)-1])
		pos := len(rem) - len(b)
		r[pos] = l
		aux := ArrayOfBigZeros(pos)
		aux1 := append(aux, l)
		aux2 := pf.Sub(rem, pf.Mul(b, aux1))
		rem = aux2[:len(aux2)-1]
	}
	return r, rem
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// Add adds two polynomials over the Finite Field
func (pf PolynomialField) Add(a, b []*BLS381.FP) []*BLS381.FP {
	r := ArrayOfBigZeros(max(len(a), len(b)))
	for i := 0; i < len(a); i++ {
		r[i] = pf.F.Add(r[i], a[i])
	}
	for i := 0; i < len(b); i++ {
		r[i] = pf.F.Add(r[i], b[i])
	}
	return r
}

// Sub subtracts two polinomials over the Finite Field
func (pf PolynomialField) Sub(a, b []*BLS381.FP) []*BLS381.FP {
	r := ArrayOfBigZeros(max(len(a), len(b)))
	for i := 0; i < len(a); i++ {
		r[i] = pf.F.Add(r[i], a[i])
	}
	for i := 0; i < len(b); i++ {
		r[i] = pf.F.Sub(r[i], b[i])
	}
	return r
}

// Eval evaluates the polinomial over the Finite Field at the given value x
func (pf PolynomialField) Eval(v []*BLS381.FP, x *BLS381.FP) *BLS381.FP {
	r := BLS381.NewFPint(0)
	for i := 0; i < len(v); i++ {
		xi := pf.F.Exp(x, BLS381.NewBIGint(i))
		elem := pf.F.Mul(v[i], xi)
		r = pf.F.Add(r, elem)
		println(fmt.Sprintf("%v", r))
	}
	return r
}

// NewPolZeroAt generates a new polynomial that has value zero at the given value
func (pf PolynomialField) NewPolZeroAt(pointPos, totalPoints int, height *BLS381.FP) []*BLS381.FP {
	fac := 1
	for i := 1; i < totalPoints+1; i++ {
		if i != pointPos {
			fac = fac * (pointPos - i)
		}
	}
	facBls := BLS381.NewFPint(fac)
	hf := pf.F.Div(height, facBls)
	r := []*BLS381.FP{hf}
	for i := 1; i < totalPoints+1; i++ {
		if i != pointPos {
			ineg := BLS381.NewFPint(-i)
			b1 := BLS381.NewFPint(1)
			r = pf.Mul(r, []*BLS381.FP{ineg, b1})
		}
	}
	return r
}

// LagrangeInterpolation performs the Lagrange Interpolation / Lagrange Polynomials operation
func (pf PolynomialField) LagrangeInterpolation(v []*BLS381.FP) []*BLS381.FP {
	// https://en.wikipedia.org/wiki/Lagrange_polynomial
	var r []*BLS381.FP
	for i := 0; i < len(v); i++ {
		r = pf.Add(r, pf.NewPolZeroAt(i+1, len(v), v[i]))
	}
	//
	return r
}

// R1CSToQAP converts the R1CS values to the QAP values
func (pf PolynomialField) R1CSToQAP(a, b, c [][]*BLS381.FP) ([][]*BLS381.FP, [][]*BLS381.FP, [][]*BLS381.FP, []*BLS381.FP) {
	aT := Transpose(a)
	bT := Transpose(b)
	cT := Transpose(c)
	var alphas [][]*BLS381.FP
	for i := 0; i < len(aT); i++ {
		alphas = append(alphas, pf.LagrangeInterpolation(aT[i]))
	}
	var betas [][]*BLS381.FP
	for i := 0; i < len(bT); i++ {
		betas = append(betas, pf.LagrangeInterpolation(bT[i]))
	}
	var gammas [][]*BLS381.FP
	for i := 0; i < len(cT); i++ {
		gammas = append(gammas, pf.LagrangeInterpolation(cT[i]))
	}
	z := []*BLS381.FP{BLS381.NewFPint(1)}
	for i := 1; i < len(aT[0])+1; i++ {
		ineg := BLS381.NewFPint(-i)
		b1 := BLS381.NewFPint(1)
		z = pf.Mul(z, []*BLS381.FP{ineg, b1})
	}
	return alphas, betas, gammas, z
}

// CombinePolynomials combine the given polynomials arrays into one, also returns the P(x)
func (pf PolynomialField) CombinePolynomials(r []*BLS381.FP, ap, bp, cp [][]*BLS381.FP) ([]*BLS381.FP, []*BLS381.FP, []*BLS381.FP, []*BLS381.FP) {
	var alpha []*BLS381.FP
	for i := 0; i < len(r); i++ {
		m := pf.Mul([]*BLS381.FP{r[i]}, ap[i])
		alpha = pf.Add(alpha, m)
	}
	var beta []*BLS381.FP
	for i := 0; i < len(r); i++ {
		m := pf.Mul([]*BLS381.FP{r[i]}, bp[i])
		beta = pf.Add(beta, m)
	}
	var gamma []*BLS381.FP
	for i := 0; i < len(r); i++ {
		m := pf.Mul([]*BLS381.FP{r[i]}, cp[i])
		gamma = pf.Add(gamma, m)
	}

	px := pf.Sub(pf.Mul(alpha, beta), gamma)
	return alpha, beta, gamma, px
}

// DivisorPolynomial returns the divisor polynomial given two polynomials
func (pf PolynomialField) DivisorPolynomial(px, z []*BLS381.FP) []*BLS381.FP {
	quo, _ := pf.Div(px, z)
	return quo
}
