package fields

import (
	"github.com/paulgoleary/go-snark/amcl"
	"github.com/paulgoleary/go-snark/amcl/BLS381"
)

var Zero = BLS381.NewFPint(0)
var One = BLS381.NewFPint(1)

func VecEquals(a, b []*BLS381.FP) bool {
	if len(a) != len(b) {
		return false
	}
	for x, _ := range a {
		if !a[x].Equals(b[x]) {
			return false
		}
	}
	return true
}

// Fq is the Z field over modulus Q
type Fq struct {
	Q *BLS381.BIG // Q
}

// NewFq generates a new Fq
func NewFq(q *BLS381.BIG) Fq {
	return Fq{
	q,
	}
}

// Zero returns a Zero value on the Fq
func (fq Fq) Zero() *BLS381.FP {
	return BLS381.NewFPint(0)
}

// One returns a One value on the Fq
func (fq Fq) One() *BLS381.FP {
	return BLS381.NewFPint(1)
}

// Add performs an addition on the Fq
func (fq Fq) Add(a, b *BLS381.FP) *BLS381.FP {
	r := BLS381.NewFPcopy(a)
	r.Add(b)
	return r
}

// Double performs a doubling on the Fq
func (fq Fq) Double(a *BLS381.FP) *BLS381.FP {
	r := BLS381.NewFPcopy(a)
	r.Add(a)
	return r
}

// Sub performs a subtraction on the Fq
func (fq Fq) Sub(a, b *BLS381.FP) *BLS381.FP {
	r := BLS381.NewFPcopy(a)
	r.Sub(b)
	return r
}

// Neg performs a negation on the Fq
func (fq Fq) Neg(a *BLS381.FP) *BLS381.FP {
	r := BLS381.NewFPcopy(a)
	r.Neg()
	return r
}

// Mul performs a multiplication on the Fq
func (fq Fq) Mul(a, b *BLS381.FP) *BLS381.FP {
	r := BLS381.NewFPcopy(a)
	r.Mul(b)
	return r
}

func (fq Fq) MulScalar(base, e *BLS381.FP) *BLS381.FP {
	return fq.Mul(base, e)
}

// Inverse returns the inverse on the Fq
func (fq Fq) Inverse(a *BLS381.FP) *BLS381.FP {
	r := BLS381.NewFPcopy(a)
	r.Inverse()
	return r
}

// Div performs the division over the finite field
func (fq Fq) Div(a, b *BLS381.FP) *BLS381.FP {
	return fq.Mul(a, fq.Inverse(b))
}

// Square performs a square operation on the Fq
func (fq Fq) Square(a *BLS381.FP) *BLS381.FP {
	m := BLS381.NewFPcopy(a)
	m.Mul(a)
	return m
}

// Exp performs the exponential over Fq
func (fq Fq) Exp(base *BLS381.FP, e *BLS381.BIG) *BLS381.FP {
	res := BLS381.NewFPcopy(base)
	res.Pow(e)
	return res
}

var amclRand = amcl.NewRAND()

func (fq Fq) Rand() (*BLS381.FP, error) {
	return BLS381.NewFPbig(BLS381.Randomnum(fq.Q, amclRand)), nil
}

func (fq Fq) IsZero(a *BLS381.FP) bool {
	return a.Equals(Zero)
}

func (fq Fq) Copy(a *BLS381.FP) *BLS381.FP {
	return BLS381.NewFPcopy(a)
}

func (fq Fq) Equal(a, b *BLS381.FP) bool {
	return a.Equals(b)
}
