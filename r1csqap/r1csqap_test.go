package r1csqap

import (
	"fmt"
	"github.com/paulgoleary/go-snark/amcl/BLS381"
	"testing"

	"github.com/paulgoleary/go-snark/fields"
	"github.com/stretchr/testify/assert"
)

func TestTranspose(t *testing.T) {
	b0 := BLS381.NewFPint(0)
	b1 := BLS381.NewFPint(1)
	bFive := BLS381.NewFPint(5)
	a := [][]*BLS381.FP{
		[]*BLS381.FP{b0, b1, b0, b0, b0, b0},
		[]*BLS381.FP{b0, b0, b0, b1, b0, b0},
		[]*BLS381.FP{b0, b1, b0, b0, b1, b0},
		[]*BLS381.FP{bFive, b0, b0, b0, b0, b1},
	}
	aT := Transpose(a)
	assert.Equal(t, aT, [][]*BLS381.FP{
		[]*BLS381.FP{b0, b0, b0, bFive},
		[]*BLS381.FP{b1, b0, b1, b0},
		[]*BLS381.FP{b0, b0, b0, b0},
		[]*BLS381.FP{b0, b1, b0, b0},
		[]*BLS381.FP{b0, b0, b1, b0},
		[]*BLS381.FP{b0, b0, b0, b1},
	})
}

func TestPol(t *testing.T) {
	b0 := BLS381.NewFPint(0)
	b1 := BLS381.NewFPint(1)
	b3 := BLS381.NewFPint(3)
	b4 := BLS381.NewFPint(4)
	b5 := BLS381.NewFPint(5)
	b6 := BLS381.NewFPint(6)
	b16 := BLS381.NewFPint(16)

	a := []*BLS381.FP{b1, b0, b5}
	b := []*BLS381.FP{b3, b0, b1}

	// new Finite Field
	r := BLS381.NewBIGints(BLS381.Modulus);
	f := fields.NewFq(r)

	// new Polynomial Field
	pf := NewPolynomialField(f)

	// polynomial multiplication
	o := pf.Mul(a, b)
	assert.True(t, fields.VecEquals(o, []*BLS381.FP{b3, b0, b16, b0, b5}))

	// polynomial addition
	o = pf.Add(a, b)
	assert.True(t, fields.VecEquals(o, []*BLS381.FP{b4, b0, b6}))

	// polynomial subtraction
	o1 := pf.Sub(a, b)
	o2 := pf.Sub(b, a)
	o = pf.Add(o1, o2)
	assert.True(t, b0.Equals(o[0])) // bytes.Equal(b0.Bytes(), o[0].Bytes()))
	assert.True(t, b0.Equals(o[1])) // bytes.Equal(b0.Bytes(), o[1].Bytes()))
	assert.True(t, b0.Equals(o[2])) // bytes.Equal(b0.Bytes(), o[2].Bytes()))

	c := []*BLS381.FP{b5, b6, b1}
	d := []*BLS381.FP{b1, b3}
	o = pf.Sub(c, d)
	assert.Equal(t, o, []*BLS381.FP{b4, b3, b1})

	// NewPolZeroAt
	o = pf.NewPolZeroAt(3, 4, b4)
	assert.Equal(t, pf.Eval(o, BLS381.NewFPint(3)), b4)
	o = pf.NewPolZeroAt(2, 4, b3)
	assert.Equal(t, pf.Eval(o, BLS381.NewFPint(2)), b3)
}

func TestLagrangeInterpolation(t *testing.T) {
	// new Finite Field
	r := BLS381.NewBIGints(BLS381.Modulus);
	f := fields.NewFq(r)
	// new Polynomial Field
	pf := NewPolynomialField(f)

	b0 := BLS381.NewFPint(0)
	b5 := BLS381.NewFPint(5)
	a := []*BLS381.FP{b0, b0, b0, b5}
	alpha := pf.LagrangeInterpolation(a)

	assert.Equal(t, pf.Eval(alpha, BLS381.NewFPint(4)), b5)
	aux := pf.Eval(alpha, BLS381.NewFPint(3))
	assert.Equal(t, aux, fields.Zero)

}

func TestR1CSToQAP(t *testing.T) {
	// new Finite Field
	r := BLS381.NewBIGints(BLS381.Modulus);
	f := fields.NewFq(r)
	// new Polynomial Field
	pf := NewPolynomialField(f)

	b0 := BLS381.NewFPint(0)
	b1 := BLS381.NewFPint(1)
	b3 := BLS381.NewFPint(3)
	b5 := BLS381.NewFPint(5)
	b9 := BLS381.NewFPint(9)
	b27 := BLS381.NewFPint(27)
	b30 := BLS381.NewFPint(30)
	b35 := BLS381.NewFPint(35)
	a := [][]*BLS381.FP{
		[]*BLS381.FP{b0, b1, b0, b0, b0, b0},
		[]*BLS381.FP{b0, b0, b0, b1, b0, b0},
		[]*BLS381.FP{b0, b1, b0, b0, b1, b0},
		[]*BLS381.FP{b5, b0, b0, b0, b0, b1},
	}
	b := [][]*BLS381.FP{
		[]*BLS381.FP{b0, b1, b0, b0, b0, b0},
		[]*BLS381.FP{b0, b1, b0, b0, b0, b0},
		[]*BLS381.FP{b1, b0, b0, b0, b0, b0},
		[]*BLS381.FP{b1, b0, b0, b0, b0, b0},
	}
	c := [][]*BLS381.FP{
		[]*BLS381.FP{b0, b0, b0, b1, b0, b0},
		[]*BLS381.FP{b0, b0, b0, b0, b1, b0},
		[]*BLS381.FP{b0, b0, b0, b0, b0, b1},
		[]*BLS381.FP{b0, b0, b1, b0, b0, b0},
	}
	alphas, betas, gammas, zx := pf.R1CSToQAP(a, b, c)
	fmt.Println(alphas)
	fmt.Println(betas)
	fmt.Println(gammas)
	fmt.Print("Z(x): ")
	fmt.Println(zx)

	w := []*BLS381.FP{b1, b3, b35, b9, b27, b30}
	ax, bx, cx, px := pf.CombinePolynomials(w, alphas, betas, gammas)
	fmt.Println(ax)
	fmt.Println(bx)
	fmt.Println(cx)
	fmt.Println(px)

	hx := pf.DivisorPolynomial(px, zx)
	fmt.Println(hx)

	// hx==px/zx so px==hx*zx
	assert.Equal(t, px, pf.Mul(hx, zx))

	// p(x) = a(x) * b(x) - c(x) == h(x) * z(x)
	abc := pf.Sub(pf.Mul(ax, bx), cx)
	assert.Equal(t, abc, px)
	hz := pf.Mul(hx, zx)
	assert.Equal(t, abc, hz)

}
