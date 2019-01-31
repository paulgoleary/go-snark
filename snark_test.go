package snark

import (
	"fmt"
	"math/big"
	"strings"
	"testing"
	"time"

	"github.com/arnaucube/go-snark/circuitcompiler"
	"github.com/arnaucube/go-snark/r1csqap"
	"github.com/stretchr/testify/assert"
)

func hornerPolyEval(poly []*big.Int, x, m *big.Int) *big.Int {

	res := big.NewInt(0)
	for i := len(poly) - 1; i >= 0; i-- {
		res.Mul(res, x)
		res.Mod(res, m)
		res.Add(res, poly[i])
		res.Mod(res, m)
	}

	return res
}

func cmpVects(a, b []*big.Int) bool {
	if len(a) != len(b) {
		return false
	}
	for x, _ := range a {
		if a[x].Cmp(b[x]) != 0 {
			return false
		}
	}
	return true
}

func checkR1CSToQAP(r1csIn [][]*big.Int, qapIn [][]*big.Int, m *big.Int) bool {

	for x, polyR := range r1csIn {
		res := make([]*big.Int, len(polyR))
		biX := big.NewInt(int64(x + 1))
		for y, polyQ := range qapIn {
			res[y] = hornerPolyEval(polyQ, biX, m)
		}
		cmp := cmpVects(res, polyR)
		if !cmp {
			return false
		}
	}
	return true
}

func TestDirectPoly(t *testing.T) {

	testPoly := make([]*big.Int, 4)
	testPoly[0], _ = new(big.Int).SetString("3", 10)
	testPoly[1], _ = new(big.Int).SetString("3648040478639879203707734290876212514758060733402672390616367364429301415931", 10)
	testPoly[2], _ = new(big.Int).SetString("10944121435919637611123202872628637544274182200208017171849102093287904247811", 10)
	testPoly[3], _ = new(big.Int).SetString("7296080957279758407415468581752425029516121466805344781232734728858602831872", 10)

	testX := big.NewInt(2)
	testEval := Utils.PF.Eval(testPoly, testX)
	testEval2 := hornerPolyEval(testPoly, testX, Utils.PF.F.Q)
	assert.Equal(t, testEval, testEval2)
}

func TestZkFromFlatCircuitCode(t *testing.T) {

	// compile circuit and get the R1CS
	flatCode := `
	func test(x):
		aux = x*x
		y = aux*x
		z = x + y
		out = z + 5
	`
	fmt.Print("\nflat code of the circuit:")
	fmt.Println(flatCode)

	// parse the code
	parser := circuitcompiler.NewParser(strings.NewReader(flatCode))
	circuit, err := parser.Parse()
	assert.Nil(t, err)
	fmt.Println("\ncircuit data:", circuit)

	b3 := big.NewInt(int64(3))
	inputs := []*big.Int{b3}
	// witness
	witness, err := circuit.CalculateWitness(inputs)
	assert.Nil(t, err)
	fmt.Println("\nwitness", witness)

	// flat code to R1CS
	fmt.Println("\ngenerating R1CS from flat code")
	a, b, c := circuit.GenerateR1CS()
	fmt.Println("\nR1CS:")
	fmt.Println("a:", a)
	fmt.Println("b:", b)
	fmt.Println("c:", c)

	// R1CS to QAP
	alphas, betas, gammas, zx := Utils.PF.R1CSToQAP(a, b, c)
	fmt.Println("qap")
	fmt.Println(alphas)
	fmt.Println(betas)
	fmt.Println(gammas)

	checkA := checkR1CSToQAP(a, alphas, Utils.PF.F.Q)
	assert.Equal(t, checkA, a)
	checkB := checkR1CSToQAP(b, betas, Utils.PF.F.Q)
	assert.Equal(t, checkB, b)
	checkC := checkR1CSToQAP(c, gammas, Utils.PF.F.Q)
	assert.Equal(t, checkC, gammas)

	ax, bx, cx, px := Utils.PF.CombinePolynomials(witness, alphas, betas, gammas)

	hx := Utils.PF.DivisorPolynomial(px, zx)

	// hx==px/zx so px==hx*zx
	assert.Equal(t, px, Utils.PF.Mul(hx, zx))

	// p(x) = a(x) * b(x) - c(x) == h(x) * z(x)
	abc := Utils.PF.Sub(Utils.PF.Mul(ax, bx), cx)
	assert.Equal(t, abc, px)
	hz := Utils.PF.Mul(hx, zx)
	assert.Equal(t, abc, hz)

	div, rem := Utils.PF.Div(px, zx)
	assert.Equal(t, hx, div)
	assert.Equal(t, rem, r1csqap.ArrayOfBigZeros(4))

	// calculate trusted setup
	setup, err := GenerateTrustedSetup(len(witness), *circuit, alphas, betas, gammas, zx)
	assert.Nil(t, err)
	fmt.Println("\nt:", setup.Toxic.T)

	// piA = g1 * A(t), piB = g2 * B(t), piC = g1 * C(t), piH = g1 * H(t)
	proof, err := GenerateProofs(*circuit, setup, hx, witness)
	assert.Nil(t, err)

	fmt.Println("\n proofs:")
	fmt.Println(proof)
	fmt.Println("public signals:", proof.PublicSignals)
	before := time.Now()
	assert.True(t, VerifyProof(*circuit, setup, proof, true))
	fmt.Println("verify proof time elapsed:", time.Since(before))
}

func TestZkFromHardcodedR1CS(t *testing.T) {
	b0 := big.NewInt(int64(0))
	b1 := big.NewInt(int64(1))
	b3 := big.NewInt(int64(3))
	b5 := big.NewInt(int64(5))
	b9 := big.NewInt(int64(9))
	b27 := big.NewInt(int64(27))
	b30 := big.NewInt(int64(30))
	b35 := big.NewInt(int64(35))
	a := [][]*big.Int{
		[]*big.Int{b0, b0, b1, b0, b0, b0},
		[]*big.Int{b0, b0, b0, b1, b0, b0},
		[]*big.Int{b0, b0, b1, b0, b1, b0},
		[]*big.Int{b5, b0, b0, b0, b0, b1},
	}
	b := [][]*big.Int{
		[]*big.Int{b0, b0, b1, b0, b0, b0},
		[]*big.Int{b0, b0, b1, b0, b0, b0},
		[]*big.Int{b1, b0, b0, b0, b0, b0},
		[]*big.Int{b1, b0, b0, b0, b0, b0},
	}
	c := [][]*big.Int{
		[]*big.Int{b0, b0, b0, b1, b0, b0},
		[]*big.Int{b0, b0, b0, b0, b1, b0},
		[]*big.Int{b0, b0, b0, b0, b0, b1},
		[]*big.Int{b0, b1, b0, b0, b0, b0},
	}
	alphas, betas, gammas, zx := Utils.PF.R1CSToQAP(a, b, c)

	// wittness = 1, 35, 3, 9, 27, 30
	w := []*big.Int{b1, b35, b3, b9, b27, b30}
	circuit := circuitcompiler.Circuit{
		NVars:    6,
		NPublic:  1,
		NSignals: len(w),
	}
	ax, bx, cx, px := Utils.PF.CombinePolynomials(w, alphas, betas, gammas)

	hx := Utils.PF.DivisorPolynomial(px, zx)

	// hx==px/zx so px==hx*zx
	assert.Equal(t, px, Utils.PF.Mul(hx, zx))

	// p(x) = a(x) * b(x) - c(x) == h(x) * z(x)
	abc := Utils.PF.Sub(Utils.PF.Mul(ax, bx), cx)
	assert.Equal(t, abc, px)
	hz := Utils.PF.Mul(hx, zx)
	assert.Equal(t, abc, hz)

	div, rem := Utils.PF.Div(px, zx)
	assert.Equal(t, hx, div)
	assert.Equal(t, rem, r1csqap.ArrayOfBigZeros(4))

	// calculate trusted setup
	setup, err := GenerateTrustedSetup(len(w), circuit, alphas, betas, gammas, zx)
	assert.Nil(t, err)

	// piA = g1 * A(t), piB = g2 * B(t), piC = g1 * C(t), piH = g1 * H(t)
	proof, err := GenerateProofs(circuit, setup, hx, w)
	assert.Nil(t, err)

	assert.True(t, VerifyProof(circuit, setup, proof, true))
}

func TestZkMultiplication(t *testing.T) {

	// compile circuit and get the R1CS
	flatCode := `
	func test(a, b):
		out = a * b
	`

	// parse the code
	parser := circuitcompiler.NewParser(strings.NewReader(flatCode))
	circuit, err := parser.Parse()
	assert.Nil(t, err)

	b3 := big.NewInt(int64(3))
	b4 := big.NewInt(int64(4))
	inputs := []*big.Int{b3, b4}
	// wittness
	w, err := circuit.CalculateWitness(inputs)
	assert.Nil(t, err)

	// flat code to R1CS
	a, b, c := circuit.GenerateR1CS()

	// R1CS to QAP
	alphas, betas, gammas, zx := Utils.PF.R1CSToQAP(a, b, c)

	ax, bx, cx, px := Utils.PF.CombinePolynomials(w, alphas, betas, gammas)

	hx := Utils.PF.DivisorPolynomial(px, zx)

	// hx==px/zx so px==hx*zx
	assert.Equal(t, px, Utils.PF.Mul(hx, zx))

	// p(x) = a(x) * b(x) - c(x) == h(x) * z(x)
	abc := Utils.PF.Sub(Utils.PF.Mul(ax, bx), cx)
	assert.Equal(t, abc, px)
	hz := Utils.PF.Mul(hx, zx)
	assert.Equal(t, abc, hz)

	div, rem := Utils.PF.Div(px, zx)
	assert.Equal(t, hx, div)
	assert.Equal(t, rem, r1csqap.ArrayOfBigZeros(1))

	// calculate trusted setup
	setup, err := GenerateTrustedSetup(len(w), *circuit, alphas, betas, gammas, zx)
	assert.Nil(t, err)

	// piA = g1 * A(t), piB = g2 * B(t), piC = g1 * C(t), piH = g1 * H(t)
	proof, err := GenerateProofs(*circuit, setup, hx, w)
	assert.Nil(t, err)

	assert.True(t, VerifyProof(*circuit, setup, proof, false))
}
