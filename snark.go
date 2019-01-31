package snark

import (
	"fmt"
	"math/big"
	"os"

	"github.com/paulgoleary/go-snark/bn128"
	"github.com/paulgoleary/go-snark/circuitcompiler"
	"github.com/paulgoleary/go-snark/fields"
	"github.com/paulgoleary/go-snark/r1csqap"
)

// Setup is the data structure holding the Trusted Setup data. The Setup.Toxic sub struct must be destroyed after the GenerateTrustedSetup function is completed
type Setup struct {
	Toxic struct {
		T      *big.Int // trusted setup secret
		Ka     *big.Int // prover
		Kb     *big.Int // prover
		Kc     *big.Int // prover
		Kbeta  *big.Int
		Kgamma *big.Int
		RhoA   *big.Int
		RhoB   *big.Int
		RhoC   *big.Int
	}

	// public
	G1T [][3]*big.Int    // t encrypted in G1 curve
	G2T [][3][2]*big.Int // t encrypted in G2 curve
	Pk  struct {         // Proving Key pk:=(pkA, pkB, pkC, pkH)
		A  [][3]*big.Int
		B  [][3][2]*big.Int
		C  [][3]*big.Int
		Kp [][3]*big.Int
		Ap [][3]*big.Int
		Bp [][3]*big.Int
		Cp [][3]*big.Int
	}
	Vk struct {
		Vka   [3][2]*big.Int
		Vkb   [3]*big.Int
		Vkc   [3][2]*big.Int
		A     [][3]*big.Int
		G1Kbg [3]*big.Int    // g1 * Kbeta * Kgamma
		G2Kbg [3][2]*big.Int // g2 * Kbeta * Kgamma
		G2Kg  [3][2]*big.Int // g2 * Kgamma
		Vkz   [3][2]*big.Int
	}
}

// Proof contains the parameters to proof the zkSNARK
type Proof struct {
	PiA           [3]*big.Int
	PiAp          [3]*big.Int
	PiB           [3][2]*big.Int
	PiBp          [3]*big.Int
	PiC           [3]*big.Int
	PiCp          [3]*big.Int
	PiH           [3]*big.Int
	PiKp          [3]*big.Int
	PublicSignals []*big.Int
}

type utils struct {
	Bn  bn128.Bn128
	FqR fields.Fq
	PF  r1csqap.PolynomialField
}

// Utils is the data structure holding the BN128, FqR Finite Field over R, PolynomialField, that will be used inside the snarks operations
var Utils = prepareUtils()

func prepareUtils() utils {
	bn, err := bn128.NewBn128()
	if err != nil {
		panic(err)
	}
	// new Finite Field
	fqR := fields.NewFq(bn.R)
	// new Polynomial Field
	pf := r1csqap.NewPolynomialField(fqR)

	return utils{
		Bn:  bn,
		FqR: fqR,
		PF:  pf,
	}
}

// GenerateTrustedSetup generates the Trusted Setup from a compiled Circuit. The Setup.Toxic sub data structure must be destroyed
func GenerateTrustedSetup(witnessLength int, circuit circuitcompiler.Circuit, alphas, betas, gammas [][]*big.Int, zx []*big.Int) (Setup, error) {
	var setup Setup
	var err error
	// generate random t value
	setup.Toxic.T, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}

	// k for calculating pi' and Vk
	setup.Toxic.Ka, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}
	setup.Toxic.Kb, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}
	setup.Toxic.Kc, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}

	// generate Kβ (Kbeta) and Kγ (Kgamma)
	setup.Toxic.Kbeta, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}
	setup.Toxic.Kgamma, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}

	// generate ρ (Rho): ρA, ρB, ρC
	setup.Toxic.RhoA, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}
	setup.Toxic.RhoB, err = Utils.FqR.Rand()
	if err != nil {
		return Setup{}, err
	}
	setup.Toxic.RhoC = Utils.FqR.Mul(setup.Toxic.RhoA, setup.Toxic.RhoB)

	// encrypt t values with curve generators
	var gt1 [][3]*big.Int
	var gt2 [][3][2]*big.Int
	for i := 0; i < witnessLength; i++ {
		tPow := Utils.FqR.Exp(setup.Toxic.T, big.NewInt(int64(i)))
		tEncr1 := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, tPow)
		gt1 = append(gt1, tEncr1)
		tEncr2 := Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, tPow)
		gt2 = append(gt2, tEncr2)
	}
	// gt1: g1, g1*t, g1*t^2, g1*t^3, ...
	// gt2: g2, g2*t, g2*t^2, ...
	setup.G1T = gt1
	setup.G2T = gt2

	setup.Vk.Vka = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Ka)
	setup.Vk.Vkb = Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, setup.Toxic.Kb)
	setup.Vk.Vkc = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Kc)

	/*
		Verification keys:
		- Vk_betagamma1: setup.G1Kbg = g1 * Kbeta*Kgamma
		- Vk_betagamma2: setup.G2Kbg = g2 * Kbeta*Kgamma
		- Vk_gamma: setup.G2Kg = g2 * Kgamma
	*/
	kbg := Utils.FqR.Mul(setup.Toxic.Kbeta, setup.Toxic.Kgamma)
	setup.Vk.G1Kbg = Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, kbg)
	setup.Vk.G2Kbg = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, kbg)
	setup.Vk.G2Kg = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, setup.Toxic.Kgamma)

	// for i := 0; i < circuit.NSignals; i++ {
	for i := 0; i < circuit.NVars; i++ {
		at := Utils.PF.Eval(alphas[i], setup.Toxic.T)
		a := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, at)
		setup.Pk.A = append(setup.Pk.A, a)
		if i <= circuit.NPublic {
			setup.Vk.A = append(setup.Vk.A, a)
		}

		bt := Utils.PF.Eval(betas[i], setup.Toxic.T)
		bg1 := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, bt)
		bg2 := Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, bt)
		setup.Pk.B = append(setup.Pk.B, bg2)

		ct := Utils.PF.Eval(gammas[i], setup.Toxic.T)
		c := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, ct)
		setup.Pk.C = append(setup.Pk.C, c)

		kt := Utils.FqR.Add(Utils.FqR.Add(at, bt), ct)
		k := Utils.Bn.G1.Affine(Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, kt))

		ktest := Utils.Bn.G1.Affine(Utils.Bn.G1.Add(Utils.Bn.G1.Add(a, bg1), c))
		if !Utils.Bn.Fq2.Equal(k, ktest) {
			os.Exit(1)
			return setup, err
		}

		setup.Pk.Ap = append(setup.Pk.Ap, Utils.Bn.G1.MulScalar(a, setup.Toxic.Ka))
		setup.Pk.Bp = append(setup.Pk.Bp, Utils.Bn.G1.MulScalar(bg1, setup.Toxic.Kb))
		setup.Pk.Cp = append(setup.Pk.Cp, Utils.Bn.G1.MulScalar(c, setup.Toxic.Kc))
		k_ := Utils.Bn.G1.MulScalar(Utils.Bn.G1.G, kt)
		setup.Pk.Kp = append(setup.Pk.Kp, Utils.Bn.G1.MulScalar(k_, setup.Toxic.Kbeta))
	}
	setup.Vk.Vkz = Utils.Bn.G2.MulScalar(Utils.Bn.G2.G, Utils.PF.Eval(zx, setup.Toxic.T))

	return setup, nil
}

// GenerateProofs generates all the parameters to proof the zkSNARK from the Circuit, Setup and the Witness
func GenerateProofs(circuit circuitcompiler.Circuit, setup Setup, hx []*big.Int, w []*big.Int) (Proof, error) {
	var proof Proof
	proof.PiA = [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}
	proof.PiAp = [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}
	proof.PiB = Utils.Bn.Fq6.Zero()
	proof.PiBp = [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}
	proof.PiC = [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}
	proof.PiCp = [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}
	proof.PiH = [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}
	proof.PiKp = [3]*big.Int{Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero(), Utils.Bn.G1.F.Zero()}

	for i := circuit.NPublic + 1; i < circuit.NVars; i++ {
		proof.PiA = Utils.Bn.G1.Add(proof.PiA, Utils.Bn.G1.MulScalar(setup.Pk.A[i], w[i]))
		proof.PiAp = Utils.Bn.G1.Add(proof.PiAp, Utils.Bn.G1.MulScalar(setup.Pk.Ap[i], w[i]))
	}

	for i := 0; i < circuit.NVars; i++ {
		proof.PiB = Utils.Bn.G2.Add(proof.PiB, Utils.Bn.G2.MulScalar(setup.Pk.B[i], w[i]))
		proof.PiBp = Utils.Bn.G1.Add(proof.PiBp, Utils.Bn.G1.MulScalar(setup.Pk.Bp[i], w[i]))

		proof.PiC = Utils.Bn.G1.Add(proof.PiC, Utils.Bn.G1.MulScalar(setup.Pk.C[i], w[i]))
		proof.PiCp = Utils.Bn.G1.Add(proof.PiCp, Utils.Bn.G1.MulScalar(setup.Pk.Cp[i], w[i]))

		proof.PiKp = Utils.Bn.G1.Add(proof.PiKp, Utils.Bn.G1.MulScalar(setup.Pk.Kp[i], w[i]))
	}

	for i := 0; i < len(hx); i++ {
		proof.PiH = Utils.Bn.G1.Add(proof.PiH, Utils.Bn.G1.MulScalar(setup.G1T[i], hx[i]))
	}
	proof.PublicSignals = w[1:2] // out signal

	return proof, nil
}

// VerifyProof verifies over the BN128 the Pairings of the Proof
func VerifyProof(circuit circuitcompiler.Circuit, setup Setup, proof Proof, printVer bool) bool {
	// e(piA, Va) == e(piA', g2)
	pairingPiaVa := Utils.Bn.Pairing(proof.PiA, setup.Vk.Vka)
	pairingPiapG2 := Utils.Bn.Pairing(proof.PiAp, Utils.Bn.G2.G)
	if !Utils.Bn.Fq12.Equal(pairingPiaVa, pairingPiapG2) {
		return false
	}
	if printVer {
		fmt.Println("✓ e(piA, Va) == e(piA', g2), valid knowledge commitment for A")
	}

	// e(Vb, piB) == e(piB', g2)
	pairingVbPib := Utils.Bn.Pairing(setup.Vk.Vkb, proof.PiB)
	pairingPibpG2 := Utils.Bn.Pairing(proof.PiBp, Utils.Bn.G2.G)
	if !Utils.Bn.Fq12.Equal(pairingVbPib, pairingPibpG2) {
		return false
	}
	if printVer {
		fmt.Println("✓ e(Vb, piB) == e(piB', g2), valid knowledge commitment for B")
	}

	// e(piC, Vc) == e(piC', g2)
	pairingPicVc := Utils.Bn.Pairing(proof.PiC, setup.Vk.Vkc)
	pairingPicpG2 := Utils.Bn.Pairing(proof.PiCp, Utils.Bn.G2.G)
	if !Utils.Bn.Fq12.Equal(pairingPicVc, pairingPicpG2) {
		return false
	}
	if printVer {
		fmt.Println("✓ e(piC, Vc) == e(piC', g2), valid knowledge commitment for C")
	}

	// Vkx, to then calculate Vkx+piA
	vkxpia := setup.Vk.A[0]
	for i := 0; i < len(proof.PublicSignals); i++ {
		vkxpia = Utils.Bn.G1.Add(vkxpia, Utils.Bn.G1.MulScalar(setup.Vk.A[i+1], proof.PublicSignals[i]))
	}

	// e(Vkx+piA, piB) == e(piH, Vkz) * e(piC, g2)
	if !Utils.Bn.Fq12.Equal(
		Utils.Bn.Pairing(Utils.Bn.G1.Add(vkxpia, proof.PiA), proof.PiB),
		Utils.Bn.Fq12.Mul(
			Utils.Bn.Pairing(proof.PiH, setup.Vk.Vkz),
			Utils.Bn.Pairing(proof.PiC, Utils.Bn.G2.G))) {
		return false
	}
	if printVer {
		fmt.Println("✓ e(Vkx+piA, piB) == e(piH, Vkz) * e(piC, g2), QAP disibility checked")
	}

	// e(Vkx+piA+piC, g2KbetaKgamma) * e(g1KbetaKgamma, piB)
	// == e(piK, g2Kgamma)
	piApiC := Utils.Bn.G1.Add(Utils.Bn.G1.Add(vkxpia, proof.PiA), proof.PiC)
	pairingPiACG2Kbg := Utils.Bn.Pairing(piApiC, setup.Vk.G2Kbg)
	pairingG1KbgPiB := Utils.Bn.Pairing(setup.Vk.G1Kbg, proof.PiB)
	pairingL := Utils.Bn.Fq12.Mul(pairingPiACG2Kbg, pairingG1KbgPiB)
	pairingR := Utils.Bn.Pairing(proof.PiKp, setup.Vk.G2Kg)
	if !Utils.Bn.Fq12.Equal(pairingL, pairingR) {
		return false
	}
	if printVer {
		fmt.Println("✓ e(Vkx+piA+piC, g2KbetaKgamma) * e(g1KbetaKgamma, piB) == e(piK, g2Kgamma)")
	}

	return true
}
