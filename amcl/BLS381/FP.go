/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/

/* Finite Field arithmetic */
/* CLINT mod p functions */

package BLS381

//import "fmt"

const NOT_SPECIAL int=0
const PSEUDO_MERSENNE int=1
const MONTGOMERY_FRIENDLY int=2
const GENERALISED_MERSENNE int=3

const MODBITS uint=381 /* Number of bits in Modulus */
const MOD8 uint=3  /* Modulus mod 8 */
const MODTYPE int=NOT_SPECIAL //NOT_SPECIAL

const FEXCESS int32=((int32(1)<<25)-1)
const OMASK Chunk= ((Chunk(-1))<<(MODBITS%BASEBITS))
const TBITS uint=MODBITS%BASEBITS // Number of active bits in top word 
const TMASK Chunk=(Chunk(1)<<TBITS)-1


type FP struct {
	x *BIG
	XES int32
}

/* Constructors */
func NewFPint(a int) *FP {
	F:=new(FP)
	F.x=NewBIGint(a)
	F.nres()
	return F
}

func NewFPbig(a *BIG) *FP {
	F:=new(FP)
	F.x=NewBIGcopy(a)
	F.nres()
	return F
}

func NewFPcopy(a *FP) *FP {
	F:=new(FP)
	F.x=NewBIGcopy(a.x)
	F.XES=a.XES
	return F
}

func (F *FP) toString() string {
	F.reduce()
	return F.redc().ToString()
}

/* convert to Montgomery n-residue form */
func (F *FP) nres() {
	if MODTYPE!=PSEUDO_MERSENNE && MODTYPE!=GENERALISED_MERSENNE {
		r:=NewBIGints(R2modp)
		d:=mul(F.x,r)
		F.x.copy(mod(d))
		F.XES=2
	} else {
		F.XES=1
	}
}

/* convert back to regular form */
func (F *FP) redc() *BIG {
	if MODTYPE!=PSEUDO_MERSENNE && MODTYPE!=GENERALISED_MERSENNE {
		d:=NewDBIGscopy(F.x)
		return mod(d)
	} else {
		r:=NewBIGcopy(F.x)
		return r
	}
}

/* reduce a DBIG to a BIG using the appropriate form of the modulus */

func mod(d *DBIG) *BIG {
	if MODTYPE==PSEUDO_MERSENNE {
		t:=d.split(MODBITS)
		b:=NewBIGdcopy(d)

		v:=t.pmul(int(MConst))

		t.add(b)
		t.norm()

		tw:=t.w[NLEN-1]
		t.w[NLEN-1]&=TMASK
		t.w[0]+=(MConst*((tw>>TBITS)+(v<<(BASEBITS-TBITS))))

		t.norm()
		return t	
	}
	if MODTYPE==MONTGOMERY_FRIENDLY {
		for i:=0;i<NLEN;i++ {
			top,bot:=muladd(d.w[i],MConst-1,d.w[i],d.w[NLEN+i-1])
			d.w[NLEN+i-1]=bot
			d.w[NLEN+i]+=top
		}
		b:=NewBIG()

		for i:=0;i<NLEN;i++ {
			b.w[i]=d.w[NLEN+i]
		}
		b.norm()
		return b		
	}

	if MODTYPE==GENERALISED_MERSENNE { // GoldiLocks only
		t:=d.split(MODBITS)
		b:=NewBIGdcopy(d)
		b.add(t);
		dd:=NewDBIGscopy(t)
		dd.shl(MODBITS/2)

		tt:=dd.split(MODBITS)
		lo:=NewBIGdcopy(dd)
		b.add(tt)
		b.add(lo)
		b.norm()
		tt.shl(MODBITS/2)
		b.add(tt)

		carry:=b.w[NLEN-1]>>TBITS
		b.w[NLEN-1]&=TMASK
		b.w[0]+=carry
			
		b.w[224/BASEBITS]+=carry<<(224%BASEBITS);
		b.norm()
		return b		
	}

	if MODTYPE==NOT_SPECIAL {
		md:=NewBIGints(Modulus)
		return monty(md,MConst,d) 
	}
	return NewBIG()
}

// find appoximation to quotient of a/m
// Out by at most 2.
// Note that MAXXES is bounded to be 2-bits less than half a word
func quo(n *BIG,m *BIG) int {
	var num  Chunk
	var den  Chunk
	hb:=uint(CHUNK)/2
	if TBITS < hb {
		sh:=hb-TBITS
		num=((n.w[NLEN-1]<<sh))|(n.w[NLEN-2]>>(BASEBITS-sh))
		den=((m.w[NLEN-1]<<sh))|(m.w[NLEN-2]>>(BASEBITS-sh))

	} else {
		num=n.w[NLEN-1]
		den=m.w[NLEN-1]
	}
	return int(num/(den+1))
}

/* reduce this mod Modulus */
func (F *FP) reduce() {
	m:=NewBIGints(Modulus)
	r:=NewBIGints(Modulus)
	var sb uint
	F.x.norm()

	if F.XES > 16 {
		q:=quo(F.x,m)
		carry:=r.pmul(q)
		r.w[NLEN-1]+=carry<<BASEBITS
		F.x.sub(r)
		F.x.norm()
		sb=2
	} else { 
		sb=logb2(uint32(F.XES-1))
	}

	m.fshl(sb)
	for sb>0 {
            sr:=ssn(r,F.x,m)
	    F.x.cmove(r,1-sr)
            sb -= 1
	}

	F.XES=1
}

/* test this=0? */
func (F *FP) iszilch() bool {
	W:=NewFPcopy(F)
	W.reduce()
	return W.x.iszilch()
}

/* copy from FP b */
func (F *FP) copy(b *FP ) {
	F.x.copy(b.x)
	F.XES=b.XES
}

/* set this=0 */
func (F *FP) zero() {
	F.x.zero()
	F.XES=1
}
	
/* set this=1 */
func (F *FP) one() {
	F.x.one(); F.nres()
}

/* normalise this */
func (F *FP) norm() {
	F.x.norm()
}

/* swap FPs depending on d */
func (F *FP) cswap(b *FP,d int) {
	c:=int32(d)
	c=^(c-1)
	t:=c&(F.XES^b.XES)
	F.XES^=t
	b.XES^=t
	F.x.cswap(b.x,d)
}

/* copy FPs depending on d */
func (F *FP) cmove(b *FP,d int) {
	F.x.cmove(b.x,d)
	c:=int32(-d)
	F.XES^=(F.XES^b.XES)&c
}

/* this*=b mod Modulus */
func (F *FP) Mul(b *FP) {

	if int64(F.XES)*int64(b.XES)>int64(FEXCESS) {F.reduce()}

	d:=mul(F.x,b.x)
	F.x.copy(mod(d))
	F.XES=2
}

/* this = -this mod Modulus */
func (F *FP) Neg() {
	m:=NewBIGints(Modulus)
	sb:=logb2(uint32(F.XES-1))

	m.fshl(sb)
	F.x.rsub(m)		

	F.XES=(1<<sb)+1
	if F.XES>FEXCESS {F.reduce()}
}


/* this*=c mod Modulus, where c is a small int */
func (F *FP) imul(c int) {
//	F.norm()
	s:=false
	if (c<0) {
		c=-c
		s=true
	}

	if MODTYPE==PSEUDO_MERSENNE || MODTYPE==GENERALISED_MERSENNE {
		d:=F.x.pxmul(c)
		F.x.copy(mod(d))
		F.XES=2
	} else {
		if F.XES*int32(c)<=FEXCESS {
			F.x.pmul(c)
			F.XES*=int32(c)
		} else {
			n:=NewFPint(c)
			F.Mul(n)
		}
	}
	if s {F.Neg(); F.norm()}
}

/* this*=this mod Modulus */
func (F *FP) sqr() {
	if int64(F.XES)*int64(F.XES)>int64(FEXCESS) {F.reduce()}	
	d:=sqr(F.x)	
	F.x.copy(mod(d))
	F.XES=2
}

/* this+=b */
func (F *FP) Add(b *FP) {
	F.x.add(b.x)
	F.XES+=b.XES
	if (F.XES>FEXCESS) {F.reduce()}
}

/* this-=b */
func (F *FP) Sub(b *FP) {
	n:=NewFPcopy(b)
	n.Neg()
	F.Add(n)
}

func (F *FP) rsub(b *FP) {
	F.Neg()
	F.Add(b)
}

/* this/=2 mod Modulus */
func (F *FP) div2() {
	if (F.x.parity()==0) {
		F.x.fshr(1)
	} else {
		p:=NewBIGints(Modulus);
		F.x.add(p)
		F.x.norm()
		F.x.fshr(1)
	}
}

// See https://eprint.iacr.org/2018/1038
// return this^(p-3)/4 or this^(p-5)/8
func (F *FP) fpow() *FP {
	ac := [11]int{1, 2, 3, 6, 12, 15, 30, 60, 120, 240,255}
	var xp []*FP
// phase 1
	xp=append(xp,NewFPcopy(F)); 
	xp=append(xp,NewFPcopy(F)); xp[1].sqr()
	xp=append(xp,NewFPcopy(xp[1])); xp[2].Mul(F)
	xp=append(xp,NewFPcopy(xp[2])); xp[3].sqr()
	xp=append(xp,NewFPcopy(xp[3])); xp[4].sqr()
	xp=append(xp,NewFPcopy(xp[4])); xp[5].Mul(xp[2])
	xp=append(xp,NewFPcopy(xp[5])); xp[6].sqr()
	xp=append(xp,NewFPcopy(xp[6])); xp[7].sqr()
	xp=append(xp,NewFPcopy(xp[7])); xp[8].sqr()
	xp=append(xp,NewFPcopy(xp[8])); xp[9].sqr()
	xp=append(xp,NewFPcopy(xp[9])); xp[10].Mul(xp[5])
	var n,c int
		
	n=int(MODBITS)
	if MODTYPE==GENERALISED_MERSENNE {    // Goldilocks ONLY
		n/=2;
	}
	if MOD8==5 {
		n-=3
		c=(int(MConst)+5)/8
	} else {
		n-=2
		c=(int(MConst)+3)/4

	}

	bw:=0; w:=1; for w<c {w*=2; bw+=1}
	k:=w-c

	i:=10; key:=NewFPint(0)
		
	if k!=0 {
		for ac[i]>k {i--}
		key.copy(xp[i])
		k-=ac[i]
	}

	for k!=0 {
		i--
		if ac[i]>k {continue}
		key.Mul(xp[i])
		k-=ac[i] 
	}
// phase 2 
	xp[1].copy(xp[2])
	xp[2].copy(xp[5])
	xp[3].copy(xp[10])

	j:=3; m:=8
	nw:=n-bw
	t:=NewFPint(0)
	for 2*m<nw {
		t.copy(xp[j]); j++
		for i=0;i<m;i++ {t.sqr()} 
		xp[j].copy(xp[j-1])
		xp[j].Mul(t)
		m*=2
	}
	lo:=nw-m
	r:=NewFPcopy(xp[j])

	for lo!=0 {
		m/=2; j--
		if lo<m {continue}
		lo-=m
		t.copy(r)
		for i=0;i<m;i++ {t.sqr()}
		r.copy(t)
		r.Mul(xp[j])
	}
// phase 3
	if bw!=0 {
		for i=0;i<bw;i++ {r.sqr()}
		r.Mul(key)
	}

	if MODTYPE==GENERALISED_MERSENNE {     // Goldilocks ONLY
		key.copy(r)
		r.sqr()
		r.Mul(F)
		for i=0;i<n+1;i++ {r.sqr()}
		r.Mul(key)
	}

	return r
}



/* this=1/this mod Modulus */
func (F *FP) Inverse() {

	if MODTYPE==PSEUDO_MERSENNE || MODTYPE==GENERALISED_MERSENNE {
		y:=F.fpow()
		if MOD8==5 {
			t:=NewFPcopy(F)
			t.sqr()
			F.Mul(t)
			y.sqr()
		} 
		y.sqr()
		y.sqr()
		F.Mul(y)
	} else {
		m2:=NewBIGints(Modulus)
		m2.dec(2); m2.norm()
		F.copy(F.Pow(m2))
	}

}

/* return TRUE if this==a */
func (F *FP) Equals(a *FP) bool {
	f:=NewFPcopy(F)
	s:=NewFPcopy(a);

	s.reduce()
	f.reduce()
	if (Comp(s.x,f.x)==0) {return true}
	return false
}


func (F *FP) Pow(e *BIG) *FP {
	var tb []*FP
	var w [1+(NLEN*int(BASEBITS)+3)/4]int8	
	F.norm()
	t:=NewBIGcopy(e)
	t.norm()
	nb:=1+(t.nbits()+3)/4

	for i:=0;i<nb;i++ {
		lsbs:=t.lastbits(4)
		t.dec(lsbs)
		t.norm()
		w[i]=int8(lsbs)
		t.fshr(4);
	}
	tb=append(tb,NewFPint(1))
	tb=append(tb,NewFPcopy(F))
	for i:=2;i<16;i++ {
		tb=append(tb,NewFPcopy(tb[i-1]))
		tb[i].Mul(F);
	}
	r:=NewFPcopy(tb[w[nb-1]])
	for i:=nb-2;i>=0;i-- {
		r.sqr()
		r.sqr()
		r.sqr()
		r.sqr()
		r.Mul(tb[w[i]])
	}
	r.reduce()
	return r
}

/* return sqrt(this) mod Modulus */
func (F *FP) sqrt() *FP {
	F.reduce();
	if MOD8==5 {
		var v *FP
		i:=NewFPcopy(F); i.x.shl(1)
		if MODTYPE==PSEUDO_MERSENNE || MODTYPE==GENERALISED_MERSENNE {
			v=i.fpow()
		} else {
			b:=NewBIGints(Modulus);
			b.dec(5); b.norm(); b.shr(3)
			v=i.Pow(b)
		}

		i.Mul(v); i.Mul(v)
		i.x.dec(1)
		r:=NewFPcopy(F)
		r.Mul(v); r.Mul(i)
		r.reduce()
		return r
	} else {
		var r *FP
		if MODTYPE==PSEUDO_MERSENNE || MODTYPE==GENERALISED_MERSENNE {
			r=F.fpow()
			r.Mul(F)
		} else {
			b:=NewBIGints(Modulus);
			b.inc(1); b.norm(); b.shr(2)
			r=F.Pow(b)
		}
		return r
	}
}

/* return jacobi symbol (this/Modulus) */
func (F *FP) jacobi() int {
	w:=F.redc();
	p:=NewBIGints(Modulus);
	return w.Jacobi(p)
}
