package BLS381

func DelegatePointFromBytes(b []byte) *ECP2 {
	return ECP2_fromBytes(b) // TODO: validation?
}
