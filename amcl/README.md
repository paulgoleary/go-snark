### AMCL Notes

The `amcl` Go package is generated from the Milagro open source crypto project: https://github.com/milagro-crypto/amcl

Under `version3/go` there is a Python utility `config64.py` which allows you to selectively generate the curve systems you want.
Right now we generate `SECP256K1` and `BLS381`.

Some cleanup of the generated code is necessary to get it to compile and run. For example:

* The package name of some imports have to be changed to point to this project.
* The file `ROM_SEC256K1_64` seems to be misnamed; it should be `ROM_SECP256K1_64`. The generation tool looks for a file with that name and
there is no such curve name anyway.

There may be some other nits I'm forgetting. The amount of code changes seem mininal and obvious so I don't think it's worth trying automate the generation and cleanup.

My main focus is to assure that the generated implementations match their C counterparts that we rely on elsewhere.
So far compatibility seems good.