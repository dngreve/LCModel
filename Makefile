.PHONY: package

test_lcm/out.ps: binaries/linux/lcmodel
	cd test_lcm && \
	../binaries/linux/lcmodel < control.file

package: binaries/linux/lcmodel.xz

binaries/linux/lcmodel: source/LCModel.f | binaries/linux/
	gfortran -std=legacy -O3 source/LCModel.f -o binaries/linux/lcmodel

%/:
	mkdir -p $@

%.xz: %
	xz -k $^
