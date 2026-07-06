.PHONY: package

test_lcm/out.ps: binaries/linux/lcmodel
	cd test_lcm && \
	../binaries/linux/lcmodel < control.file

test_lcm/multi-voxel/multi-voxel.csv: binaries/linux/lcmodel
	cd test_lcm/multi-voxel && \
	../../binaries/linux/lcmodel < control.file

test_lcm/multi-voxel-10/multi-voxel.csv: binaries/linux/lcmodel
	cd test_lcm/multi-voxel-10 && \
	../../binaries/linux/lcmodel < control.file && \
	diff test-reference-multi-voxel.csv multi-voxel.csv 

package: binaries/linux/lcmodel.xz

binaries/linux/lcmodel: source/LCModel.f | binaries/linux/
	gfortran -std=legacy -O3 source/LCModel.f -o binaries/linux/lcmodel

%/:
	mkdir -p $@

%.xz: %
	xz -k $^
