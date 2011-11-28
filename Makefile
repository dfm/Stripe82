all: likelihood

clean:
	rm -rf calibration/_likelihood.so calibration/gp/_gp.so build

cython: calibration/gp/_gp.pyx
	cython calibration/gp/_gp.pyx

likelihood: calibration/_likelihood.c calibration/gp/_gp.c
	python setup.py build_ext --inplace
