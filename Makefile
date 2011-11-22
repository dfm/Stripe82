all: likelihood

clean:
	rm -f calibration/_likelihood.so

source: calibration/gp/_gp.pyx
	cython calibration/gp/_gp.pyx

likelihood: calibration/_likelihood.c calibration/gp/_gp.c
	python setup.py build_ext --inplace
