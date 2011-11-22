all: likelihood

clean:
	rm -f calibration/_likelihood.so

likelihood: calibration/_likelihood.c calibration/gp/_gp.c
	python setup.py build_ext --inplace
