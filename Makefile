all: likelihood

clean:
	rm -f stripe82/calibration/_likelihood.so

likelihood: stripe82/calibration/_likelihood.c
	python setup.py build_ext --inplace
