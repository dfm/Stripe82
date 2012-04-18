default:

db:
	bin/populate_db.py --stars --fields --clobber --get

pieces:
	bin/

clean:
	rm -rf calibration/_likelihood.so calibration/gp/_gp.so build

