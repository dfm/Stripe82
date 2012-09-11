import os
import sys

S82_ROOT = "/home/dfm265/projects/s82"

sys.path.insert(0, os.path.join(S82_ROOT))
sys.path.insert(0, os.path.join(S82_ROOT, "webapp"))

activate_this = os.path.join(S82_ROOT, "venv/bin/activate_this.py")
execfile(activate_this, dict(__file__=activate_this))

from app import app as application
