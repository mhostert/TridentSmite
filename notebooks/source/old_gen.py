import numpy as np
import subprocess
import os

from . import const
from . import parser

# guarantee that script is exec
def make_executable(path):
	mode = os.stat(path).st_mode
	mode |= (mode & 0o444) >> 2
	os.chmod(path, mode)

# run shell commands from notebook
def subprocess_cmd(command, verbose=2):
	process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout,stderr = process.communicate()
	if verbose=2:	
		print(stdout.decode("utf-8"))
		print(stderr.decode("utf-8"))
	elif verbose=1:	
		print(stderr.decode("utf-8"))
#
def run_MC(exp,gprime,mzprime, verbose=True):
	make_executable(f"../exps/{exp}.sh")
	command = [f"cd ../exps/; ./{exp}.sh {gprime} {mzprime}"]
	subprocess_cmd(command, verbose=verbose)

def generate_tridents(EXP,gprime,mzprime,materials, verbose=True, old_dic=False):
	run_MC(EXP,gprime,mzprime, verbose=verbose)
	dic=parser.load_files_to_dic(f'../data/{EXP}/',gprime,mzprime, gaprime=0.0)
	new_dic=parser.combine_MC_outputs(dic,materials)
	if old_dic:
		return [dic,new_dic]
	else:
		return new_dic


