import numpy as np
import scipy
import re
import subprocess
import os

from . import Cfv
from . import const
from . import parser
from . import kin

# guarantee that script is exec
def make_executable(path):
	mode = os.stat(path).st_mode
	mode |= (mode & 0o444) >> 2
	os.chmod(path, mode)

# run shell commands from notebook
def subprocess_cmd(command, verbose=2):
	process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout,stderr = process.communicate()
	if verbose==2:	
		print(command)
		print(stdout.decode("utf-8"))
		print(stderr.decode("utf-8"))
	elif verbose==1:
		if len(stderr.decode("utf-8"))>2:
			print(command)
			print('n',stderr.decode("utf-8"),'m')


class exp_setup:

	def __init__(self, name, Z=None, A=None, channels=None, data_path=""):

		self.name = name
		self.path_to_eventfiles = f"{data_path}data/{self.name}/"
		self.Z = Z
		self.A = A
		self.channels = channels

		self.Erange = np.array([0.1,20.0])
		self.fluxfile=None

		self.new_physics=False

	def set_new_physics(self, gprimeV=0, gprimeA=0, mzprime=0):

		self.new_physics=True

		self.gprimeV=gprimeV
		self.gprimeA=gprimeA
		self.mzprime=mzprime

	def get_number_events(self):
		return np.sum(self.all_events["weights"])*self.norm

	def get_kin(self):
		return self.kin

	def run_MC(self, nevents=1e4, include_pel=True, include_nel=True, verbose=2):
		
		self.nevents=nevents

		if self.new_physics:
			self.global_arguments = fr""" --path="{self.path_to_eventfiles}" """\
								fr"""--fluxfile="{self.fluxfile}" """\
								f"--emin={self.Erange[0]} --emax={self.Erange[1]} "\
								f"-N {int(self.nevents)} "\
								f"--gprimeV={self.gprimeV} "\
								f"--gprimeA={self.gprimeA} "\
								f"--mzprime={self.mzprime} "
		else:
			self.global_arguments = fr"""--path="{self.path_to_eventfiles}" """\
								fr"""--fluxfile="{self.fluxfile}" """\
								f"--emin={self.Erange[0]} --emax={self.Erange[1]} "\
								f"-N {int(self.nevents)} "

		if not (self.Z is None):
			for c in self.channels:
				for i in range(np.size(self.Z)):
					

					# nucleus or hydrogen
					command = [f"cd ../; "\
								f"./gen_SM "\
								f"-c {c} "\
								f"-z {self.Z[i]} "\
								f"-a {self.A[i]} "\
								+self.global_arguments]
					subprocess_cmd(command, verbose=verbose)

					
					# protons in nucleus				
					if self.Z[i]>1 and include_pel:
						command_proton = [f"cd ../; "\
									f"./gen_SM "\
									f"--pb "\
									f"-c {c} "\
									f"-z {1} "\
									f"-a {1} "\
									+self.global_arguments]
						subprocess_cmd(command_proton, verbose=verbose)
					
					
					# neutrons in nucleus				
					if self.Z[i]>1 and include_nel:
						command_neutron = [f"cd ../; "\
									f"./gen_SM "\
									f"--pb "\
									f"-c {c} "\
									f"-z {0} "\
									f"-a {1} "\
									+self.global_arguments]
						subprocess_cmd(command_neutron, verbose=verbose)


			# now load all files and combine them into a single dictionary with approprate weights
			self.dic 		= parser.load_files_to_dic(self)
			self.all_events = parser.combine_MC_outputs(self.dic, self)
			self.kin        = kin.kin(self.all_events, const.m_mu, const.m_mu, Emin=self.Erange[0], Emax=self.Erange[1])

			print(f"Succesfully obtained events for {self.name}.")
		else:
			print(f"Error! First set scattering targets for {self.name}")




