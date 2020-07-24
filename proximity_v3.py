#!usr/bin/python

#Author: ShengShee Thor

#Prgm Name: proximity_v3.py --> Proximity Protocol Version 2

#General Info:
#This prgm scans through all the residues of a user-inputted pdb and determines the nearest neighbors (within <dist_crit> angstroms) of each residue. Residues within a specified region are excluded from all counts. If any atoms of two residues are found to be within <dist_crit> A of each other, centroid side-chains are then used as the reference structures to calculate the "average" side-chain interaction distance between two residues. The # of nearest neighbors and amount of diversity for each residue are also recorded.

#===========================================================================================================================
#First import all the necessary modules
print "\nLoading PyRosetta Modules...\n"

# PyRosetta | Python-2.7 Linux Scientific/Red Hat Linux (64-bit) release
from rosetta import *
init()

from math import *

import os

#===========================================================================================================================
#Create dictionary of all atoms in each of the 20 aa
#Keys = name of aa, values = # of atoms in aa
#Note that the number of atoms in each residue depends on its position along the chain
#N-terminal has 2 more atoms than when in the middle, C-terminal has 1 more atom than when in the middle
res_atom_dict = {'GLY_p:NtermProteinFull':9, 'ALA_p:NtermProteinFull':12, 'VAL_p:NtermProteinFull':18, 'LEU_p:NtermProteinFull':21, 'ILE_p:NtermProteinFull':21, 'PRO_p:NtermProteinFull':17, 'MET_p:NtermProteinFull':19, 'CYS_p:NtermProteinFull':13, 'SER_p:NtermProteinFull':13, 'THR_p:NtermProteinFull':16, 'ASP_p:NtermProteinFull':14, 'GLU_p:NtermProteinFull':17, 'ASN_p:NtermProteinFull':16, 'GLN_p:NtermProteinFull':19, 'LYS_p:NtermProteinFull':24, 'ARG_p:NtermProteinFull':26, 'HIS_p:NtermProteinFull':19, 'PHE_p:NtermProteinFull':22, 'TYR_p:NtermProteinFull':23, 'TRP_p:NtermProteinFull':26, 'GLY_p:CtermProteinFull':8, 'ALA_p:CtermProteinFull':11, 'VAL_p:CtermProteinFull':17, 'LEU_p:CtermProteinFull':20, 'ILE_p:CtermProteinFull':20, 'PRO_p:CtermProteinFull':16, 'MET_p:CtermProteinFull':18, 'CYS_p:CtermProteinFull':12, 'SER_p:CtermProteinFull':12, 'THR_p:CtermProteinFull':15, 'ASP_p:CtermProteinFull':13, 'GLU_p:CtermProteinFull':16, 'ASN_p:CtermProteinFull':15, 'GLN_p:CtermProteinFull':18, 'LYS_p:CtermProteinFull':23, 'ARG_p:CtermProteinFull':25, 'HIS_p:CtermProteinFull':18, 'PHE_p:CtermProteinFull':21, 'TYR_p:CtermProteinFull':22, 'TRP_p:CtermProteinFull':25, 'GLY':7, 'ALA':10, 'VAL':16, 'LEU':19, 'ILE':19, 'PRO':15, 'MET':17, 'CYS':11, 'SER':11, 'THR':14, 'ASP':12, 'GLU':15, 'ASN':14, 'GLN':17, 'LYS':22, 'ARG':24, 'HIS':17, 'PHE':20, 'TYR':21, 'TRP':24}

res_name_dict = {'GLY_p:NtermProteinFull':'GLY', 'ALA_p:NtermProteinFull':'ALA', 'VAL_p:NtermProteinFull':'VAL', 'LEU_p:NtermProteinFull':'LEU', 'ILE_p:NtermProteinFull':'ILE', 'PRO_p:NtermProteinFull':'PRO', 'MET_p:NtermProteinFull':'MET', 'CYS_p:NtermProteinFull':'CYS', 'SER_p:NtermProteinFull':'SER', 'THR_p:NtermProteinFull':'THR', 'ASP_p:NtermProteinFull':'ASP', 'GLU_p:NtermProteinFull':'GLU', 'ASN_p:NtermProteinFull':'ASN', 'GLN_p:NtermProteinFull':'GLN', 'LYS_p:NtermProteinFull':'LYS', 'ARG_p:NtermProteinFull':'ARG', 'HIS_p:NtermProteinFull':'HIS', 'PHE_p:NtermProteinFull':'PHE', 'TYR_p:NtermProteinFull':'TYR', 'TRP_p:NtermProteinFull':'TRP', 'GLY_p:CtermProteinFull':'GLY', 'ALA_p:CtermProteinFull':'ALA', 'VAL_p:CtermProteinFull':'VAL', 'LEU_p:CtermProteinFull':'LEU', 'ILE_p:CtermProteinFull':'ILE', 'PRO_p:CtermProteinFull':'PRO', 'MET_p:CtermProteinFull':'MET', 'CYS_p:CtermProteinFull':'CYS', 'SER_p:CtermProteinFull':'SER', 'THR_p:CtermProteinFull':'THR', 'ASP_p:CtermProteinFull':'ASP', 'GLU_p:CtermProteinFull':'GLU', 'ASN_p:CtermProteinFull':'ASN', 'GLN_p:CtermProteinFull':'GLN', 'LYS_p:CtermProteinFull':'GLN', 'ARG_p:CtermProteinFull':'ARG', 'HIS_p:CtermProteinFull':'HIS', 'PHE_p:CtermProteinFull':'PHE', 'TYR_p:CtermProteinFull':'TYR', 'TRP_p:CtermProteinFull':'TRP'}

#===========================================================================================================================
#Create class
#===========================================================================================================================
#This class creates single objects that contain the following attributes and methods

#Attributes: regions_exclude_status, regions_exclude_dict, bb_exclude_status, atom_dist_calc_start, dist_crit_status, dist_crit, codons_status, codons

#Methods: update_regions_exclude_status, update_regions_exclude_dict, update_bb_exclude_status, update_atom_dist_calc_start, update_dist_crit_status, update_codons_status, update_codons

#An object of this class contains all the info needed to scan a PDB as desired

class pdb_scan_options():
	#Initialize with these attributes
	def __init__(self, name, regions_exclude_status = 'None', regions_exclude_dict = 'None', bb_exclude_status = 'exclude bb', atom_dist_calc_start = 5, dist_crit_status = 'default 2', dist_crit = 2, codons_status = 'default 32', codons = 32):

		self.name = name
		self.regions_exclude_status = regions_exclude_status
		self.regions_exclude_dict = regions_exclude_dict #--> dictionary of regions
		self.bb_exclude_status = bb_exclude_status
		self.atom_dist_calc_start = atom_dist_calc_start
		self.dist_crit_status = dist_crit_status
		self.dist_crit = dist_crit

		self.codons_status = codons_status
		self.codons = codons



	#Update attribute methods
	def update_regions_exclude_status(self, status):
		self.regions_exclude_status = status #status = "exclude" or "None"

	def update_regions_exclude_dict(self, exclude_dict):
		self.regions_exclude_dict = exclude_dict #Enter dictionary of regions to exclude

	def update_bb_exclude_status(self, status):
		self.bb_exclude_status = status #status = "exclude bb" or "count bb"

	def update_atom_dist_calc_start(self, start):
		self.atom_dist_calc_start = start #specify starting atom 1 to include bb, 5 to exclude bb

	def update_dist_crit_status(self, status):
		self.dist_crit_status = status #status = "ready" or "empty"
	
	def update_dist_crit(self, number):
		self.dist_crit = number #number = distance criterion

	def update_codons_status(self, status):
		self.codons_status = status #status = "ready" or "default 32"

	def update_codons(self, number):
		self.codons = number #specify number of codons to consider for library diversity calculation

#===========================================================================================================================
#Define all functions used for User Interface
#===========================================================================================================================

#---------------------------------------------------------------------------------------------------------------------------
#Main menu display

def main_menu(input_pdb_status, input_dir_status, scan_options_status):
    print """
========================================
     Proximity Protocol Version 3
           proximity_v3.py
Identify Nearest Neighbor (NN) residues
========================================
Main Menu"""
    print "1. Input PDB file (status: " + input_pdb_status + ")"
    print "2. Input external data file directory (status: " + input_dir_status + ")" 
    print "3. Scan options (status: " + scan_options_status + ")" #---> need to specify this status
    print "4. Start scanning proteins"
    print "5. Instructions"
    print "6. Exit"

#---------------------------------------------------------------------------------------------------------------------------
#Instructions

def instructions():
	os.system('clear')
	print """
=========================================================
Instructions
=========================================================
1. Input PDB file
	Enter the name of the cleaned PDB file that you wish 
	to analyze. This file should be in the same folder as
	the external data file and the prgm proximity_v2.py
	itself. Python is case-sensitive so make sure you type
	the PDB filename correctly as shown below.
	
	ex/ Name of PDB: 1FNA.pdb

2. Input external data file directory
	An external data file does not need to be created
	beforehand. Just enter an available directory and
	the prgm will automatically create in the specified folder 
	two data files: one for a viewer-friendly data table and 
	one for an exportable data table. Note that if the same
	directory is given as one from a previous run, the program 
	will overwrite the previous text files with new data. 
	(Python is case-sensitive so make sure you type the 
	directory correctly.)

	Linux ex: /home/thors/PyRosetta.ScientificLinux-r49071.64Bit/proximity_data.txt

3. Scan options
	This will take you to the following options where you can customize
	the scanning process:
		1) Use same option set for all PDBs
			-Whatever you specify here will be applied to all the scans
			regardless of PDB file. Generally, use this for similar PDBs
			of the same protein or when you're only scanning 1 PDB.
		2) Customize option set for each PDB"
			-Select the PDB that you want customized scan options for.
			Options for the other PDBs that are not adjusted will take
			on the default settings. Generally, use this feature when 
			you are scanning multiple PDBs of different proteins and need
			to specify regions to exclude.
	
	Note that the scan will use the option set that you last accessed. For example,
	if you last accessed 1), then this option set will be used. If you last accessed
	2), then the customized option set for each PDB will be used.

	Default scan options are listed below:
		-Count all residues (exclude no regions)
		-Exclude bacbone atoms from scan
		-Distance criterion = 2
		-# of codons to consider for library = 32

4. Start protein system scan
	The program starts scanning the protein system for
	nearest neigbors. Note that if a cleaned PDB and an
	external data file directory were not provided, the 
	program will not run.

5. Instructions
	Takes user to these instructions.

6. Exit
	Exit/terminate the prgm. PDBs, directory, and scan
	options info will not be saved for future use.
"""
	raw_input('Hit Enter to return to Main Menu: ')


#---------------------------------------------------------------------------------------------------------------------------
#Input PDBs

def input_pdb_operation(pdbs_options_objects, loaded_pdbs_names):

	#pdbs_options_objects = [] #initialize list of pdb options objects
	enter_pdb = 'notdone'
	while enter_pdb != 'done':
		os.system('clear')
		print """
==============================================================
Cleaned PDB file name
==============================================================
Loaded PDBs: """ + str(loaded_pdbs_names)

		print """
Enter 'clear' to clear PDB list
Enter 'r' to return to Main Menu
"""
		pdb_name = raw_input("Name of PDB (ex: 1FNA.pdb): ")
		if pdb_name == 'r':
			enter_pdb = 'done'

		elif pdb_name == 'clear':
			loaded_pdbs_names = []
			pdbs_options_objects = []
			input_pdb_status = 'empty'
			print "PDB list cleared"

		else:
			pdb_option = pdb_scan_options(pdb_name) #create options object for pdb

			try:
				print "\nChecking PDB...\n"
				protein_system = 0
				protein_system = pose_from_pdb(pdb_option.name) #Initialize the PDB to check	
				
			except:				
				print "\nInvalid PDB file detected: " + pdb_option.name
				print "\n  --> Refer to instructions for further help\n"
				raw_input('Hit Enter to continue: ')
			else:
				print "\nValid PDB inputted\n"
				pdbs_options_objects.append(pdb_option)
				loaded_pdbs_names.append(pdb_option.name)
				
			
	input_pdb_status = str(len(pdbs_options_objects)) + ' loaded'
	input_pdb_update = [input_pdb_status, loaded_pdbs_names, pdbs_options_objects]
	return input_pdb_update #list [status, list of options objects for each pdb]

#---------------------------------------------------------------------------------------------------------------------------
#Input File Directory
#*****Make sure to test if data still records if you put a random data_dir name and run scan*****

def input_file_directory(input_dir_status = 'empty', data_dir = 'None', data_dir_exp = 'None'):
	os.system('clear')
	dummy = 0
	print """
==============================================================
External data file directory
==============================================================
Current data file directory: """ + data_dir

	print """
Make sure that you chose an easily accessed directory otherwise
it will be hard to find your external data text file. Note that
Python is case-sensitive.

Linux ex: /home/msi/thors/PyRosetta.ScientificLinux-r49071.64Bit/proximity_data.txt

Enter 'r' to return to Main Menu
"""
	done = 'no'
	while done != 'yes':
		try:
			entry = raw_input("External data file's directory: ")
		except:
			print "Invalid entry --> Please enter valid directory"
			raw_input('\nHit Enter to continue: ')

		else:
			if entry == 'r':
				if input_dir_status == 'empty':
					done = 'no'
					print "No directory set up"
				else:
					done = 'yes'
			else:
				done = 'yes'
				data_dir = entry
				input_dir_status = 'ready'
				data_dir_exp = data_dir.split('.txt')[0] + '_export.txt' #create file for exportable data

	data_files = [input_dir_status, data_dir, data_dir_exp] #[status, reader-friendly file, exportable file]
	return data_files

#-----------------------------------------------------
#Scan options menu display

def scan_options_menu(scan_options_status, pdbs_options_objects, regions_exclude_status, regions_exclude_dict, bb_exclude_status, atom_dist_calc_start, dist_crit_status, dist_crit, codons_status, codons):
	dummy = 0
	choice = 0
	while choice != 3:
		os.system('clear')
		print """"
========================================
Scan Options Menu
========================================
Status: """ + scan_options_status

		print "\n1. Use same option set for all PDBs"
		print "2. Customize option set for each PDB"
		print "3. Return to main menu"

		try:
			choice = int(raw_input('\nPlease enter the desired option #: '))
		except:
			print "Invalid entry --> Please enter the number of the desired option"
			raw_input('Hit Enter to continue: ')
		else:
			if choice == 3:
				dummy = dummy/1
				if scan_options_status == 'default': #check to see if options were changed
					scan_options = pdbs_options_objects #since options not changed, objects are still at default values
				else:
					dummy = dummy/1 #options have been changed, scan_options is ready to go

			elif choice == 1: #Use same option set for all pdbs
				if len(pdbs_options_objects) == 0:
					print "No PDBs have been loaded"
					raw_input("Hit Enter to continue: ")

				else:
					options_list = scan_options_menu_display_1(regions_exclude_status, regions_exclude_dict, bb_exclude_status, atom_dist_calc_start, dist_crit_status, dist_crit, codons_status, codons)
					#update variables
					regions_exclude_status = options_list[0]
					regions_exclude_dict = options_list[1]
					bb_exclude_status = options_list[2]
					atom_dist_calc_start = options_list[3]
					dist_crit_status = options_list[4]
					dist_crit = options_list[5]
					codons_status = options_list[6]
					codons = options_list[7]

					scan_options = options_list #copy list for return
					scan_options_status = 'same option set used for all PDBs'

			elif choice == 2:
				if len(pdbs_options_objects) == 0:
					print "No PDBs have been loaded"
					raw_input("Hit Enter to continue: ")
					
				else:
					pdbs_options_objects = scan_options_menu_display_2(pdbs_options_objects) #update options for each pdb
					scan_options = pdbs_options_objects #copy objects list for return
					scan_options_status = 'customized for each PDB'

			else:
				print "Option not recognized"
				raw_input('Hit Enter to continue: ')
				
	output = [scan_options_status, scan_options] 
	return output #[status, either a list of defaults for all pdbs or list of customized options for each pdb]

#options_list = [regions_exclude_status, regions_exclude_dict, bb_exclude_status, atom_dist_calc_start, 	dist_crit_status, dist_crit, codons_status, codons]
#---------------------------------------------------------------------------------------------------------------------------
#Scan options menu display for option 1
def scan_options_menu_display_1(regions_exclude_status, regions_exclude_dict, bb_exclude_status, atom_dist_calc_start, dist_crit_status, dist_crit, codons_status, codons):
	dummy = 0
	choice = 0
	while choice != 5:
		os.system('clear')
		print "Same option settings for all PDBs"
		print "==================================="
		print "Current Option Settings\n"
		print "1. Region Exclusion: " + regions_exclude_status
		print "   Regions to exclude --> " + str(regions_exclude_dict)
		print "2. Backbone Exclusion: " + bb_exclude_status
		print "3. Distance Criterion: " + str(dist_crit) + " A"
		print "4. Codons: " + str(codons)
		print "5. Return to Scan Options Menu"

		try:
		    choice = int(raw_input('\nEnter desired option # to change: '))
		except:
		    print "Invalid entry --> Please enter the number of the desired option you want to change"
		    raw_input('\nHit Enter to continue: ')
		else:
		    if choice == 5:
			dummy = dummy/1
	
		    elif choice == 1:
			 regions_exclusion_update = exclude_regions_operation(regions_exclude_status, regions_exclude_dict)
			 regions_exclude_status = regions_exclusion_update[0]
			 regions_exclude_dict = regions_exclusion_update[1]
	
		    elif choice == 2:
			 bb_exclusion_update = exclude_bb_operation()
			 bb_exclude_status = bb_exclusion_update[0]
			 atom_dist_calc_start = bb_exclusion_update[1]
	
		    elif choice == 3:
			 dist_crit_update = set_distance_criterion_operation()
			 dist_crit_status = dist_crit_update[0]
			 dist_crit = dist_crit_update[1]
	
		    elif choice == 4:
			 codons_update = set_codons_operation()
			 codons_status = codons_update[0]
			 codons = codons_update[1]
	
		    else:
			print "Option not recognized"
			raw_input('Hit Enter to continue: ')
	
	options_list = [regions_exclude_status, regions_exclude_dict, bb_exclude_status, atom_dist_calc_start, 	dist_crit_status, dist_crit, codons_status, codons]

	return options_list


#---------------------------------------------------------------------------------------------------------------------------
#Scan options menu display for option 2
def scan_options_menu_display_2(pdb_options_objects):
	dummy = 0

	choice = 0
	customize = 'notdone'
	while customize != 'done':
		os.system('clear')
		print """
==============================================================
Customize options set for each PDB
==============================================================
Enter 'r' to return to previous menu
"""
		print "PDB(s) you have inputted:"
		for i in range(0, len(pdb_options_objects)):
			print pdb_options_objects[i].name


		try:
			choice = raw_input('\nEnter name of inputted PDB to customize (ex/ 1FNA.pdb): ')
		except:
			print "Invalid PDB name --> check spelling and case"
			raw_input('Hit Enter to continue: ')
		else:
			if choice == 'r':
				customize = 'done'
			else:
				find = 'no'
				for i in range(0, len(pdb_options_objects)):
					pdb_name = pdb_options_objects[i].name
					if choice != pdb_name:
						dummy = dummy/1
					elif choice == pdb_name:
						find = 'yes'
						pdb_options = customize_pdb_options(pdb_options_objects[i]) #customize object
						pdb_options_objects[i] = pdb_options #replace old object with updated object

				if find == 'no':
					print "No such PDB loaded"
					raw_input("Hit Enter to continue: ")
				else:
					dummy = dummy/1

	return pdb_options_objects #return updated customized options objects			
			
#---------------------------------------------------------------------------------------------------------------------------
#Customization menu used for Scan options menu display option 2
def customize_pdb_options(a):
	choice = 0
	while choice != 5:
		os.system('clear')

		print "Select which settings to customize for " + a.name
	
		print "\n1. Region Exclusion: " + a.regions_exclude_status
		print "   Regions to exclude --> " + str(a.regions_exclude_dict) 
		print "2. Backbone Exclusion: " + a.bb_exclude_status
		print "3. Distance Criterion: " + str(a.dist_crit) + " A"
		print "4. Codons: " + str(a.codons)
		print "5. Return to Customized Options Menu"

		try:
		    choice = int(raw_input('\nEnter desired option to change: '))

		except:
		    print "Invalid entry --> Please enter the number of the desired option you want to change"
		    raw_input('\nHit Enter to continue: ')

		else:
		    if choice == 5:
			os.system('clear')
	
		    elif choice == 1:
			 regions_exclusion_update = exclude_regions_operation(a.regions_exclude_status, a.regions_exclude_dict)
			 regions_exclude_status = regions_exclusion_update[0]
			 regions_exclude_dict = regions_exclusion_update[1]
			 a.update_regions_exclude_status(regions_exclude_status)
			 a.update_regions_exclude_dict(regions_exclude_dict)
	
		    elif choice == 2:
			 bb_exclusion_update = exclude_bb_operation()
			 bb_exclude_status = bb_exclusion_update[0]
			 atom_dist_calc_start = bb_exclusion_update[1]
			 a.update_bb_exclude_status(bb_exclude_status)
			 a.update_atom_dist_calc_start(atom_dist_calc_start)
	
		    elif choice == 3:
			 dist_crit_update = set_distance_criterion_operation()
			 dist_crit_status = dist_crit_update[0]
			 dist_crit = dist_crit_update[1]
			 a.update_dist_crit_status(dist_crit_status)
			 a.update_dist_crit(dist_crit)
	
		    elif choice == 4:
			 codons_update = set_codons_operation()
			 codons_status = codons_update[0]
			 codons = codons_update[1]
			 a.update_codons_status(codons_status)
			 a.update_codons(codons)
	
		    else:
			print "Option not recognized"
			raw_input('Hit Enter to continue')
	
	return a #updated customized options for single pdb		
		
#---------------------------------------------------------------------------------------------------------------------------
#Update # of codons to consider for library diversity calculations
def set_codons_operation():
	set_codons = 'not done'
	while set_codons != 'done':
		try:
		   codons = input('\nPlease set the desired <codons>: ')
		except:
		   set_codons = 'not done'
		   print "Invalid entry --> Please enter a numerical positive number\n"
		else:
		   if 0 < codons <= 64 :
		      set_codons = 'done'
		      codons_status = 'ready'
		   else:
			set_codons = 'not done'
			print "Invalid # --> Please enter a positive # between 0 and 64"   
			raw_input('Hit Enter to continue: ')

	codons_update = [codons_status, codons]
	return codons_update #list of codon specification [status, # of codons]

#---------------------------------------------------------------------------------------------------------------------------
#Update distance criterion
def set_distance_criterion_operation():
	set_distance = 'not done'
	while set_distance != 'done':
		try:
		  dist_crit = input('\nPlease set the desired <dist_crit>: ')
		except:
		  set_distance = 'not done'
		  print "Invalid entry --> Please enter a numerical positive distance"
		else:
		  if dist_crit > 0:
			set_distance = 'done'
			dist_crit_status = str(dist_crit)
		  else:
			set_distance = 'not done'
			print "Invalid distance --> Please enter a positive distance"
			raw_input('Hit Enter to continue: ')

	dist_crit_update = [dist_crit_status, dist_crit]
	return dist_crit_update #list of distance criterion update [status, distance cutoff]

#---------------------------------------------------------------------------------------------------------------------------
#Update bb exclusion option
def exclude_bb_operation():
	exclude = 0
	while exclude == 0:
		try:
		   exclude = raw_input("\nExclude backbone (y/n)? ")
		except:
		   exclude = 0
		   print "Invalid entry --> Please enter 'y' or 'n'"
		   raw_input('Hit Enter to continue: ')
		else:
		   if exclude == 'y':
			bb_exclude_status = 'exclude bb'
			atom_dist_calc_start = 5
		   elif exclude == 'n':
			bb_exclude_status = 'count bb'
			atom_dist_calc_start = 1
		   else:
			exclude = 0
			print "Invalid entry --> Please enter 'y' or 'n'"
			raw_input('Hit Enter to continue: ')
	
	bb_exclusion_update = [bb_exclude_status, atom_dist_calc_start]
	return bb_exclusion_update # list of bb exclusion update [status, start atom #]

#---------------------------------------------------------------------------------------------------------------------------
#Update regions exclusion option
#***Regions can still overlap***FIX!!!
def exclude_regions_operation(regions_exclude_status, regions_exclude_dict):
	os.system('clear')
	print """
Region Exclusion Settings
-------------------------------------
Current Regions to exclude: """ + str(regions_exclude_dict)
	print "Hit 'r' to return to menu (entered data won't be saved)"

	dummy = 0
	regions_exclude_status = 'disregard region(s)'
	how_many_regions = 'None'#options_list = [regions_exclude_status, regions_exclude_dict, bb_exclude_status, atom_dist_calc_start, 	dist_crit_status, dist_crit, codons_status, codons]
	while how_many_regions == 'None':
		try:
			how_many_regions = int(raw_input("\nHow many unwanted regions? "))
		except:
			how_many_regions = 'None'
			print "Invalid entry --> Please enter a positive integer"
		else:
			if how_many_regions < 0:
				how_many_regions = 'None'
				print "Invalid entry --> Please enter a positive integer"
			elif how_many_regions == 0:
				regions_exclude_status = 'count all residues'
			else:
				dummy = dummy/1

	regions_exclude_dict = {}
	regions_counter = 0
	while regions_counter < how_many_regions:
		regions_counter = regions_counter + 1
		statement_start = "\n Enter the start residue number for Region " + str(regions_counter) + "? "
		statement_end = "\n Enter the end residue number for Region " + str(regions_counter) + "? "
		region_res_start = 0
		while region_res_start == 0:
			try:
				region_res_start = int(raw_input(statement_start))
			except:
				region_res_start = 0
				print "Invalid entry --> Please enter a positive integer"
			else:
				if region_res_start < 0:
					region_res_start = 0
					print "Invalid entry --> Please enter a positive integer"
				else:
					dummy = dummy/1
 
		region_res_end = 0
		while region_res_end == 0:
			try:
				region_res_end = int(raw_input(statement_end))
			except:
				region_res_end = 0
				print "Invalid entry --> Please enter a positive integer"
			else:
				if region_res_end < 0:
					region_res_end = 0
					print "Invalid entry --> Please enter a positive integer"
				elif region_res_end < region_res_start:
					region_res_end = 0
					print "\nInvalid value: "
					print "Ending residue # of region should be greater than starting residue #"
				else:
					dummy = dummy/1
     		  
		regions_exclude_dict.update({regions_counter:[region_res_start, region_res_end]})

#	for i in range(1, len(regions_exclude_dict)):
#		for j in range(1, len(regions_exclude_dict)):
			
    
	region_exclusion_options = [regions_exclude_status, regions_exclude_dict]
	return region_exclusion_options #return list of region exclusion update [status, dictionary]

#---------------------------------------------------------------------------------------------------------------------------



#===========================================================================================================================
#Define all calculation functions that will be used
#---------------------------------------------------------------------------------------------------------------------------

#The first and last residues have weird name labels. This function serves to rename them to the usual 3 letter names
#name.find() method --> if you do find the text _p:, it returns a number not -1. If you don't find _p:, it returns -1
def rename(name):
    dummy = name.find("_p:")
    if dummy == -1:
	return name
    else:
	newname = res_name_dict[name]
	return newname



#---------------------------------------------------------------------------------------------------------------------------
#Check if residue is within unwanted region
def check_if_unwanted(i, j, regions_exclude_status, regions_exclude_dict):
	dummy = 0
	if regions_exclude_status == 'None':
		return

	elif regions_exclude_status == 'exclude':
		for counter in range(1, len(regions_exclude_dict)+1):
			region = regions_exclude_dict[counter]
			region_start = region[0]
			region_end = region[1]
			if region_start <= i <= region_end or region_start <= j <= region_end:
				count_or_not = "exclude"
				return count_or_not
	
			else:
				dummy = dummy/1 #Not in current region

#---------------------------------------------------------------------------------------------------------------------------
#Start New table for next PDB
def data_new_table(scan_options_status, scan_options, pdb_name, options_object, data_dir, data_dir_exp):
	if scan_options_status == 'default' or scan_options_status == 'same option set used for all PDBs': #scan options is list of settings
		regions_exclude_status = scan_options[0]
		bb_exclude_status = scan_options[2]
		dist_crit_status = scan_options[4]
		codons = scan_options[7]		

		#Start the viewer-friendly table
		f = open(data_dir, 'a')
		f.write("Proximity Code Version 3 (proximity_v3.py) Reader-Friendly data file\n====================================================================================\n")
		f.write("OPTION SETTINGS\n\n")
		f.write("PDB File: " + pdb_name + "\n")
		f.write("Data File Directory: " + data_dir + "\n")
		f.write("Region Exclusion Option: " + regions_exclude_status + "\n")
		f.write("Backbone Exclusion Option: " + bb_exclude_status + "\n")
		f.write("Distance Criterion: " + dist_crit_status + "\n")
		f.write("# of codons = " + str(codons) + "\n")
		f.write("====================================================================================\n\n")
		f.write("Selected Residue;\t" + "      NN residue;\t" + "      Centroid Side Chain Distance;\n")
		f.close()

		#Start the exportable data table
		f = open(data_dir_exp, 'a')
		f.write("\n\n\nProximity Protocol Version 3 (proximity_v3.py) Exportable data file \n====================================================================================\n")
		f.write("OPTION SETTINGS\n\n")
		f.write("PDB File: " + pdb_name + "\n")
		f.write("Data File Directory: " + data_dir_exp + "\n")
		f.write("Region Exclusion Option: " + regions_exclude_status + "\n")
		f.write("Backbone Exclusion Option: " + bb_exclude_status + "\n")
		f.write("Distance Criterion: " + dist_crit_status + "\n")
		f.write("# of codons = " + str(codons) + "\n")
		f.write("====================================================================================\n\n")
		f.write("Selected Res;\t" + "Selected Res #;\t" + "NN residue;\t" + "NN residue #;\t" + "Centroid Side Chain Distance;\n")
		f.close()

	else: #scan_options is list of setting objects for each pdb
		regions_exclude_status = options_object.regions_exclude_status
		bb_exclude_status = options_object.bb_exclude_status
		dist_crit_status = options_object.dist_crit_status
		codons = options_object.codons
		
		#Start the viewer-friendly table
		f = open(data_dir, 'a')
		f.write("Proximity Code Version 3 (proximity_v3.py) Reader-Friendly data file\n====================================================================================\n")
		f.write("OPTION SETTINGS\n\n")
		f.write("PDB File: " + pdb_name + "\n")
		f.write("Data File Directory: " + data_dir + "\n")
		f.write("Region Exclusion Option: " + regions_exclude_status + "\n")
		f.write("Backbone Exclusion Option: " + bb_exclude_status + "\n")
		f.write("Distance Criterion: " + dist_crit_status + "\n")
		f.write("# of codons = " + str(codons) + "\n")
		f.write("====================================================================================\n\n")
		f.write("Selected Residue;\t" + "      NN residue;\t" + "      Centroid Side Chain Distance;\n")
		f.close()

		#Start the exportable data table
		f = open(data_dir_exp, 'a')
		f.write("\n\n\nProximity Protocol Version 3 (proximity_v3.py) Exportable data file \n====================================================================================\n")
		f.write("OPTION SETTINGS\n\n")
		f.write("PDB File: " + pdb_name + "\n")
		f.write("Data File Directory: " + data_dir_exp + "\n")
		f.write("Region Exclusion Option: " + regions_exclude_status + "\n")
		f.write("Backbone Exclusion Option: " + bb_exclude_status + "\n")
		f.write("Distance Criterion: " + dist_crit_status + "\n")
		f.write("# of codons = " + str(codons) + "\n")
		f.write("====================================================================================\n\n")
		f.write("Selected Res;\t" + "Selected Res #;\t" + "NN residue;\t" + "NN residue #;\t" + "Centroid Side Chain Distance;\n")
		f.close()
		

#options_list = [regions_exclude_status, regions_exclude_dict, bb_exclude_status, atom_dist_calc_start, 	dist_crit_status, dist_crit, codons_status, codons]
#---------------------------------------------------------------------------------------------------------------------------

#Data recording function for reader-friendly (rf) file
def data_rec_rf(protein_system, i, j, dist, data_dir):
	dummy = 0
	res_select = rename(protein_system.residue(i).name())
	res_select_IDnum = str(i)
	if j == 0:
		res_neighbor = str(j)
		res_neighbor_IDnum = str(j)
	else:
		res_neighbor = rename(protein_system.residue(j).name())
		res_neighbor_IDnum = str(j)

	f = open(data_dir, 'a')
	f.write(res_select + "  " + res_select_IDnum + '   \t\t\t' + res_neighbor + "  " + res_neighbor_IDnum + '   \t\t   ' + dist + '\n' )
	f.close()

#---------------------------------------------------------------------------------------------------------------------------

#Data recording function to make exportable file
def datarec_export(protein_system, i, j, dist, data_dir_exp):
	res_select = rename(protein_system.residue(i).name())
	res_select_IDnum = str(i)
	if j == 0:
		res_neighbor = str(j)
		res_neighbor_IDnum = str(j)
	else:
		res_neighbor = rename(protein_system.residue(j).name())	
		res_neighbor_IDnum = str(j)

	f = open(data_dir_exp, 'a')
	f.write(res_select + ";\t\t" + res_select_IDnum + ";\t\t" + res_neighbor + ";\t\t" + res_neighbor_IDnum + ";\t\t" + dist + "\n")
	f.close()

#---------------------------------------------------------------------------------------------------------------------------

#Function to start new row for next residue in reader-friendly file
def data_start_new_res(num_of_nn, i_diversity, data_dir):
    num_of_nn = str(num_of_nn)
    i_diversity_sn = str('%e' % i_diversity) #turn into scientific notation
    i_diversity = str(i_diversity)
    f = open(data_dir, 'a')
    f.write("\n")
    f.write("# of NN = " + num_of_nn + "	Diversity = " + i_diversity + " (" + i_diversity_sn + ")\n")
    f.write("------------------------------------------------------------------------------------\n")

#---------------------------------------------------------------------------------------------------------------------------

#Centroid distance calculator          
def cen_dist_calc(i, j, pose):
    dummy = 0
    check = pose.is_centroid()
    if check is False:
    	to_centroid = SwitchResidueTypeSetMover('centroid')
    	to_centroid.apply(pose)
    elif check is True:
	dummy = dummy/1

    cen_i_coord = pose.residue(i).xyz('CEN')
    cen_j_coord = pose.residue(j).xyz('CEN')
    dist = sqrt((cen_i_coord.x - cen_j_coord.x)**2 + (cen_i_coord.y - cen_j_coord.y)**2 + (cen_i_coord.z - cen_j_coord.z)**2)
    dist_rec = str(dist)
    return dist_rec    


#---------------------------------------------------------------------------------------------------------------------------	

#Atom-atom distance calculator

#Input the number of atoms for both residues
#Calculate distances between all possible combo of atoms between the two residue
#If there is a distance < 2 A, stop calculations and return the distance
#If distance > 2 A, check to see if all atoms have been compared
#If yes, then terminate process and output that there is no atom pairs with distance < 2 A
def atom_dist_calc(pose, res_i_num, res_i_atoms, res_j_num, res_j_atoms, dist_crit, atom_dist_calc_start):
	dummy = 0
	check = pose.is_fullatom() #check if pose is in full atom mode
	if check is False:
		to_fullatom = SwitchResidueTypeSetMover('fa_standard')
		to_fullatom.apply(pose)
	elif check is True:
		dummy = dummy/1
	
	for i_atom in range(atom_dist_calc_start, res_i_atoms+1): #select atom of residue i
		i_atom_coord = pose.residue(res_i_num).xyz(i_atom) #retrieve atomic coordinates of i_atom
		for j_atom in range(atom_dist_calc_start, res_j_atoms+1): #select atom of residue j
			j_atom_coord = pose.residue(res_j_num).xyz(j_atom) #retrieve atomic coordinates of j_atom
			#compute the distance between i_atom and j_atom using distance formula
			distance = sqrt((i_atom_coord.x - j_atom_coord.x)**2 + (i_atom_coord.y - j_atom_coord.y)**2 + 						(i_atom_coord.z - j_atom_coord.z)**2)
			if distance > dist_crit: 
				#if distance is greater than the cutoff distance, check that you're not at the last atoms of both residues
				if i_atom == res_i_atoms and j_atom == res_j_atoms:
					return
					
			elif distance < dist_crit:
				return distance

#===========================================================================================================================
def single_options_set_scan(input_pdb_status, loaded_pdbs_names, scan_options, scan_options_status, data_dir, data_dir_exp, res_atom_dict):
	dummy = 0
	if input_pdb_status == 'empty':
		print "No PDBs loaded --> Please load PDBs and set options first"
		raw_input('Hit Enter to continue: ')
	else:
		#Load relevant scan settings
		regions_exclude_status = scan_options[0]
		regions_exclude_dict = scan_options[1] 
		atom_dist_calc_start = scan_options[3]
		dist_crit = scan_options[5]
		codons = scan_options[7]

		for pdb in range(0, len(loaded_pdbs_names)):
			pdb_name = loaded_pdbs_names[pdb] #object containing scan options for particular pdb
			print "Starting scan on " + pdb_name
			protein_system = pose_from_pdb(pdb_name) #load pose
			total = protein_system.total_residue()
			total_diversity = 0
			data_new_table(scan_options_status, scan_options, pdb_name, dummy, data_dir, data_dir_exp)

			#start scan for single pdb
			for i in range(1, total+1):
				res_i_name = rename(protein_system.residue(i).name())
				res_i_atoms = res_atom_dict[res_i_name]
				num_of_nn = 0
				i_diversity = 0
				for j in range(1, total+1):
					res_j_name = rename(protein_system.residue(j).name())
					res_j_atoms = res_atom_dict[res_j_name]
					if i == j:
						print "res_i = res_j: no comparison made"
					elif i != j:
						count_or_not = check_if_unwanted(i, j, regions_exclude_status, regions_exclude_dict)
						if count_or_not == "exclude":
							print res_i_name, i, "+", res_j_name, j, "pair not counted --> one or both found in a region"
						elif count_or_not is None:
							atom_dist = atom_dist_calc(protein_system, i, res_i_atoms, j, res_j_atoms, dist_crit, atom_dist_calc_start)
							if atom_dist is None:
								dummy = dummy/1
							elif atom_dist < dist_crit:
								print "residue", res_i_name, i, "and residue", res_j_name, j, "are close neighbors"
								num_of_nn = num_of_nn + 1
								cen_dist = cen_dist_calc(i, j, protein_system)
								data_rec_rf(protein_system, i, j, cen_dist, data_dir)
								datarec_export(protein_system, i, j, cen_dist, data_dir_exp)

				i_diversity = codons**(num_of_nn + 1)
				total_diversity = total_diversity + i_diversity
				data_start_new_res(num_of_nn, i_diversity, data_dir)


			total_diversity_sn = str('%e' % total_diversity)
			total_diversity = str(total_diversity)
			f = open(data_dir, 'a')
			f.write("\n" + "Total diversity = " + total_diversity + " (" + total_diversity_sn + ")")
			f.write("\n\n=======================================================================================\n")
			f.close()			
			
			os.system('clear')
			print "\n Scan completed for " + pdb_name + "\n"
			print "Total Diversity = " + total_diversity + " (" + total_diversity_sn + ")\n"


def customized_option_sets_scan(input_pdb_status, loaded_pdbs_names, scan_options, scan_options_status, data_dir, data_dir_exp, res_atom_dict):
	dummy = 0
	if input_pdb_status == 'empty':
		print "No PDBs loaded --> Please load PDBs and set options first"
		raw_input('Hit Enter to continue: ')
	else:

		for pdb in range(0, len(scan_options)):
			options_object = scan_options[pdb] #object containing scan options for particular pdb
			pdb_name = options_object.name
			regions_exclude_status = options_object.regions_exclude_status
			regions_exclude_dict = options_object.regions_exclude_dict 
			atom_dist_calc_start = options_object.atom_dist_calc_start
			dist_crit = options_object.dist_crit
			codons = options_object.codons

			print "Starting scan on " + pdb_name
			protein_system = pose_from_pdb(pdb_name) #load pose
			total = protein_system.total_residue()
			total_diversity = 0
			data_new_table(scan_options_status, scan_options, pdb_name, options_object, data_dir, data_dir_exp)

			#start scan for single pdb
			for i in range(1, total+1):
				res_i_name = rename(protein_system.residue(i).name())
				res_i_atoms = res_atom_dict[res_i_name]
				num_of_nn = 0
				for j in range(1, total+1):
					res_j_name = rename(protein_system.residue(j).name())
					res_j_atoms = res_atom_dict[res_j_name]
					if i == j:
						print "res_i = res_j: no comparison made"
					elif i != j:
						count_or_not = check_if_unwanted(i, j, regions_exclude_status, regions_exclude_dict)
						if count_or_not == "exclude":
							print res_i_name, i, "+", res_j_name, j, "pair not counted --> one or both found in a region"
						elif count_or_not is None:
							atom_dist = atom_dist_calc(protein_system, i, res_i_atoms, j, res_j_atoms, dist_crit, atom_dist_calc_start)
							if atom_dist is None:
								dummy = dummy/1
							elif atom_dist < dist_crit:
								print "residue", res_i_name, i, "and residue", res_j_name, j, "are close neighbors"
								num_of_nn = num_of_nn + 1
								cen_dist = cen_dist_calc(i, j, protein_system)
								data_rec_rf(protein_system, i, j, cen_dist, data_dir)
								datarec_export(protein_system, i, j, cen_dist, data_dir_exp)

				i_diversity = codons**(num_of_nn + 1)
				total_diversity = total_diversity + i_diversity
				data_start_new_res(num_of_nn, i_diversity, data_dir)

			total_diversity_sn = str('%e' % total_diversity)
			total_diversity = str(total_diversity)
			f = open(data_dir, 'a')
			f.write("\n" + "Total diversity = " + total_diversity + " (" + total_diversity_sn + ")")
			f.write("\n\n=======================================================================================\n")
			f.close()			
			
			f = open(data_dir_exp, 'a')
			f.write("\n\n=======================================================================================\n")
			f.close()

			os.system('clear')
			print "\n Scan completed for " + pdb_name + "\n"
			print "Total Diversity = " + total_diversity + " (" + total_diversity_sn + ")\n"

#===========================================================================================================================
#Main Body Script
#===========================================================================================================================

#Default values
dummy = 0
input_pdb_status = 'empty'
loaded_pdbs_names = []
pdbs_options_objects = []
input_dir_status = 'empty'
data_dir = 'None'
data_dir_exp = 'None'
regions_exclude_status = 'count all residues'
regions_exclude_dict = {}
bb_exclude_status = 'exclude'
atom_dist_calc_start = 5
dist_crit_status = 'default 2'
dist_crit = 2
codons_status = 'default 32'
codons = 32
scan_options_status = 'default'
scan_options = [regions_exclude_status, regions_exclude_dict, bb_exclude_status, atom_dist_calc_start, dist_crit_status, dist_crit, codons_status, codons]
choice = 0


#options_list = [regions_exclude_status, regions_exclude_dict, bb_exclude_status, atom_dist_calc_start, 	dist_crit_status, dist_crit, codons_status, codons]

while choice != 4:
	os.system('clear')
	main_menu(input_pdb_status, input_dir_status, scan_options_status)
	try:
		choice = int(raw_input('\nPlease enter the desired option #: '))
	except:
		print "Invalid entry --> Please enter the number of the desire option"
		raw_input('\nHit Enter to continue: ')
	else:
	 #================================================================================

		if choice == 6:
			quit()

		elif choice == 1: #input pdb option
			input_pdb_update = input_pdb_operation(pdbs_options_objects, loaded_pdbs_names)
			input_pdb_status = input_pdb_update[0]
			loaded_pdbs_names = input_pdb_update[1]
			pdbs_options_objects = input_pdb_update[2] #pull out the pdbs_options_list

		elif choice == 2: #input directory option
			data_files = input_file_directory(input_dir_status, data_dir, data_dir_exp)
			input_dir_status = data_files[0]
			data_dir = data_files[1]
			data_dir_exp = data_files[2]
			try:
				f = open(data_dir, 'w')
				f.write('\n')
				f.close()
				f = open(data_dir_exp, 'w')
				f.write('\n')
				f.close()
			except:
				print "Invalid directory probably given"
				print "Check directory entry again"
				raw_input("Hit Enter to continue: ")
			else:
				dummy = dummy/1
 
	
		elif choice == 3:
			output = scan_options_menu(scan_options_status, pdbs_options_objects, regions_exclude_status, regions_exclude_dict, bb_exclude_status, atom_dist_calc_start, dist_crit_status, dist_crit, codons_status, codons)
			scan_options_status = output[0] #goes back to main menu input
			scan_options = output[1] #goes down to scan script below
	
		elif choice == 5:
			instructions()

		elif choice == 4:
			if input_pdb_status == 'empty':
				choice = 0
				print "No PDBs loaded --> Please load PDBs and set options first"
				raw_input('Hit Enter to continue: ')

			elif input_dir_status == 'empty':
				choice = 0
				print "No data file directory loaded --> Please load data file directory first"
				raw_input('Hit Enter to continue: ')				

			else:

				if scan_options_status == 'default' or scan_options_status == 'same option set used for all PDBs': #scan options is list of settings
					single_options_set_scan(input_pdb_status, loaded_pdbs_names, scan_options, scan_options_status, data_dir, data_dir_exp, res_atom_dict)
			
				else:
					customized_option_sets_scan(input_pdb_status, loaded_pdbs_names, scan_options, scan_options_status, data_dir, data_dir_exp, res_atom_dict)

		else:
			print "Option not recognized"
			raw_input('Hit Enter to continue: ')


print """
Scanning complete
Please see the following directories for the data:
"""
print "reader-friendly: " + data_dir
print "exportable: " + data_dir_exp
print "\n"





