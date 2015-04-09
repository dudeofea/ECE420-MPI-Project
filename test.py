#!/usr/bin/python

#	Python Script for Testing
#	MPI version against serial

from subprocess import check_output

#Run serial
check_output(['./serialcalc'], shell=True)
check_output(['mv data_output serial_output'], shell=True)

#Run parallel
check_output(['make run'], shell=True)
check_output(['mv data_output mpi_output'], shell=True)

#Compare files
diffs = 0
with open('serial_output', 'r') as sfile:
	with open('mpi_output', 'r') as mfile:
		while 1:
			sline = sfile.readline().replace('\n', '')
			mline = mfile.readline().replace('\n', '')
			if not sline or not mline:
				break
			if sline != mline:
				diffs += 1
				#print "diff:", sline, mline
if diffs == 0:
	print "All good"
else:
	print diffs, "differences found"

#Delete files
check_output(['rm mpi_output serial_output'], shell=True)