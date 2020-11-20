
import os.path


p = 100000000003
nthreads = [1,2,4,8,12]


Npoint_var = 0
for nthread in nthreads:
	fname="outputs/lmpcmt_p"+str(p)+"_nthread"+str(nthread)
	if os.path.isfile(fname):
		ff = open(fname,'r')
		for line in ff:
			line = line.rstrip()
			npoints, time = line.split(' ')

			#verify that the number of points is the same
			if Npoint_var == 0:
				Npoint_var = npoints
			else:
				assert Npoint_var == npoints

			# TODO: compute the average runtime

			
		ff.close()
	else:
		print(fname,' does not exist')

