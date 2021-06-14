
import os.path


p = 16752341
nthreads = [1,2,4,6,8,12]


Npoint_var = 0
for nthread in nthreads:
	fname="outputs/lmpcmt_p"+str(p)+"_nthread"+str(nthread)
	average_time = 0
	if os.path.isfile(fname):
		ff = open(fname,'r')
		ntrials = 0
		for line in ff:
			line = line.rstrip()
			npoints, time = line.split(' ')

			#verify that the number of points is the same
			if Npoint_var == 0:
				Npoint_var = npoints
			else:
				assert Npoint_var == npoints

			# TODO: compute the average runtime
			average_time += float(time)
			ntrials+=1

		average_time = (average_time)/ntrials
		print('for p = ', p, ' on nthreads = ', nthread, ' average_time = ', average_time) 
		ff.close()
	else:
		print(fname,' does not exist')

