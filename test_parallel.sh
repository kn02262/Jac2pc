#!/bin/bash


Ntrials=10
allThreads=(1 2 4 8 12)

#p = 100000000003
#coeff0 = 67183483819
#coeff1 = 72263994727
#coeff2 = 14731848372
#coeff3 = 35549991545
#coeff4 = 44951725753

p=16752341
coeff0=69777942,13833081,11033555,1335360,15030548
#coeff1=13833081
#coeff2=11033555
#coeff3=1335360
#coeff4=15030548

for nthread in ${allThreads[@]}; do
	echo ${nthread}
	filename_out=lmpcmt_p${p}_nthread${nthread}

	if [ -f $filename_out ] ; then
    	rm $filename_out
	fi

	for i in {1..10}
	do	
		(time ./LMPMCT -p ${p}  -f ${coeff0} --nthreads $nthread) 2>>$filename_out
	done
done