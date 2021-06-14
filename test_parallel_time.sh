#!/bin/bash


Ntrials=50
allThreads=(8 4 2)

p=100000000003
coeff0=67183483819,72263994727,14731848372,35549991545,44951725753
#p = 100000000003
#coeff0 = 67183483819
#coeff1 = 72263994727
#coeff2 = 14731848372
#coeff3 = 35549991545
#coeff4 = 44951725753

#p=16752341
#coeff0=69777942,13833081,11033555,1335360,15030548
#coeff1=13833081
#coeff2=11033555
#coeff3=1335360
#coeff4=15030548

for nthread in ${allThreads[@]}; do
	echo ${nthread}
	filename_out=outputs/lmpcmt_p${p}_nthread${nthread}



	for i in $(seq 1 $Ntrials)
	do	
		(./LMPMCT -p ${p}  -f ${coeff0} --nthreads $nthread -T 4) 1>>$filename_out
	done

	for i in {1..10}
	do	
		(./LMPMCT -p ${p}  -f ${coeff0} --nthreads $nthread -T 8) 1>>$filename_out
	done

	for i in {1..10}
	do	
		(./LMPMCT -p ${p}  -f ${coeff0} --nthreads $nthread -T 12) 1>>$filename_out
	done
done
