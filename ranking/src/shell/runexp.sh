repeattime=10
trainsrk=job_509470.json.train.srk
traingt=fullres.out.train.filtered
testsrk=job_509470.json."test".srk
testgt=fullres.out."test".filtered


for mtrain in 600 700 800 900 1000 #100 200 300 400 500
do
	echo "mtrain=$mtrain"
	expdir=m_"$mtrain"
	mkdir $expdir
	cp Tester.jar sample_train.py $trainsrk $testsrk $traingt $testgt $expdir/
	cd $expdir
	for (( iter=1;  iter<=$repeattime; iter++ ))
	do
		iterexp=iter_$iter
		echo "=====>>>>>> ITER $iter <<<<<<====="
		mkdir $iterexp
		cp Tester.jar sample_train.py $trainsrk $testsrk $traingt $testgt $iterexp/
		cd $iterexp
			python sample_train.py $trainsrk $traingt $mtrain
			for method in naive naive_trained pl worker
			do
				echo "====METHOD $method======"
				java -jar Tester.jar -data $testsrk -gt $testgt -traindata $trainsrk.sample -traingt $traingt.sample -result $method.result -method $method -eval prf
			done
		cd ..
	done
	cd ..
done
