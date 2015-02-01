repeattime=10
ntrain=1000
ntest=5000

lambda=0.5
rho=2

method=worker

for nmratio in 2 5 10 20
do
	#echo model_"$lambda"_"$rho".txt
	mtrain=$(( $ntrain * $nmratio ))
	mtest=$(( $ntest * $nmratio ))
	echo "ntrain=$ntrain, mtrain=$mtrain;  ntest=$ntest, mtest=$mtest"
	expdir=train_n_"$ntrain"_m_"$mtrain"_test_n_"$ntest"_m_"$mtest"
#	mkdir $expdir
#	cp Tester.jar datagen_beta.py $expdir/
	cd $expdir
	for (( iter=1;  iter<=$repeattime; iter++ ))
	do
		iterexp=iter_$iter
		echo "=====>>>>>> ITER $iter <<<<<<====="
#		mkdir $iterexp
#		cp Tester.jar datagen_beta.py $iterexp
		cd $iterexp
#			python datagen_beta.py $lambda $rho model_"$lambda"_"$rho".txt $ntrain $mtrain rkdata_"$lambda"_"$rho"_train.txt gtlabel_"$lambda"_"$rho"_train.txt gtscore_"$lambda"_"$rho"_train.txt
#			python datagen_beta.py $lambda $rho model_"$lambda"_"$rho".txt $ntest  $mtest  rkdata_"$lambda"_"$rho"_test.txt  gtlabel_"$lambda"_"$rho"_test.txt  gtscore_"$lambda"_"$rho"_test.txt   
			for method in naive
                        do
				echo "====METHOD $method======"
				java -jar Tester.jar -data rkdata_"$lambda"_"$rho"_test.txt -gt gtlabel_"$lambda"_"$rho"_test.txt -traindata rkdata_"$lambda"_"$rho"_train.txt -traingt gtscore_"$lambda"_"$rho"_train.txt -result $method.result -method $method -eval prf
			done
		cd ..
	done
	cd ..
done
