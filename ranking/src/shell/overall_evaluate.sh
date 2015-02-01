

repeattime=10
mkdir tmpagresult


for setting in train_n_*_m_*_test_n_*_m_*
do
	echo "============$setting=============="
	for method in naive naive_trained pl worker
	do
		echo "+++++$method+++++++"
		agresult=tmpagresult/$setting.$method.aggregated.result
		for (( iter=1;  iter<=$repeattime; iter++ ))
		do
			cat $setting/iter_$iter/$method.result >> $agresult
		        echo >> $agresult
                done
		python avgeval_readable.py $agresult
	done
done
