

repeattime=10


for dim in 0 3
do
	mkdir tmpagresult.valueonly
	echo "========DIM $dim========="
	for method in naive naive_trained pl worker
	do		
		for setting in m_100 m_200 m_300 m_400 m_500 m_600 m_700 m_800 m_900 m_1000
		do
			agresult=tmpagresult.valueonly/$setting.$method.aggregated.result
			for (( iter=1;  iter<=$repeattime; iter++ ))
			do
				cat $setting/iter_$iter/$method.result >> $agresult
			        echo >> $agresult
	                done
			python avgeval.py $agresult $dim
		done
		echo
	done
	rm -r tmpagresult.valueonly
done
