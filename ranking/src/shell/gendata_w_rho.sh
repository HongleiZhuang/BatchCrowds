lambda=0.5
for rho in 1 1.5 2 2.5 3
do
	echo $rho
	echo model_"$lambda"_"$rho".txt
	python datagen_beta.py $lambda $rho model_"$lambda"_"$rho".txt 1000 10000 rkdata_"$lambda"_"$rho".txt gtlabel_"$lambda"_"$rho".txt gtscore_"$lambda"_"$rho".txt   
done
