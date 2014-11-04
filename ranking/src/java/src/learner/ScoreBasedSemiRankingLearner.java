package learner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map.Entry;

import model.RankingDataSet;
import model.SemiRankingDataSet;



public abstract class ScoreBasedSemiRankingLearner {
	SemiRankingDataSet rkdata;
	
	double[] scores;
	
	abstract void trainRankings();
	
	public ArrayList<Integer> getInferredPositiveIdList() {
		ArrayList<Integer> ret = new ArrayList<Integer>();
		for (int i = 0; i < scores.length; ++i) 
			if (scores[i] > 0.5) ret.add(i);
		return ret;
	}
	
	public void evaluate(String gtFileName) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(new File(gtFileName)));
		HashSet<Integer> gtSet = new HashSet<Integer>();
		String s;
		while ((s = br.readLine()) != null) {
			s = s.trim();
			if (s == "") continue;
			gtSet.add(rkdata.name2Id.get(s));
		}
		br.close();
		System.out.println(gtSet.size());
		
		int tp = 0;
		ArrayList<Integer> inferredList = getInferredPositiveIdList(); 
		for (int id : inferredList) if (gtSet.contains(id)) ++tp;
		int fp = inferredList.size() - tp;
		int fn = gtSet.size() - tp;
		double precision = (double) tp / (tp + fp);
		double recall    = (double) tp / (tp + fn);
		double f1        = 2 * precision * recall / (precision + recall);
		System.out.println("TP = " + tp + ", FP = " + fp + ", FN = " + fn);
		System.out.println("Precision = " + precision);
		System.out.println("Recall    = " + recall);
		System.out.println("F1-score  = " + f1);
	}
}
