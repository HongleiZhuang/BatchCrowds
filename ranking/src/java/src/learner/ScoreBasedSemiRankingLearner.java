package learner;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import util.MapSorter;
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
		HashSet<Integer> gtSet = readGt(gtFileName);
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
	
	public void evaluateByROC(String gtFileName) throws Exception {
		HashSet<Integer> gtSet = readGt(gtFileName);
		HashMap<Integer, Double> name2Score = new HashMap<Integer, Double>();
		for (int i = 0; i < scores.length; ++i) {
			name2Score.put(i, scores[i]);
		}
		ArrayList<Entry<Integer, Double>> sortedList = MapSorter.sortByValue(name2Score, false);
		double AUC = 0.0;
		int tp = 0;
		for (Entry<Integer, Double> e : sortedList) {
			if (gtSet.contains(e.getKey()))	++tp;
			else AUC += (double) tp / (gtSet.size() * (sortedList.size() - gtSet.size())) ;
		}
		System.out.println("AUC = " + AUC);
	}
	
	private HashSet<Integer> readGt(String gtFileName) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(new File(gtFileName)));
		HashSet<Integer> gtSet = new HashSet<Integer>();
		String s;
		while ((s = br.readLine()) != null) {
			s = s.trim();
			if (s == "") continue;
			gtSet.add(rkdata.name2Id.get(s));
		}
		br.close();
		return gtSet;
	}
	
	private HashMap<Integer, Double> readGtScores(String gtScoreFileName) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(new File(gtScoreFileName)));
		HashMap<Integer, Double> ret = new HashMap<Integer, Double>();
		String s;
		while ((s = br.readLine()) != null) {
			s = s.trim();
			if (s == "") continue;
			String[] slist = s.split("\\s+");
			ret.put(rkdata.name2Id.get(slist[0]), Double.parseDouble(slist[1]));
		}
		br.close();
		return ret;
	}
	
	public void outputGtScoresAccordingToGivenOrder(String orderFileName, String outputFile) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(new File(orderFileName)));
		ArrayList<Integer> orderedIdList= new ArrayList<Integer>();
		String s;
		while ((s = br.readLine()) != null) {
			s = s.trim();
			if (s == "") continue;
			String[] slist = s.split("\\s+");
			orderedIdList.add(rkdata.name2Id.get(slist[0]));
		}
		br.close();
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputFile)));
		for (Integer id : orderedIdList) {
			if (id == null) 
				bw.write("-1\t" + 0.0 + "\n");
			else 
				bw.write(rkdata.id2Name.get(id) + "\t" + scores[id] + "\n");
		}
		bw.flush();
		bw.close();
	}
}
