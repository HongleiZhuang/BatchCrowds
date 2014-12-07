package learner;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import util.MapSorter;
import model.RankingDataSet;
import model.SemiRankingDataSet;

public class NaiveLearner extends ScoreBasedSemiRankingLearner {

	
	static public ScoreBasedSemiRankingLearner createNaiveLearner(SemiRankingDataSet rkdata) {
		return new NaiveLearner(rkdata);
	}
	
	private NaiveLearner(SemiRankingDataSet rkdata) {
		this.rkdata = rkdata;
		this.scores = new double[rkdata.id2Name.size()];
	}
	
	@Override
	void trainRankings() {
		int[] poscnt = new int[scores.length];
		int[] totcnt = new int[scores.length];
		for (ArrayList<int[]> rankedList : rkdata.semiRankingLists) {
			for (int k : rankedList.get(0)) ++poscnt[k];
			for (int[] l : rankedList) {
				for (int k : l) ++totcnt[k];
			}
		}
		
		for (int i = 0; i < scores.length; ++i) {
			scores[i] = (double) poscnt[i] / totcnt[i];
			System.out.println(poscnt[i] + "/" + totcnt[i]);
		}
		
		
		HashMap<String, Double> name2Score = new HashMap<String, Double>();
		for (int i = 0; i < scores.length; ++i) {
			name2Score.put(rkdata.id2Name.get(i), scores[i]);
		}
		ArrayList<Entry<String, Double>> sortedList = MapSorter.sortByValue(name2Score, false);
		int top = 0;
		for (Entry<String, Double> e: sortedList) {
			System.out.println(e.getKey() + ":" + e.getValue() + "\t");
//			System.out.print(e.getValue() + "\t");
//			if (++top >= 20) break;
		}
		System.out.println();
	}

	static public void main(String[] args) throws Exception {
		SemiRankingDataSet dataSet = new SemiRankingDataSet();
		dataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/synthetic/rankedlists.txt");
		ScoreBasedSemiRankingLearner learner = NaiveLearner.createNaiveLearner(dataSet);
		learner.trainRankings();
		learner.evaluate("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label.txt");
		learner.evaluateByROC("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label.txt");

		
		/*
		SemiRankingDataSet dataSet = new SemiRankingDataSet();
		dataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/job_509470.json.test.srk");
		ScoreBasedSemiRankingLearner learner = NaiveLearner.createNaiveLearner(dataSet);
		learner.trainRankings();
		learner.evaluate("/Users/hzhuang/Work/beta/ranking/data/fullres.out.test.filtered");
		learner.evaluateByROC("/Users/hzhuang/Work/beta/ranking/data/fullres.out.test.filtered");
		/**/
	}
	
}
