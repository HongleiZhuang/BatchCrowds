package learner;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import util.MapSorter;
import model.RankingDataSet;
import model.SemiRankingDataSet;

public class NaiveLearner extends ScoreBasedSemiRankingLearner {

	public double threshold;
	
	static public ScoreBasedSemiRankingLearner createNaiveLearner(SemiRankingDataSet rkdata) {
		return new NaiveLearner(rkdata);
	}
	
	private NaiveLearner(SemiRankingDataSet rkdata) {
		this.rkdata = rkdata;
		this.scores = new double[rkdata.id2Name.size()];
		this.threshold = 0.5;
	}
	
	@Override
	public ArrayList<Integer> getInferredPositiveIdList() {
		ArrayList<Integer> ret = new ArrayList<Integer>();
		for (int i = 0; i < scores.length; ++i) 
			if (scores[i] > this.threshold) ret.add(i);
		return ret;
	}
	
	void trainWorkerModelParams(SemiRankingDataSet trainRkdata, HashMap<Integer, Double> rawGroundTruthScores) {
		HashMap<Integer, Integer> pos = new HashMap<Integer, Integer>(), tot = new HashMap<Integer, Integer>();
		for (int j = 0; j < trainRkdata.semiRankingLists.size(); ++j) {
			ArrayList<int[]> semiRankList = trainRkdata.semiRankingLists.get(j);
			for (int id : semiRankList.get(0)) {
				pos.put(id, pos.containsKey(id) ? pos.get(id) + 1 : 0);
			}
			for (int[] sess : semiRankList) {
				for (int id : sess) {
					tot.put(id, tot.containsKey(id) ? tot.get(id) + 1 : 0);
				}
			}
		}
		HashMap<Integer, Double> trainScores = new HashMap<Integer, Double>();
		for (int id : pos.keySet()) trainScores.put(id, (double) pos.get(id) / tot.get(id));
		ArrayList<Entry<Integer, Double>> sortedList = MapSorter.sortByValue(trainScores, false);
		HashSet<Integer> posIdSet = new HashSet<Integer>();
		for (Entry<Integer, Double> e : rawGroundTruthScores.entrySet()) if (e.getValue() > 0.5) posIdSet.add(e.getKey());
		
		int tp = 0, fp = 0, fn = posIdSet.size();
		double maxF1 = 0.0;
		for (Entry<Integer, Double> e : sortedList) {
			if (posIdSet.contains(e.getKey())) {
				++tp;
				--fn;
			}
			else ++fp;
			double prc = (double) tp / (tp + fp);
			double rcl = (double) tp / (tp + fn);
			double f1 = 2 * prc * rcl / (prc + rcl);
			if (maxF1 < f1) {
				maxF1 = f1;
				threshold = e.getValue();
			}
		}
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
		/*
		SemiRankingDataSet trainDataSet = new SemiRankingDataSet();
		trainDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/job_509470.json.train.srk");
		HashMap<Integer, Double> trainScores = trainDataSet.readGtScores("/Users/hzhuang/Work/beta/ranking/data/fullres.out.train.filtered");
		System.out.println(trainDataSet.id2Name.size());
		
		
		SemiRankingDataSet testDataSet = new SemiRankingDataSet();
		testDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/job_509470.json.test.srk");
		System.out.println(testDataSet.id2Name.size());

		
		NaiveLearner learner = new NaiveLearner(testDataSet);
		learner.trainWorkerModelParams(trainDataSet, trainScores);
		
		learner.trainRankings();
		
		learner.evaluate("/Users/hzhuang/Work/beta/ranking/data/fullres.out.test.filtered");
		learner.evaluateByROC("/Users/hzhuang/Work/beta/ranking/data/fullres.out.test.filtered");
		/**/
		
		
		SemiRankingDataSet trainDataSet = new SemiRankingDataSet();
		trainDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/synthetic/rankedlists.txt");
		HashMap<Integer, Double> trainScores = trainDataSet.readGtScores("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_score.txt");
		System.out.println(trainDataSet.id2Name.size());
		
		
		SemiRankingDataSet testDataSet = new SemiRankingDataSet();
		testDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/synthetic/rankedlists_test.txt");
		System.out.println(testDataSet.id2Name.size());

		
		NaiveLearner learner = new NaiveLearner(testDataSet);
		learner.trainWorkerModelParams(trainDataSet, trainScores);
		System.out.println("Threshold = " + learner.threshold);
		
		learner.trainRankings();
		
		learner.evaluate("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label_test.txt");
		learner.evaluateByROC("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label_test.txt");
		/**/
		
		/*
		SemiRankingDataSet dataSet = new SemiRankingDataSet();
		dataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/synthetic/rankedlists_test.txt");
		ScoreBasedSemiRankingLearner learner = NaiveLearner.createNaiveLearner(dataSet);
		learner.trainRankings();
		learner.evaluate("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label_test.txt");
		learner.evaluateByROC("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label_test.txt");
		/**/
		
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
