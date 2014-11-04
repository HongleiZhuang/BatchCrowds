package learner;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import util.MapSorter;
import model.RankingDataSet;
import model.SemiRankingDataSet;

public class PlackettLuceMMLearner extends ScoreBasedRankingLearner {

	int maxIter = 2000;      // Maximum iteration 
	double epsilon = 1e-2;  // Stop criterion for convergence
	
	
	int curRow;
	double[][] tempScores; 
	
	
	static public ScoreBasedRankingLearner createPlackettLuceMMLearner(RankingDataSet rkdata) {
		return new PlackettLuceMMLearner(rkdata);
	}
	
	private PlackettLuceMMLearner(RankingDataSet rkdata) {
		this.rkdata = rkdata;
		this.tempScores = new double[2][rkdata.id2Name.size()];
		this.curRow = 0;
		this.scores = tempScores[curRow];
	}
	
	@Override
	void trainRankings() {
		double[] w = new double[scores.length];
		for (ArrayList<Integer> rankedList : rkdata.rankingLists) {
			for (int k = 0; k < rankedList.size() - 1; ++k) w[rankedList.get(k)] += 1.0;
		}
		
		HashMap<String, Double> name2Score = new HashMap<String, Double>();
		
		for (int i = 0; i < tempScores[curRow].length; ++i) tempScores[curRow][i] = 1.0;
		for (int iter = 0; iter < maxIter; ++iter) {
			int nextRow = curRow ^ 1;
			for (int i = 0; i < scores.length; ++i) tempScores[nextRow][i] = 0.0;
			for (int j = 0; j < rkdata.rankingLists.size(); ++j) {
				ArrayList<Integer> rankedList = rkdata.rankingLists.get(j);
				double[] suffixSum = new double[rankedList.size()];
				suffixSum[rankedList.size() - 1] = tempScores[curRow][rankedList.get(rankedList.size() - 1)]; 
				for (int k = rankedList.size() - 2; k >= 0 ; --k) {
					suffixSum[k] = suffixSum[k + 1] + tempScores[curRow][rankedList.get(k)];
				}
				double prefixSum = 0.0;
				for (int k = 0; k < rankedList.size(); ++k) {
					prefixSum += (double) 1.0 / suffixSum[k];
					tempScores[nextRow][rankedList.get(k)] += prefixSum;
				}
			}
			for (int i = 0; i < scores.length; ++i) tempScores[nextRow][i] = (double) w[i] / tempScores[nextRow][i];
			curRow ^= 1;
			scores = tempScores[curRow];
			
			
			for (int i = 0; i < scores.length; ++i) {
				name2Score.put(rkdata.id2Name.get(i), scores[i]);
			}
			ArrayList<Entry<String, Double>> sortedList = MapSorter.sortByValue(name2Score, false);
			int top = 0;
			for (Entry<String, Double> e: sortedList) {
				System.out.println(e.getKey() + ":" + e.getValue() + "\t");
//				System.out.print(e.getValue() + "\t");
				if (++top >= 20) break;
			}
			System.out.println();
//			System.out.println("================");
		}
		
		// ================Debug========================
		normalizedByLInfinite(scores);
		for (int i = 0; i < scores.length; ++i) {
			name2Score.put(rkdata.id2Name.get(i), scores[i]);
		}
		ArrayList<Entry<String, Double>> sortedList = MapSorter.sortByValue(name2Score, false);
		int top = 0;
		for (Entry<String, Double> e: sortedList) {
			System.out.println(e.getKey() + ":" + e.getValue() + "\t");
			if (++top >= 20) break;
		}
		System.out.println();
		// =============================================
	}
	
	private double getLInfinite(double[] v) {
		double max = Math.abs(v[0]);
		for (int i = 1; i < v.length; ++i) max = max > Math.abs(v[i]) ? max : Math.abs(v[i]);
		return max;
	}
	
	private void normalizedByLInfinite(double[] v) {
		double lInfinite = getLInfinite(v);
		if (lInfinite > 0) for (int i = 0; i < v.length; ++i) v[i] /= lInfinite;
	}
	
	
	static public void main(String[] args) throws Exception {
//		RankingDataSet rkdata = new RankingDataSet();
//		rkdata.readRankingLists("C:\\Coursework\\CS598Aditya\\project\\crowdsource\\exp\\1013_try\\rankedlists_bin_200.txt");
		SemiRankingDataSet dataSet = new SemiRankingDataSet();
		dataSet.readSemiRankingLists("C:\\Coursework\\CS598Aditya\\project\\crowdsource\\exp\\1013_try\\rankedlists_beta.txt");
		ScoreBasedRankingLearner learner = PlackettLuceMMLearner.createPlackettLuceMMLearner(dataSet.cloneAsRankingDataSet());
		learner.trainRankings();
		
		
	}

}
