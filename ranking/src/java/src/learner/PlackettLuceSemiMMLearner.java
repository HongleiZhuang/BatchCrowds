package learner;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import util.MapSorter;
import model.SemiRankingDataSet;

public class PlackettLuceSemiMMLearner extends ScoreBasedSemiRankingLearner {

	int maxIter = 100;      // Maximum iteration 
	double epsilon = 1e-2;  // Stop criterion for convergence
	
	int curRow;
	double[][] tempScores; 
	ArrayList<ArrayList<ArrayList<int[]>>> permLists;
	
	static public ScoreBasedSemiRankingLearner createPlacketBasedSemiRankingLearner(SemiRankingDataSet rkdata) {
		return new PlackettLuceSemiMMLearner(rkdata);
	}
	
	private PlackettLuceSemiMMLearner(SemiRankingDataSet rkdata) {
		this.rkdata = rkdata;
		this.tempScores = new double[2][rkdata.id2Name.size()];
		this.curRow = 0;
		this.scores = tempScores[curRow];
	}
	
	@Override
	void trainRankings() {
		double[] w = new double[scores.length];
		for (ArrayList<int[]> rankedList : rkdata.semiRankingLists) {
			for (int s = 0; s < rankedList.size() - 1; ++s) {
				for (int k : rankedList.get(s)) w[k] += 1.0;
			}
		}
				
		for (int i = 0; i < tempScores[curRow].length; ++i) tempScores[curRow][i] = 0.001;
		
		permLists = new ArrayList<ArrayList<ArrayList<int[]>>>();
		for (int j = 0; j < rkdata.semiRankingLists.size(); ++j) {
			ArrayList<int[]> semiRankedList = rkdata.semiRankingLists.get(j);
			ArrayList<ArrayList<int[]>> semiPermList = new ArrayList<ArrayList<int[]>>();
			for (int s = 0; s < semiRankedList.size() - 1; ++s) {  //No need to generate permutations for the last session
				if (semiRankedList.get(s).length == 0) {
					semiPermList.add(new ArrayList<int[]>());
					continue;
				}
				ArrayList<int[]> perms = getPermutation(semiRankedList.get(s));  
				semiPermList.add(perms);
			}
			permLists.add(semiPermList);
		}
		
		HashMap<String, Double> name2Score = new HashMap<String, Double>();
		
		for (int iter = 0; iter < maxIter; ++iter) {
			int nextRow = curRow ^ 1;
			for (int i = 0; i < scores.length; ++i) tempScores[nextRow][i] = 0.0;
			
			for (int j = 0; j < rkdata.semiRankingLists.size(); ++j) {
				ArrayList<int[]> semiRankedList = rkdata.semiRankingLists.get(j);
				
				double[] suffixSessionSum = new double[semiRankedList.size()];
				for (int s = 0; s < semiRankedList.size(); ++s) {
					for (int id : semiRankedList.get(s)) {
						suffixSessionSum[s] += tempScores[curRow][id]; 
					}
				}
				for (int s = suffixSessionSum.length - 2; s >= 0; --s) {
					suffixSessionSum[s] += suffixSessionSum[s + 1];
				}
				
				double prefixSum = 0.0;
				for (int s = 0; s < semiRankedList.size(); ++s) {
					for (int id : semiRankedList.get(s)) tempScores[nextRow][id] += prefixSum;
					if (s == semiRankedList.size() - 1) break;
					
					double sessionSum = 0.0;
					double[] suffixSubsessionSum = new double[semiRankedList.get(s).length];
					int cnt = permLists.get(j).get(s).size();
					if (cnt == 0) continue;					
					for (int[] perm : permLists.get(j).get(s)) {
						
						suffixSubsessionSum[suffixSubsessionSum.length - 1] = tempScores[curRow][perm[perm.length - 1]];
						for (int k = suffixSubsessionSum.length - 2; k >= 0; --k) suffixSubsessionSum[k] = suffixSubsessionSum[k + 1] + tempScores[curRow][perm[k]];
						
						double permPrefixSum = 0.0;
						for (int k = 0; k < perm.length; ++k) {
							permPrefixSum += (double) 1.0 / (suffixSubsessionSum[k] + suffixSessionSum[s + 1]);
							tempScores[nextRow][perm[k]] += (double) permPrefixSum / cnt;
						}
						sessionSum += (double)  permPrefixSum / cnt;
					}
					prefixSum += sessionSum;
				}
				
 			}

			
			for (int i = 0; i < scores.length; ++i) tempScores[nextRow][i] = (double) w[i] / tempScores[nextRow][i];
			curRow ^= 1;
			scores = tempScores[curRow];
//			normalizedByL1(scores);
			
			
			
			// ================Debug========================
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
	}
	
	private void normalizedByL1(double[] v) {
		double sum = 0.0;
		for (double d : v) sum += Math.abs(d);
		if (sum > 0) 
			for (int i = 0; i < v.length; ++i) v[i] /= sum;
	}
	
	ArrayList<int[]> getPermutation(int[] elemList) {
		ArrayList<int[]> ret = new ArrayList<int[]>();
		generatePermutation(elemList, 0, ret);
		return ret;
	}
	
	private void generatePermutation(int[] elemList, int cur, ArrayList<int[]> ret) {
		if (cur == elemList.length - 1) {
			ret.add(elemList.clone());
			return;
		}
		int curElem = elemList[cur];
		for (int i = cur; i < elemList.length; ++i) {
			elemList[cur] = elemList[i];
			elemList[i] = curElem;
			generatePermutation(elemList, cur + 1, ret);
			elemList[i] = elemList[cur];
		}
		elemList[cur] = curElem;
	}
	
	static public void main(String[] args) throws Exception {
		SemiRankingDataSet dataSet = new SemiRankingDataSet();
		dataSet.readSemiRankingLists("C:\\Coursework\\CS598Aditya\\project\\crowdsource\\exp\\1013_try\\rankedlists_bin_200.txt");
		PlackettLuceSemiMMLearner learner = new PlackettLuceSemiMMLearner(dataSet);
		learner.trainRankings();
	}
	

}
