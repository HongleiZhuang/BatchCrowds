package learner;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import util.MapSorter;
import model.SemiRankingDataSet;

public class PlackettLuceSemiMMLearner extends ScoreBasedSemiRankingLearner {

	int maxIter = 200;        // Maximum iteration 
	double epsilon = 1e-2;    // Stop criterion for convergence
	double tempLambda;  // Probability that a worker guesses the answer
	
	
	int curRow;
	double[][] tempScores; 
	ArrayList<ArrayList<ArrayList<int[]>>> permLists;
	
	
	double betaDistAlpha = 2, betaDistBeta = 2;
	
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
				
		for (int i = 0; i < tempScores[curRow].length; ++i) tempScores[curRow][i] = 1;
		
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

			
			for (int i = 0; i < scores.length; ++i) {
				double A = tempScores[nextRow][i];
				double B = -tempScores[nextRow][i] - w[i] - betaDistAlpha - betaDistBeta + 2;
				double C = w[i] + betaDistAlpha - 1;
//				System.out.println("i = " + i + ", A = " + A + ", B = " + B + ", C =" + C);
				if (A == 0) tempScores[nextRow][i] = -C / B;
				else tempScores[nextRow][i] = (-B - Math.sqrt(B * B - 4 * A * C)) / (2 * A);
//				tempScores[nextRow][i] = (double) w[i] / tempScores[nextRow][i];
			}
			curRow ^= 1;
			scores = tempScores[curRow];
//			normalizedByL1(scores);
			
			
			
			// ================Debug========================
			System.out.println("==================Iter = " + iter + "====================");
			for (int i = 0; i < scores.length; ++i) {
				name2Score.put(rkdata.id2Name.get(i), scores[i]);
			}
			ArrayList<Entry<String, Double>> sortedList = MapSorter.sortByValue(name2Score, false);
			int top = 0;
			for (Entry<String, Double> e: sortedList) {
				System.out.println(e.getKey() + ":" + e.getValue() + "\t");
//				System.out.print(e.getValue() + "\t");
//				if (++top >= 20) break;
			}
			System.out.println();
//			System.out.println("================");
			
		}
	}
	
	
	public void rescaleScores(double b, double rho, double tempLambda)  {
		this.tempLambda = tempLambda;
		double maxScore = getLInfinite(scores);
		double alphaUpperBound = 1.0 / maxScore;
		HashMap<Integer, Double> endPointsMap = new HashMap<Integer, Double>();
		endPointsMap.put(-1, alphaUpperBound);
		System.out.println(rkdata.semiRankingLists.size());
		for (int j = 0; j < rkdata.semiRankingLists.size(); ++j) {
			ArrayList<int[]> semiRankedList = rkdata.semiRankingLists.get(j);
			double u = semiRankedList.get(0).length == 0 ? 1.0 : getMin(getScoreVector(semiRankedList.get(0)));
			double l = semiRankedList.get(1).length == 0 ? 1.0 : getLInfinite(getScoreVector(semiRankedList.get(1)));
//			System.out.println("u=" + u +  ", l="  + l);
//			System.out.println("u_alpha=" + (0.5 / u) +  ", l_alpha="  + (0.5 / l));
			if (0.5 < u * alphaUpperBound) endPointsMap.put(j * 2, 0.5 / u);
			if (0.5 < l * alphaUpperBound) endPointsMap.put(j * 2 + 1, 0.5 / l);
		}
		System.out.println(endPointsMap.size());
		double maxL = - Double.MAX_VALUE, bestAlpha = alphaUpperBound;
		for (Entry<Integer, Double> e : endPointsMap.entrySet()) {
			double tempAlpha = e.getValue();
			System.out.println("when alpha = " + tempAlpha);
			double L = 0.0;
			for (int j = 0; j < rkdata.semiRankingLists.size(); ++j) {
				double lambda = this.tempLambda;
				ArrayList<int[]> semiRankedList = rkdata.semiRankingLists.get(j);
				double u = semiRankedList.get(0).length == 0 ? 1.0 : getMin(getScoreVector(semiRankedList.get(0)));
				double l = semiRankedList.get(1).length == 0 ? 1.0 : getLInfinite(getScoreVector(semiRankedList.get(1)));
				double p = getRandomPickingProb(b, rho, semiRankedList.get(0).length, semiRankedList.get(0).length + semiRankedList.get(1).length);
				if (0.5 >= l * tempAlpha && 0.5 < u * tempAlpha) p = p * (1 - lambda) + lambda;
				L += Math.log(p);
			}
			System.out.println("L = " + L);
			if (L > maxL) {
				bestAlpha = tempAlpha;
				maxL = L;
			}
		}
		for (int i = 0; i < scores.length; ++i) scores[i] *= bestAlpha;
		
		// ================Debug========================
		HashMap<String, Double> name2Score = new HashMap<String, Double>();
		for (int i = 0; i < scores.length; ++i) {
			name2Score.put(rkdata.id2Name.get(i), scores[i]);
		}
		ArrayList<Entry<String, Double>> sortedList = MapSorter.sortByValue(name2Score, false);
		int top = 0;
		for (Entry<String, Double> e: sortedList) {
			System.out.println(e.getKey() + ":" + e.getValue() + "\t");
//			if (++top >= 20) break;
		}
		System.out.println();
		// =============================================
	}
	
	private double getRandomPickingProb(double b, double pho, int k, int n) {
		//Assuming P(picking the k from [0,...,n]) = C(k + b)^(-pho) 
		double[] randomPickProb = new double[n + 1];
		for (int i = 0; i <= n; ++i) randomPickProb[i] = Math.exp(- pho * Math.log(i + b));
		normalizedByL1(randomPickProb);
		return randomPickProb[k];
	}
	
	private double[] getScoreVector(int[] idVector) {
		double[] ret = new double[idVector.length];
		for (int i = 0; i < idVector.length; ++i) ret[i] = scores[idVector[i]];
		return ret;
	}

	private double getMin(double[] v) {
		double min = Math.abs(v[0]);
		for (int i = 1; i < v.length; ++i) min = min < Math.abs(v[i]) ? min : Math.abs(v[i]);
		return min;		
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

	
	private void normalizedByL1(double[] v) {
		double sum = 0.0;
		for (double d : v) sum += Math.abs(d);
		if (sum > 0) for (int i = 0; i < v.length; ++i) v[i] /= sum;
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
		dataSet.readSemiRankingLists("C:\\Coursework\\CS598Aditya\\project\\crowdsource\\exp\\1013_try\\rankedlists_beta.txt");
		PlackettLuceSemiMMLearner learner = new PlackettLuceSemiMMLearner(dataSet);
		learner.trainRankings();
		learner.rescaleScores(1, 1, 0.8);
		learner.evaluate("C:\\Coursework\\CS598Aditya\\project\\crowdsource\\exp\\1013_try\\ground_truth.txt");
	}
	
}
