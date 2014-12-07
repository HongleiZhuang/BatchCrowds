package learner;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import model.SemiRankingDataSet;
import util.MapSorter;

public class CollectiveWorkerModel extends ScoreBasedSemiRankingLearner {

	int maxIter = 2000;        // Maximum iteration 
	double epsilon = 1e-3;    // Stop criterion for convergence
	double tempLambda;  // Probability that a worker guesses the answer
	
	int trainWorkerMaxIter = 100;
	double deltaProb = 1e-2;  //A small probability prevent the learned parameters to be zero.
	
	double dampingFactor = 0.5;
	int curRow;
	double[][] tempScores; 
	ArrayList<ArrayList<ArrayList<int[]>>> permLists;
	double[] piProbs;
	
	double[] pickingProb;
	double   correctProb;
	
	
	double betaDistAlpha = 1.5, betaDistBeta = 1.5;
	
	
	
	
	private CollectiveWorkerModel(SemiRankingDataSet rkdata) {
		this.rkdata = rkdata;
		this.tempScores = new double[2][rkdata.id2Name.size()];
		this.curRow = 0;
		this.scores = tempScores[curRow];
		this.piProbs = new double[rkdata.semiRankingLists.size()];
	}
	

	
	@Override
	void trainRankings() {

		for (int i = 0; i < tempScores[curRow].length; ++i) tempScores[curRow][i] = 0.1;
		int maxLength = 0;
		for (int j = 0; j < rkdata.semiRankingLists.size(); ++j) {
			ArrayList<int[]> semiRankList = rkdata.semiRankingLists.get(j);
			// Get maximum semi-ranking list length
			int l = 0; 
			for (int[] sess : semiRankList) l += sess.length;
			if (l > maxLength) maxLength = l;
		}

		double[][] pickingProbHat = new double[2][maxLength + 1];
		double[] correctProbHat = new double[2];
		for (int i = 0; i < maxLength + 1; ++i) pickingProbHat[curRow][i] = Math.random();
		normalizedByL1(pickingProbHat[curRow]);
		correctProbHat[curRow] = 0.5;
		this.correctProb = correctProbHat[curRow];
		this.pickingProb = pickingProbHat[curRow];

		
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
		
		
		double oldL = -Double.MAX_VALUE;
		int stopCnt = 0;
		for (int iter = 0; iter < maxIter; ++iter) {
			int nextRow = curRow ^ 1;
			for (int i = 0; i < scores.length; ++i) tempScores[nextRow][i] = 0.0;
			calcPiProbs(tempScores[curRow]);
			
			double[] w = new double[scores.length];
			double[] negW = new double[scores.length];
			
			for (int j = 0; j < rkdata.semiRankingLists.size(); ++j) {				
				ArrayList<int[]> semiRankedList = rkdata.semiRankingLists.get(j);
//				double correctFlag = 1.0;
//				for (int id : semiRankedList.get(0)) if (tempScores[curRow][id] < 0.5) correctFlag = 0.0;
//				for (int id : semiRankedList.get(1)) if (tempScores[curRow][id] >= 0.5) correctFlag = 0.0;
				double correctFlag = 0.0;
				for (int id : semiRankedList.get(0)) {
					correctFlag += Math.log(tempScores[curRow][id]);
				}
				for (int id : semiRankedList.get(1)) {
					correctFlag += Math.log(1 - tempScores[curRow][id]);
				}
				
//				double pZX = (correctProb * correctFlag) / (correctProb * correctFlag + (1 - correctProb) * piProbs[j] * pickingProb[semiRankedList.get(0).length]);
				double pZX = (correctProb * Math.exp(correctFlag)) / (correctProb * Math.exp(correctFlag) + (1 - correctProb) * piProbs[j] * pickingProb[semiRankedList.get(0).length]);
				correctProbHat[nextRow] += pZX;
				pickingProbHat[nextRow][semiRankedList.get(0).length] += 1 - pZX;

				
				for (int s = 0; s < semiRankedList.size() - 1; ++s) {
//					for (int k : semiRankedList.get(s)) w[k] += 1 - pZX;
					for (int k : semiRankedList.get(s)) w[k] += 1 - pZX + pZX ;
				}
				for (int k : semiRankedList.get(semiRankedList.size() - 1)) negW[k] += pZX; 
				
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
					for (int id : semiRankedList.get(s)) tempScores[nextRow][id] += (1 - pZX) * prefixSum;
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
							tempScores[nextRow][perm[k]] += (1 - pZX)  * (double) permPrefixSum / cnt;
						}
						sessionSum += (double)  permPrefixSum / cnt;
					}
					prefixSum += sessionSum;
				}
				
 			}

			correctProbHat[nextRow] /= rkdata.semiRankingLists.size();
			normalizedByL1(pickingProbHat[nextRow]);

			for (int i = 0; i < scores.length; ++i) {
				
				double A = tempScores[nextRow][i];
//				double B = -tempScores[nextRow][i] - w[i] - betaDistAlpha - betaDistBeta + 2;
				double B = -tempScores[nextRow][i] - w[i] - negW[i] - betaDistAlpha - betaDistBeta + 2;
				double C = w[i] + betaDistAlpha - 1;
				
				if (A == 0) tempScores[nextRow][i] = tempScores[curRow][i] * (dampingFactor) + (1 - dampingFactor) * (-C / B);
				else {
					tempScores[nextRow][i] = tempScores[curRow][i] * (dampingFactor)
							+ (1 - dampingFactor) * ((-B - Math.sqrt(B * B - 4 * A * C)) / (2 * A));
//					if (tempScores[nextRow][i] < 0) tempScores[nextRow][i] =  (-B + Math.sqrt(B * B - 4 * A * C)) / (2 * A);
				}
//				System.out.println("i = " + i + ", A = " + A + ", B = " + B + ", C =" + C + ", B^2-4AC=" + (B * B - 4 * A * C) + ", score=" + tempScores[nextRow][i]);
 				/**/ 
				//tempScores[nextRow][i] = (double) w[i] / tempScores[nextRow][i];
			}
			
			curRow ^= 1;
			scores = tempScores[curRow];
			this.correctProb = correctProbHat[curRow];
			this.pickingProb = pickingProbHat[curRow];
			
			
			
			
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
			/**/
			
			System.out.println("lambda = " + correctProb);
			System.out.println("picking prob = ");
			for (int i = 0; i < maxLength + 1; ++i) System.out.print(pickingProb[i] + "\t");
			System.out.println();

			
//			calcLogLikelihood();
			double L = calcLogLikelihood();
//			for (int i = 0; i < piProbs.length; ++i) L += Math.log(piProbs[i]);
			
			System.out.println("L = " + L);
			if (L < oldL || Math.abs((L - oldL) / oldL) < epsilon) {
				++stopCnt;
				oldL = L;
				if (stopCnt > 5) {
					System.out.println("Converge");
					break;
				}
			}
			else {
				stopCnt = 0;
				oldL = L;
			}
		}
		
//		normalizedByLInfinite(scores);
	}
	
	
	public double calcLogLikelihood() {
		double L = 0.0;
		calcPiProbs(tempScores[curRow]);
		for (int j = 0; j < rkdata.semiRankingLists.size(); ++j) {
			ArrayList<int[]> semiRankedList = rkdata.semiRankingLists.get(j);
			double correctFlag = 1.0;
			for (int id : semiRankedList.get(0)) if (tempScores[curRow][id] < 0.5) correctFlag = 0.0;
			for (int id : semiRankedList.get(1)) if (tempScores[curRow][id] >= 0.5) correctFlag = 0.0;
			double prob = (1 - correctProb) * piProbs[j] * pickingProb[semiRankedList.get(0).length];
			if (correctFlag == 1.0) prob += correctFlag;
			L += Math.log(prob);
		}
		return L;
	}
	
	public void calcPiProbs(double[] scores) {
		for (int recId = 0; recId <  rkdata.semiRankingLists.size(); ++recId) {
			ArrayList<int[]> semiRankingList = rkdata.semiRankingLists.get(recId);
			int cnt = 0;
			for (int j = 0; j < semiRankingList.size(); ++j) cnt += semiRankingList.get(j).length;
			int[] fullList = new int[cnt];
			int k = semiRankingList.get(0).length;
			for (int j = 1; j < semiRankingList.size(); ++j) 
				for (int l = 0; l < semiRankingList.get(j).length; ++l) fullList[k++] = semiRankingList.get(j)[l];
			piProbs[recId] = 0.0;
			if (permLists.get(recId).get(0).size() == 0) piProbs[recId] = 1;
			else {
				for (int[] permList : permLists.get(recId).get(0)) {
					for (int j = 0; j < permList.length; ++j) fullList[j] = permList[j];
					double prob = 0.0;
	//				for (int i : fullList) System.out.print(i + ",");
	//				System.out.println();
					for (int i = 0; i < permList.length; ++i) {
						double sum = 0.0;
						for (int j = i; j < fullList.length; ++j)
							sum += scores[fullList[j]];
	//					System.out.println("score=" + scores[fullList[i]] + ", sum=" + sum);
						prob += Math.log(scores[fullList[i]]) - Math.log(sum);
					}
					piProbs[recId] += Math.exp(prob);
				}
			}
//			System.out.println(likelihoods[recId]);
		}
	}
	
	
	
	public double getProb(int[] orderedList) {
		double ret = 1.0;
		for (int i = 0; i < orderedList.length; ++i) {
			double sum = 0;
			for (int j = i; j < orderedList.length; ++j)
				sum += scores[orderedList[j]];
			ret *= scores[orderedList[i]] / sum;
		}
		return ret;
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
		/*
		SemiRankingDataSet trainDataSet = new SemiRankingDataSet();
//		trainDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/job_509470.json.train.srk");
		trainDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/synthetic/rankedlists.txt");
		System.out.println(trainDataSet.id2Name.size());
//		HashMap<Integer, Integer> trainLabels = trainDataSet.readGtLabel("/Users/hzhuang/Work/beta/ranking/data/fullres.out.train.filtered");
		HashMap<Integer, Integer> trainLabels = trainDataSet.readGtLabel("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label.txt");
		/**/
		
		SemiRankingDataSet testDataSet = new SemiRankingDataSet();
		testDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/job_509470.json.test.srk");
		System.out.println(testDataSet.id2Name.size());
		/**/
		
		
//		SimpleWorkerModelLearner learner = new SimpleWorkerModelLearner(testDataSet);
		CollectiveWorkerModel learner = new CollectiveWorkerModel(testDataSet);
		
//		learner.trainWorkerModelParams(trainDataSet, trainLabels);
		learner.trainRankings();
		learner.evaluate("/Users/hzhuang/Work/beta/ranking/data/fullres.out.test.filtered");
		learner.evaluateByROC("/Users/hzhuang/Work/beta/ranking/data/fullres.out.test.filtered");
	}

}
