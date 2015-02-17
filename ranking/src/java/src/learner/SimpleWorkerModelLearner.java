package learner;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import util.MapSorter;
import model.SemiRankingDataSet;

public class SimpleWorkerModelLearner extends ScoreBasedSemiRankingLearner {
	
	int maxIter = 2000;        // Maximum iteration 
	double epsilon = 1e-3;    // Stop criterion for convergence
	double tempLambda;  // Probability that a worker guesses the answer
	
	int trainWorkerMaxIter = 100;
	double deltaProb = 1e-2;  //A small probability prevent the learned parameters to be zero.
	double deltaScore = 1e-5; //A small constant prevent the scores to be 0 or 1. 
	
	double dampingFactor = 0;
	int curRow;
	double[][] tempScores; 
	ArrayList<ArrayList<ArrayList<int[]>>> permLists;
	double[] piProbs;
	
	double[] pickingProb;
	double   correctProb;
	
	
	double betaDistAlpha = 2, betaDistBeta = 2;
	
	String logFileString = null;
	
	static public ScoreBasedSemiRankingLearner createSimpleWorkerModelLearner(SemiRankingDataSet rkdata, SemiRankingDataSet trainRkdata, HashMap<Integer, Double> trainGtScores) throws IOException {
		return new SimpleWorkerModelLearner(rkdata, trainRkdata, trainGtScores);
	}
	
	private SimpleWorkerModelLearner(SemiRankingDataSet rkdata, SemiRankingDataSet trainRkdata, HashMap<Integer, Double> trainGtScores) throws IOException {
		this.rkdata = rkdata;
		this.tempScores = new double[2][rkdata.id2Name.size()];
		this.curRow = 0;
		this.scores = tempScores[curRow];
		this.piProbs = new double[rkdata.semiRankingLists.size()];
		this.trainWorkerModelParams(trainRkdata, trainGtScores);
	}
	
	public void saveModel(String modelFileString) throws Exception {
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(modelFileString)));
		bw.write(this.correctProb + "\n");
		bw.write(this.pickingProb.length + "\n");
		for (double ptau : this.pickingProb) bw.write(ptau + "\n");
		bw.flush();
		bw.close();
	}
	
	
	void trainWorkerModelParams(SemiRankingDataSet trainRkdata, HashMap<Integer, Double> rawGroundTruthScores) throws IOException {
		int maxLength = 0;
		double[] correctFlag = new double[trainRkdata.semiRankingLists.size()];
		double[] plProb = new double[trainRkdata.semiRankingLists.size()];
		HashMap<Integer, Double> groundTruthScores = new HashMap<Integer, Double>();
		for (Entry<Integer, Double> e : rawGroundTruthScores.entrySet()) {
			groundTruthScores.put(e.getKey(), (e.getValue() + deltaScore) / (1 + 2 * deltaScore));
		}
		for (int j = 0; j < trainRkdata.semiRankingLists.size(); ++j) {
			ArrayList<int[]> semiRankList = trainRkdata.semiRankingLists.get(j);
			// Get maximum semi-ranking list length
			int l = 0;
			double sumScores = 0.0;
			for (int[] sess : semiRankList) {
				l += sess.length;
				for (int id : sess) sumScores += groundTruthScores.get(id);
			}
			if (l > maxLength) maxLength = l;
			// Initialize whether a semi-ranking list is exactly coherent with the ground truth
			for (int id : semiRankList.get(0)) {
				correctFlag[j] += Math.log(groundTruthScores.get(id));
			}
			for (int id : semiRankList.get(1)) {
				correctFlag[j] += Math.log(1 - groundTruthScores.get(id));
			}
			
			if (semiRankList.get(0).length == 0) plProb[j] = 1.0;
			else {
				ArrayList<int[]> perms = getPermutation(semiRankList.get(0));
				for (int[] sess : perms) {
					double lProb = 0.0;
					double tempSumScores = sumScores;
					for (int t = 0; t < sess.length; ++t) {
						lProb += Math.log(groundTruthScores.get(sess[t])) - Math.log(tempSumScores);
						tempSumScores -= groundTruthScores.get(sess[t]);
					}
					plProb[j] += Math.exp(lProb);
				}
			}
		}
		
		// Initialize worker model parameters to learn
		double[][] pickingProbHat = new double[2][maxLength + 1];
		double[] correctProbHat = new double[2];
		for (int i = 0; i < maxLength + 1; ++i) pickingProbHat[0][i] = Math.random();
		normalizedByL1(pickingProbHat[0]);
		correctProbHat[0] = 0.5;
		int cur = 0;
		this.correctProb = correctProbHat[cur];
		this.pickingProb = pickingProbHat[cur];
		
		
		BufferedWriter logWriter = null;
		if (logFileString != null) {
			logWriter = new BufferedWriter(new FileWriter(new File(logFileString)));
		}
		
		for (int iter = 0; iter < trainWorkerMaxIter; ++iter) {
			int next = cur ^ 1;
			correctProbHat[next] = 0.0; // Temp variable for calculating next "correctProb"
			for (int i = 0; i < maxLength + 1; ++i) pickingProbHat[next][i] = deltaProb; // Temp variable for calculating next "pickingProb"
			for (int j = 0; j < trainRkdata.semiRankingLists.size(); ++j) {
				ArrayList<int[]> semiRankList = trainRkdata.semiRankingLists.get(j);
				double pZX = correctProb * Math.exp(correctFlag[j]) / (correctProb * Math.exp(correctFlag[j]) + (1 - correctProb) * plProb[j] * pickingProb[semiRankList.get(0).length]);
				correctProbHat[next] += pZX;
				pickingProbHat[next][semiRankList.get(0).length] += 1 - pZX;
			}
			correctProbHat[next] /= trainRkdata.semiRankingLists.size();
			normalizedByL1(pickingProbHat[next]);
			cur ^= 1;
			this.correctProb = correctProbHat[cur];
			this.pickingProb = pickingProbHat[cur];
//			System.out.println("lambda = " + correctProb);
//			System.out.println("picking prob = ");
//			for (int i = 0; i < maxLength + 1; ++i) System.out.print(pickingProb[i] + "\t");
//			System.out.println();
//			
			double L = 0.0;
			for (int j = 0; j < trainRkdata.semiRankingLists.size(); ++j) {
				ArrayList<int[]> semiRankList = trainRkdata.semiRankingLists.get(j);
				L += Math.log(correctProb * Math.exp(correctFlag[j]) + (1 - correctProb) * plProb[j] * pickingProb[semiRankList.get(0).length]);
			}
//			System.out.println("L = " + L);
//			if (logFileString != null) {
//				logWriter.write(L + "\n");
//				logWriter.flush();
//			}
		}
		
		System.out.println("lambda = " + correctProb);
		System.out.println("picking prob = ");
		for (int i = 0; i < maxLength + 1; ++i) System.out.print(pickingProb[i] + "\t");
		System.out.println();
		
		
		if (logFileString != null) logWriter.close();
	}
	
	/*
	void trainWorkerModelParams(SemiRankingDataSet trainRkdata, HashMap<Integer, Integer> groundTruth) {
		int maxLength = 0;
		double[] correctFlag = new double[trainRkdata.semiRankingLists.size()];
		for (int j = 0; j < trainRkdata.semiRankingLists.size(); ++j) {
			ArrayList<int[]> semiRankList = trainRkdata.semiRankingLists.get(j);
			// Get maximum semi-ranking list length
			int l = 0; 
			for (int[] sess : semiRankList) l += sess.length;
			if (l > maxLength) maxLength = l;
			// Initialize whether a semi-ranking list is exactly coherent with the ground truth 
			correctFlag[j] = 1.0;
			for (int id : semiRankList.get(0)) if (groundTruth.get(id) == 0) correctFlag[j] = 0.0;
			for (int id : semiRankList.get(1)) if (groundTruth.get(id) == 1) correctFlag[j] = 0.0;
		}
		
		// Initialize worker model parameters to learn
		double[][] pickingProbHat = new double[2][maxLength + 1];
		double[] correctProbHat = new double[2];
		for (int i = 0; i < maxLength + 1; ++i) pickingProbHat[0][i] = Math.random();
		normalizedByL1(pickingProbHat[0]);
		correctProbHat[0] = 0.5;
		int cur = 0;
		this.correctProb = correctProbHat[cur];
		this.pickingProb = pickingProbHat[cur];
		
		for (int iter = 0; iter < trainWorkerMaxIter; ++iter) {
			int next = cur ^ 1;
			correctProbHat[next] = 0.0; // Temp variable for calculating next "correctProb"
			for (int i = 0; i < maxLength + 1; ++i) pickingProbHat[next][i] = deltaProb; // Temp variable for calculating next "pickingProb"
			for (int j = 0; j < trainRkdata.semiRankingLists.size(); ++j) {
				ArrayList<int[]> semiRankList = trainRkdata.semiRankingLists.get(j);
				// Notice that we are assuming the probability of generating a correct ranking here is 1.
				// Because we do not actually know the eta value for ground-truth data.  
				double pZX = correctProb * correctFlag[j] / (correctProb * correctFlag[j] + (1 - correctProb) * pickingProb[semiRankList.get(0).length]);
				correctProbHat[next] += pZX;
				pickingProbHat[next][semiRankList.get(0).length] += 1 - pZX;
			}
			correctProbHat[next] /= trainRkdata.semiRankingLists.size();
			normalizedByL1(pickingProbHat[next]);
			cur ^= 1;
			this.correctProb = correctProbHat[cur];
			this.pickingProb = pickingProbHat[cur];
			System.out.println("lambda = " + correctProb);
			System.out.println("picking prob = ");
			for (int i = 0; i < maxLength + 1; ++i) System.out.print(pickingProb[i] + "\t");
			System.out.println();
		}
	}
	/**/

	
	@Override
	void trainRankings() {

		for (int i = 0; i < tempScores[curRow].length; ++i) tempScores[curRow][i] = 0.1;
		
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
					
					
					double[] conditional = new double[permLists.get(j).get(s).size()];
					for (int permId = 0; permId < permLists.get(j).get(s).size(); ++permId) {
						int[] perm = permLists.get(j).get(s).get(permId);
						suffixSubsessionSum[suffixSubsessionSum.length - 1] = tempScores[curRow][perm[perm.length - 1]];
						for (int k = suffixSubsessionSum.length - 2; k >= 0; --k) 	suffixSubsessionSum[k] = suffixSubsessionSum[k + 1] + tempScores[curRow][perm[k]];
						for (int k = 0; k < perm.length; ++k) conditional[permId] -= Math.log(suffixSubsessionSum[k] + suffixSessionSum[s + 1]);
						conditional[permId] = Math.exp(conditional[permId]);
					}	
					normalizedByL1(conditional);
					
					
					for (int permId = 0; permId < permLists.get(j).get(s).size(); ++permId) {
						int[] perm = permLists.get(j).get(s).get(permId);
						suffixSubsessionSum[suffixSubsessionSum.length - 1] = tempScores[curRow][perm[perm.length - 1]];
						for (int k = suffixSubsessionSum.length - 2; k >= 0; --k) 	suffixSubsessionSum[k] = suffixSubsessionSum[k + 1] + tempScores[curRow][perm[k]];

						double permPrefixSum = 0.0;
						for (int k = 0; k < perm.length; ++k) {
							permPrefixSum += (double) 1.0 / (suffixSubsessionSum[k] + suffixSessionSum[s + 1]);
//							tempScores[nextRow][perm[k]] += (1 - pZX)  * (double) permPrefixSum / cnt;
							tempScores[nextRow][perm[k]] += (1 - pZX)  * (double) permPrefixSum * conditional[permId];
						}
//						sessionSum += (double)  permPrefixSum / cnt;
						sessionSum += (double)  permPrefixSum * conditional[permId];
					}
					
					
					prefixSum += sessionSum;
				}
				
 			}

			
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
			
			
			
			/*
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
			
			
//			calcLogLikelihood();
			double L = calcLogLikelihood();
//			for (int i = 0; i < piProbs.length; ++i) L += Math.log(piProbs[i]);
			
			System.out.println("L = " + L);
			if (L < oldL || Math.abs((L - oldL) / oldL) < epsilon) {
				++stopCnt;
				oldL = L;
				if (stopCnt > 3) {
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
			double prob = 0.0;
			for (int id : semiRankedList.get(0)) prob += Math.log(scores[id]);
			for (int id : semiRankedList.get(1)) prob += Math.log(1 - scores[id]);
			prob += Math.log(correctProb);
			prob = (1 - correctProb) * piProbs[j] * pickingProb[semiRankedList.get(0).length] + Math.exp(prob);
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
		// ptau testing
		/*
		for (String rho : new String[]{"1", "1.5", "2", "2.5", "3"}) {
			SemiRankingDataSet trainDataSet = new SemiRankingDataSet();
			trainDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/exp/synthetic_params/rkdata_0.5_" + rho + ".txt");
			HashMap<Integer, Double> trainScores = trainDataSet.readGtScores("/Users/hzhuang/Work/beta/ranking/exp/synthetic_params/gtscore_0.5_" + rho + ".txt");
			System.out.println(trainDataSet.id2Name.size());
			
			
	//		SemiRankingDataSet testDataSet = new SemiRankingDataSet();
	//		testDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/synthetic_param/rankedlists_test.txt");
	//		System.out.println(testDataSet.id2Name.size());
	
			
			SimpleWorkerModelLearner learner = new SimpleWorkerModelLearner(trainDataSet);
			learner.trainWorkerModelParams(trainDataSet, trainScores);
			
			learner.saveModel("/Users/hzhuang/Work/beta/ranking/exp/synthetic_params/learnedmodel_0.5_" + rho + ".txt");
			
//			learner.trainRankings();
			
//			learner.evaluate("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label_test.txt");
//			learner.evaluateByROC("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label_test.txt");
		}
		/**/
		
		/*  //lambda testing
		for (int i = 1; i < 10; ++i) {
			SemiRankingDataSet trainDataSet = new SemiRankingDataSet();
			trainDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/exp/synthetic_params/rkdata_0." + i + "_2.txt");
			HashMap<Integer, Double> trainScores = trainDataSet.readGtScores("/Users/hzhuang/Work/beta/ranking/exp/synthetic_params/gtscore_0." + i + "_2.txt");
			System.out.println(trainDataSet.id2Name.size());
			
			
	//		SemiRankingDataSet testDataSet = new SemiRankingDataSet();
	//		testDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/synthetic_param/rankedlists_test.txt");
	//		System.out.println(testDataSet.id2Name.size());
	
			
			SimpleWorkerModelLearner learner = new SimpleWorkerModelLearner(trainDataSet);
			learner.trainWorkerModelParams(trainDataSet, trainScores);
			
			learner.saveModel("/Users/hzhuang/Work/beta/ranking/exp/synthetic_params/learnedmodel_0." + i + "_2.txt");
			
//			learner.trainRankings();
			
//			learner.evaluate("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label_test.txt");
//			learner.evaluateByROC("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label_test.txt");
		}
		/**/


		
		/*
		// SYNTHETIC data set
		SemiRankingDataSet trainDataSet = new SemiRankingDataSet();
		trainDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/synthetic/rankedlists.txt");
		HashMap<Integer, Double> trainScores = trainDataSet.readGtScores("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_score.txt");
		System.out.println(trainDataSet.id2Name.size());
		
		
		SemiRankingDataSet testDataSet = new SemiRankingDataSet();
		testDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/synthetic/rankedlists_test.txt");
		System.out.println(testDataSet.id2Name.size());

		
		SimpleWorkerModelLearner learner = new SimpleWorkerModelLearner(testDataSet);
		learner.trainWorkerModelParams(trainDataSet, trainScores);
		
		learner.trainRankings();
		
		learner.evaluate("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label_test.txt");
		learner.evaluateByROC("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label_test.txt");
		/**/

		// SYNTHETIC2 data set
		SemiRankingDataSet trainDataSet = new SemiRankingDataSet();
		trainDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/exp/0216_synthetic/temp/rkdata.txt");
		HashMap<Integer, Double> trainScores = trainDataSet.readGtScores("/Users/hzhuang/Work/beta/ranking/exp/0216_synthetic/temp/gtscore.txt");
		System.out.println(trainDataSet.id2Name.size());
		
		
		SemiRankingDataSet testDataSet = new SemiRankingDataSet();
		testDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/exp/0216_synthetic/temp/rkdata.test.txt");
		System.out.println(testDataSet.id2Name.size());

		
		SimpleWorkerModelLearner learner = new SimpleWorkerModelLearner(testDataSet, trainDataSet, trainScores);
//		learner.trainWorkerModelParams(trainDataSet, trainScores);
		
		learner.trainRankings();
		
		learner.evaluate("/Users/hzhuang/Work/beta/ranking/exp/0216_synthetic/temp/gtlabel.test.txt");
		learner.evaluateByROC("/Users/hzhuang/Work/beta/ranking/exp/0216_synthetic/temp/gtlabel.test.txt");
		/**/
		
		
		
		/*
		// LNKD data set
		SemiRankingDataSet trainDataSet = new SemiRankingDataSet();
		trainDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/job_509470.json.train.srk");
		HashMap<Integer, Double> trainScores = trainDataSet.readGtScores("/Users/hzhuang/Work/beta/ranking/data/fullres.out.train.filtered");
		System.out.println(trainDataSet.id2Name.size());
		
		
		SemiRankingDataSet testDataSet = new SemiRankingDataSet();
		testDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/job_509470.json.test.srk");
		System.out.println(testDataSet.id2Name.size());

		
		SimpleWorkerModelLearner learner = new SimpleWorkerModelLearner(testDataSet);
//		learner.logFileString = "/Users/hzhuang/Work/beta/ranking/exp/lnkd/loglikelihood.txt";
		learner.trainWorkerModelParams(trainDataSet, trainScores);
		
		learner.trainRankings();
		
		learner.evaluate("/Users/hzhuang/Work/beta/ranking/data/fullres.out.test.filtered");
		learner.evaluateByROC("/Users/hzhuang/Work/beta/ranking/data/fullres.out.test.filtered");
		/**/
		
		/*
		SemiRankingDataSet trainDataSet = new SemiRankingDataSet();
//		trainDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/job_509470.json.train.srk");
		trainDataSet.readSemiRankingLists("/Users/hzhuang/Work/beta/ranking/data/synthetic/rankedlists.txt");

		System.out.println(trainDataSet.id2Name.size());
//		HashMap<Integer, Integer> trainLabels = trainDataSet.readGtLabel("/Users/hzhuang/Work/beta/ranking/data/fullres.out.train.filtered");
//		HashMap<Integer, Integer> trainLabels = trainDataSet.readGtLabel("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label.txt");
		HashMap<Integer, Double> trainScores  = trainDataSet.readGtScores("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_score.txt");


		SimpleWorkerModelLearner learner = new SimpleWorkerModelLearner(trainDataSet);
		
		learner.logFileString = "/Users/hzhuang/Work/beta/ranking/exp/synthetic/loglikelihood.txt";
		learner.trainWorkerModelParams(trainDataSet, trainScores);
		
//		learner.trainWorkerModelParams(trainDataSet, trainLabels);
//		learner.trainRankings();
//		learner.evaluate("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label.txt");
//		learner.evaluateByROC("/Users/hzhuang/Work/beta/ranking/data/synthetic/ground_truth_label.txt");
		

		/**/
	}
	
	
}
