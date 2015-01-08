package learner;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;

import model.SemiRankingDataSet;

public class SemiRankingLearnerTester {
	String trainRkdataFileName, trainGtFileName;
	String testGtFileName;
	String rkdataFileName, resultsFileName, scoresFileName;
	int method = -1, eval = 0;
	static final public int NAIVE = 0, PL = 1, WORKER = 2, NAIVE_TRAINED = 3;
	static final public int PRF = 1, AUC = 2;
	
	public SemiRankingLearnerTester(String[] args) {
		loadParameter(args);
	}
	
	private void loadParameter(String[] args) {
		for (int i = 0; i < args.length; i++) {
			if (args[i].equalsIgnoreCase("-data")) {
				this.rkdataFileName = args[++i];
			} else if (args[i].equalsIgnoreCase("-result")) {
				this.resultsFileName = args[++i];
			} else if (args[i].equalsIgnoreCase("-gt")) {
				this.testGtFileName = args[++i];
			} else if (args[i].equalsIgnoreCase("-scrout")) {
				this.scoresFileName = args[++i];
			} else if (args[i].equalsIgnoreCase("-traindata")) {
				this.trainRkdataFileName = args[++i];
			} else if (args[i].equalsIgnoreCase("-traingt")) {
				this.trainGtFileName = args[++i];
			} else if (args[i].equalsIgnoreCase("-method")) {
				++i;
				if (args[i].equalsIgnoreCase("naive")) method = NAIVE;
				else if (args[i].equalsIgnoreCase("naive_trained")) method = NAIVE_TRAINED;
				else if (args[i].equalsIgnoreCase("pl")) method = PL;
				else if (args[i].equalsIgnoreCase("worker")) method = WORKER;
			} else if (args[i].equalsIgnoreCase("-eval")) {
				++i;
				if (args[i].equalsIgnoreCase("prf")) eval = PRF;
				else if (args[i].equalsIgnoreCase("auc")) eval = AUC;
			}
		}
	}
	
	public void generateResults() throws Exception {
		ScoreBasedSemiRankingLearner learner = null;
		SemiRankingDataSet rkdata = new SemiRankingDataSet();
		rkdata.readSemiRankingLists(this.rkdataFileName);

		if (method == NAIVE) {
			learner = NaiveLearner.createNaiveLearner(rkdata);
		}
		else if (method == PL) {
			learner = PlackettLuceSemiMMLearner.createPlacketBasedSemiRankingLearner(rkdata);
		}
		else if (method == WORKER) {
			SemiRankingDataSet trainDataSet = new SemiRankingDataSet();
			trainDataSet.readSemiRankingLists(this.trainRkdataFileName);
			HashMap<Integer, Double> trainScores = trainDataSet.readGtScores(this.trainGtFileName);
			learner = SimpleWorkerModelLearner.createSimpleWorkerModelLearner(rkdata, trainDataSet, trainScores);
		}
		else if (method == NAIVE_TRAINED) {
			SemiRankingDataSet trainDataSet = new SemiRankingDataSet();
			trainDataSet.readSemiRankingLists(this.trainRkdataFileName);
			HashMap<Integer, Double> trainScores = trainDataSet.readGtScores(this.trainGtFileName);
			learner = NaiveLearner.createNaiveLearnerWithThresholdTraining(rkdata, trainDataSet, trainScores);
		}
		else {
			System.out.println("Method specified error!");
			return;
		}
		
		learner.trainRankings();
		
		if (eval > 0) {
			double[] result = null;
			if (eval == PRF) {
				result = learner.evaluate(this.testGtFileName);
			}
			else if (eval == AUC) {
				result = learner.evaluateByROC(this.testGtFileName);				
			}
			BufferedWriter resBw = new BufferedWriter(new FileWriter(new File(this.resultsFileName)));
			if (result == null) {
				System.out.println("No results to output!");
				resBw.flush();
				resBw.close();
				return;
			}
			else {
				for (double val : result) resBw.write(val + "\t");
			}
			resBw.flush();
			resBw.close();
		}
	}
	
	static public void main(String[] args) throws Exception {
		SemiRankingLearnerTester tester = new SemiRankingLearnerTester(args);
		tester.generateResults();
	}
}
