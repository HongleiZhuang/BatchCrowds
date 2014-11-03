package learner;

import model.SemiRankingDataSet;



public abstract class ScoreBasedSemiRankingLearner {
	SemiRankingDataSet rkdata;
	
	double[] scores;
	
	abstract void trainRankings();
}
