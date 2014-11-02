package learner;

import model.RankingDataSet;

public abstract class ScoreBasedRankingLearner {
	RankingDataSet rkdata;
	
	double[] scores;
	
	abstract void trainRankings();
}
