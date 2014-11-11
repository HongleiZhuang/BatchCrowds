package learner;

public class SemiRankingLearnerTester {

	String rkdataFileName, resultsFileName, scoresFileName, groundTruthFileName;
	int method;
	static final public int NAIVE = 0, PL = 1;
	static final public int PRF = 1, AUC = 2;
	
	private void loadParameter(String[] args) {
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-data")) {
				this.rkdataFileName = args[++i];
			} else if (args[i].equals("-result")) {
				this.resultsFileName = args[++i];
			} else if (args[i].equals("-gt")) {
				this.groundTruthFileName = args[++i];
			} else if (args[i].equals("-scrout")) {
				this.scoresFileName = args[++i];
			} else if (args[i].equals("-method")) {
				++i;
				if (args[i].equalsIgnoreCase("naive")) method = NAIVE;
				else if (args[i].equalsIgnoreCase("pl")) method = PL;
			}
		}
	}
	
	
	
	
	
}
