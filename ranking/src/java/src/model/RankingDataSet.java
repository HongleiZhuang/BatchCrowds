package model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

public class RankingDataSet {
	public ArrayList<ArrayList<Integer>> rankingLists ;
	public HashMap<String, Integer> name2Id; 
	public ArrayList<String> id2Name;
	
		
	public void addRanking(ArrayList<Integer> rankedList) {
		rankingLists.add(rankedList);
	}
	
	
	public void readRankingLists(String fileName) throws Exception  {
		name2Id = new HashMap<String, Integer>();
		id2Name = new ArrayList<String>();
		rankingLists = new ArrayList<ArrayList<Integer>>();
		BufferedReader br = new BufferedReader (new FileReader( new File(fileName)));
		String s;
		while ((s = br.readLine()) != null) {
			s = s.trim();
			ArrayList<Integer> rankedList = new ArrayList<Integer>();
			String[] slist = s.split("\\s+");
			for (String e : slist) {
				if (!name2Id.containsKey(e)) {
					name2Id.put(e, id2Name.size());
					id2Name.add(e);
				}
				int id = name2Id.get(e);
				rankedList.add(id);
			}
			if (rankedList.size() > 0) {
				addRanking(rankedList);
			}
		}
		br.close();
	}
	
}
