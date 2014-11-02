package util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;



public class MapSorter {
	
	static public <K, V extends Comparable<? super V>> ArrayList<Entry<K, V>> sortByValue(Map<K, V> mapToSort) {
		return MapSorter.sortByValue(mapToSort, true);
	}
	
	static public <K, V extends Comparable<? super V>> ArrayList<Entry<K, V>> sortByValue(Map<K, V> mapToSort, boolean naturalOrder) {
		ArrayList<Entry<K, V>> entryList = new ArrayList<Map.Entry<K,V>>(mapToSort.entrySet());
		if (naturalOrder) 
	        Collections.sort(entryList, new Comparator<Map.Entry<K, V>>() {
	            public int compare( Map.Entry<K, V> e1, Map.Entry<K, V> e2 ){
	                return (e1.getValue()).compareTo(e2.getValue());
	            }
	        } );
		else 
	        Collections.sort(entryList, new Comparator<Map.Entry<K, V>>() {
	            public int compare( Map.Entry<K, V> e1, Map.Entry<K, V> e2 ){
	                return -(e1.getValue()).compareTo(e2.getValue());
	            }
	        } );			
		return entryList;
	}
	
	static public <K extends Comparable<? super K>, V> ArrayList<Entry<K, V>> sortByKey(Map<K, V> mapToSort) {
		return MapSorter.sortByKey(mapToSort, true);
	}
	
	static public <K extends Comparable<? super K>, V> ArrayList<Entry<K, V>> sortByKey(Map<K, V> mapToSort, boolean naturalOrder) {
		ArrayList<Entry<K, V>> entryList = new ArrayList<Map.Entry<K,V>>(mapToSort.entrySet());
		if (naturalOrder)
	        Collections.sort(entryList, new Comparator<Map.Entry<K, V>>() {
	            public int compare( Map.Entry<K, V> e1, Map.Entry<K, V> e2 ){
	                return (e1.getKey()).compareTo(e2.getKey());
	            }
	        } );
		else
	        Collections.sort(entryList, new Comparator<Map.Entry<K, V>>() {
	            public int compare( Map.Entry<K, V> e1, Map.Entry<K, V> e2 ){
	                return -(e1.getKey()).compareTo(e2.getKey());
	            }
	        } );
		return entryList;
	}
	
	
	static public void main(String[] args) {
		HashMap<Integer, Double> testMap = new HashMap<Integer, Double>();
		for (int i = 0; i < 10; ++i) {
			testMap.put(i, Math.random());
			System.out.println(i + ": " +  testMap.get(i));
		}
		System.out.println("================");
		ArrayList<Entry<Integer, Double>> sortedList = MapSorter.sortByValue(testMap, false);
		for (int i = 0; i < sortedList.size(); ++i) {
			Entry<Integer, Double> e = sortedList.get(i);	
			System.out.println(e.getKey() + ": " + e.getValue());
		}
	}
}
