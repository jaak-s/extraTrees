package org.extratrees;

import java.util.ArrayList;
import java.util.Random;

import org.extratrees.data.Matrix;
import org.junit.Test;

public class BenchmarkRange {

	@Test
	public void testRange() {
		int N = 10*1000*1000;
		int col = 1;
		Matrix m = new Matrix(N, 2);
		for (int i=0; i<N; i++) {
			m.set(i, col, Math.random());
		}
		int[] ids = AbstractTrees.seq(N);
		double[] range;
		Timer.tic();
		range = AbstractTrees.getRange(ids, col, m);
		Timer.toc("AbstractTrees.getRange(10M sequential)");
	}
	
	public static int[] randIds(int size, int max) {
		ArrayList<Integer> allIds = AbstractTrees.arrayToList( AbstractTrees.seq( max ) );
		Timer.tic();
		ShuffledIterator<Integer> shuffle = new ShuffledIterator<Integer>(allIds, new Random());
		
		int[] subset = new int[ size ];
		for (int i=0; i < subset.length; i++) {
			subset[i] = shuffle.next();
		}
		Timer.toc("shuffle");
		return subset;
	}
	
	@Test
	public void testRangeOnRandom() {
		int N = 1*1000*1000;
		int col = 1;
		Matrix m = new Matrix(N, 2);
		for (int i=0; i<N; i++) {
			m.set(i, col, Math.random());
		}
		
		int[] ids = randIds(N/2, N);
		double[] range;
		Timer.tic();
		range = AbstractTrees.getRange(ids, col, m);
		Timer.toc("AbstractTrees.getRange(500k random-ids)");
	}

}
