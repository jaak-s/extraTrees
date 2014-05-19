package org.extratrees;

import org.extratrees.data.Matrix;
import org.junit.Test;

public class BenchmarkRange {

	@Test
	public void testRange() {
		int N = 10*1000*1000;
		int col = 2;
		Matrix m = new Matrix(N, 5);
		for (int i=0; i<N; i++) {
			m.set(i, col, Math.random());
		}
		int[] ids = AbstractTrees.seq(N);
		double[] range;
		Timer.tic();
		range = AbstractTrees.getRange(ids, col, m);
		Timer.toc("AbstractTrees.getRange");
	}

}
