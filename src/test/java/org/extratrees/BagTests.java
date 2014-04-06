package org.extratrees;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class BagTests {
	public static ExtraTrees getSampleData(int ndata, int ndim) {
		double[] output = new double[ndata];
		double[] v = new double[ndata*ndim];
		for (int i=0; i<v.length; i++) {
			v[i] = Math.random();
		}
		Matrix m = new Matrix(v, ndata, ndim);
		// generate values for all outputs
		for (int row=0; row<output.length; row++) {
			m.set(row, 2, 0.5);
			output[row] = m.get(row, 1)+0.2*m.get(row, 3);
		}
		ExtraTrees et = new ExtraTrees(m, output);
		return et;
	}

	
	@Test
	public void test() {
		ExtraTrees et = getSampleData(100, 10);
		assertTrue( et.bagElems == null );
		assertTrue( et.bagSizes == null );
		
		et.setBagging(90);
		assertArrayEquals(new int[]{90}, et.bagSizes);
		
		int[] ids = et.getInitialSamples();
		assertTrue("Random bag should contain 90 samples.", ids.length == 90);
		
		int[] bagLabels = new int[100];
		for (int i = 0; i<bagLabels.length; i++) {
			bagLabels[i] = (i < 30) ?0 :1;
		}
		int[] bagSizes = new int[]{9, 50};
		et.setBagging( bagSizes, bagLabels );
		
		assertArrayEquals(bagSizes, et.bagSizes);
		assertArrayEquals(AbstractTrees.seq(30), et.bagElems[0]);
		assertArrayEquals(AbstractTrees.seq(30, 100), et.bagElems[1]);

		ids = et.getInitialSamples();
		assertTrue("Random bag should contain 59 samples.", ids.length == AbstractTrees.sum(bagSizes) );
		//System.out.println( AbstractTrees.arrayToList(ids));

		et.learnTrees(5, 3, 10);
	}

}
