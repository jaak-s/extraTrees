package org.extratrees;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;

import java.util.Random;

import org.extratrees.data.Matrix;
import org.junit.Test;

public class SubsetTests {
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
		assertTrue( et.subsetElems == null );
		assertTrue( et.subsetSizes == null );
		
		et.setSubsetting(90);
		assertArrayEquals(new int[]{90}, et.subsetSizes);
		
		int[] ids = et.getInitialSamples(new Random());
		assertTrue("Random subset should contain 90 samples.", ids.length == 90);
		
		int[] subsetGroups = new int[100];
		for (int i = 0; i<subsetGroups.length; i++) {
			subsetGroups[i] = (i < 30) ?0 :1;
		}
		int[] subsetSizes = new int[]{9, 50};
		et.setSubsetting( subsetSizes, subsetGroups );
		
		assertArrayEquals(subsetSizes, et.subsetSizes);
		assertArrayEquals(AbstractTrees.seq(30), et.subsetElems[0]);
		assertArrayEquals(AbstractTrees.seq(30, 100), et.subsetElems[1]);

		ids = et.getInitialSamples(new Random());
		assertTrue("Random subset should contain 59 samples.", ids.length == AbstractTrees.sum(subsetSizes) );

		et.learnTrees(5, 3, 10);
	}

}
