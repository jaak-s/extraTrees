package org.extratrees;

import static org.junit.Assert.*;

import org.junit.Test;

public class ExtraTreeTests {

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

	/** testing {@code et.getAllValues(input)} */
	@Test
	public void testGetAll() {
		ExtraTrees et = getSampleData(50, 5);
		et.learnTrees(5, 4, 10);
		// get all predictions by trees:
		Matrix all = et.getAllValues(et.input);
		double[] yhat = et.getValues(et.input);
		// check if their mean is equal to extraTree predictions:
		//System.out.println(all);
		for (int row=0; row<yhat.length; row++) {
			double sum = 0;
			for (int j=0; j<all.ncols; j++) {
				sum += all.get(row, j);
			}
			assertEquals("row="+row, yhat[row], sum/all.ncols, 1e-6);
		}
	}
	
	@Test
	public void testSelectTrees() {
		int ndata = 50;
		int ndim  = 5;
		ExtraTrees et = getSampleData(ndata, ndim);
		int ntrees = 10;
		et.learnTrees(5, 4, ntrees);
		boolean[] selection = new boolean[ntrees];
		int[] treeIds = {0, 5};
		for (int i=0; i<selection.length; i++) {
			selection[ i ] = false;
		}
		for (int i=0; i<treeIds.length; i++) {
			selection[ treeIds[i] ] = true;
		}
		// making sub-extraTree based on selection:
		ExtraTrees et2 = et.selectTrees(selection);
		
		Matrix m  = et.getAllValues(et.input);
		Matrix m2 = et2.getAllValues(et2.input);
		
		// checking that the chosen trees are the same:
		assertEquals( treeIds.length, et2.trees.size() );
		assertEquals( treeIds.length, m2.ncols );
		for (int i=0; i<treeIds.length; i++) {
			assertTrue( et2.trees.get(i)==et.trees.get(treeIds[i]) );
			// checking the getAll: (first)
			assertEquals( m2.get(0, i), m.get(0, treeIds[i]), 1e-6);
		}
		
	}

}
