package org.extratrees;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.extratrees.AbstractTrees.CutResult;
import org.extratrees.data.Matrix;
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
		Matrix m2 = et2.getAllValues(et.input);
		
		// checking that the chosen trees are the same:
		assertEquals( treeIds.length, et2.trees.size() );
		assertEquals( treeIds.length, m2.ncols );
		for (int i=0; i<treeIds.length; i++) {
			assertTrue( et2.trees.get(i)==et.trees.get(treeIds[i]) );
			// checking the getAll: (first)
			assertEquals( m2.get(0, i), m.get(0, treeIds[i]), 1e-6);
		}
		
	}
	
	public static Matrix getSampleData(int ndata) {
		int ndim = 3;
		Matrix m = new Matrix(ndata, ndim);
		for (int i=0; i<ndata; i++) {
			m.set(i, 0, i);
			m.set(i, 1, i<ndata/2 ?1.0 :2.0);
			m.set(i, 2, -1);
		}
		// generate values for all outputs
		return m;
	}
	
	@Test
	public void testSplitIds() {
		Matrix x = getSampleData(20);
		int[] ids = new int[10];
		for (int i=0; i<ids.length; i++) {
			ids[i] = i*2;
		}
		int[][] split = AbstractTrees.splitIds(x, ids, 0, 9.9);
		assertArrayEquals(new int[]{0,2,4,6,8}, split[0] );
		assertArrayEquals(new int[]{10,12,14,16,18}, split[1] );
	}
	
	@Test
	public void testConstantY() {
		Matrix input = getSampleData(20);
		final double fixedOutput = 1.0;
		double[] output = new double[input.nrows];
		for (int i=0; i<output.length; i++) {
			output[i] = fixedOutput;
		}
		int ntrees = 10;
		ExtraTrees et = new ExtraTrees(input, output);
		et.learnTrees(5, 4, ntrees);
		double[] outputHat = et.getValues(input);
		for (int i=0; i<output.length; i++) {
			assertEquals(fixedOutput, outputHat[i], 1e-10);
		}
	}

	// data for testing weights:
	public static ExtraTrees getSampleDataW(boolean weights) {
		int ndata = 10;
		int ndim = 2;
		double[] output = new double[ndata];
		double[] v = new double[ndata*ndim];
		double[] w = new double[ndata];
		Matrix m = new Matrix(v, ndata, ndim);
		for (int i=0; i < m.nrows; i++) {
			m.set(i, 0, i);
			m.set(i, 1, (i < m.nrows/2) ?1.0 :0.0 );
			if (i < 5) {
				output[i] = 1.0;
			} else {
				output[i] = 0.0;
			}
			w[i] = (i < m.nrows/2) ?1.0 :0.2;
		}
		ExtraTrees et = new ExtraTrees(m, output);
		if (weights) {
			et.setWeights(w);
		}
		return et;
	}


	@Test
	public void testWeights() {
		ExtraTrees et = getSampleDataW(false);
		ExtraTrees etw = getSampleDataW(true);
		
		assertTrue( ! et.useWeights );
		assertTrue( etw.useWeights );
		
		int[] all = AbstractTrees.seq(etw.input.nrows());
		CutResult resultw = new CutResult();
		CutResult result  = new CutResult();
		etw.calculateCutScore(all, 0, 5.5, resultw);
		et.calculateCutScore(all, 0, 5.5, result);
		// test sanity checks:
		assertEquals(6, resultw.countLeft);
		assertEquals(4, resultw.countRight);
		assertTrue(resultw.rightConst);
		assertTrue( ! resultw.leftConst);
		
		double mean   = 5.0 / 6.0;
		double meanw  = 5.0 / 5.2;
		// score for case of no weights:
		double score  = (Math.pow(1.0-mean,  2)*5 + mean*mean);
		// score for case of weights:
		double scorew = (Math.pow(1.0-meanw, 2)*5 + meanw*meanw*0.2);
		assertEquals(score,  result.score, 1e-7);
		assertEquals(scorew, resultw.score, 1e-7);
	}

	@Test
	public void testWeightTraining() {
		ExtraTrees etw = getSampleDataW(true);
		etw.learnTrees(4, 1, 20);
		etw.getValues(etw.input);
	}

}
