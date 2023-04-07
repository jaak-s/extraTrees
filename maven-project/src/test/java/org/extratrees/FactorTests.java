package org.extratrees;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Random;

import org.extratrees.AbstractTrees.CutResult;
import org.extratrees.data.Matrix;
import org.junit.Test;


public class FactorTests {

	@Test
	public void testGini() {
		double[] counts;
		counts = new double[]{0,0,5,0,0};
		assertEquals(0, FactorExtraTrees.getGiniIndex(counts), 1e-6 );
		
		counts = new double[]{0,0,1,0,0};
		assertEquals(0, FactorExtraTrees.getGiniIndex(counts), 1e-6 );
		
		counts = new double[]{1,0,1,0,0};
		assertEquals(0.5, FactorExtraTrees.getGiniIndex(counts), 1e-6 );
		
		counts = new double[]{1,0,1,0,1};
		assertEquals(2.0/3.0, FactorExtraTrees.getGiniIndex(counts), 1e-6 );
		
		counts = new double[]{7,0,7,0,0};
		assertEquals(0.5, FactorExtraTrees.getGiniIndex(counts), 1e-6 );
		
		counts = new double[]{7,7,7,7,7};
		assertEquals(0.8, FactorExtraTrees.getGiniIndex(counts), 1e-6 );
		
		counts = new double[]{4,0,0,0,1};
		assertEquals(0.32, FactorExtraTrees.getGiniIndex(counts), 1e-6 );
	}

	/**
	 * @param ndata
	 * @param ndim
	 * @return FactorExtraTrees, with factors of 0, 1, 2.
	 */
	public static FactorExtraTrees getSampleData(int ndata, int ndim) {
		int[] output = new int[ndata];
		double[] v   = new double[ndata*ndim];
		for (int i=0; i<v.length; i++) {
			v[i] = Math.random();
		}
		Matrix m = new Matrix(v, ndata, ndim);
		// generate values for all outputs
		for (int row=0; row<output.length; row++) {
			m.set(row, 2, 0.5);
			output[row] = (int) Math.floor(m.get(row, 1)+2*m.get(row, 3));
		}
		FactorExtraTrees et = new FactorExtraTrees(m, output);
		return et;
	}

	@Test
	public void testFactorExtraTrees() {
		FactorExtraTrees fet = getSampleData(100, 5);
		assertEquals("number of factors",  3, fet.getnFactors()); 
		fet.learnTrees(2, 4, 10);
	}
	
	@Test
	public void testBug() {
		int ndata=20;
		int[] output = new int[ndata];
		Matrix m = new Matrix(ndata, 1);
		for (int i=0; i<ndata; i++) {
			m.set(i, 0, i/(double)ndata);
			output[i] = i*2<ndata ?0 :1;
		}
		FactorExtraTrees fet = new FactorExtraTrees(m, output);
		fet.setNumRandomCuts(4);
		fet.learnTrees(1, 1, 1);
	}
	
	/** testing {@code et.getAllValues(input)} */
	@Test
	public void testGetAll() {
		int ndata = 40;
		int nTrees = 10;
		FactorExtraTrees et = getSampleData(ndata, 5);
		et.learnTrees(5, 4, nTrees);
		// get all predictions by trees:
		Matrix all = et.getAllValues(et.input);
		assertEquals(et.input.nrows(), all.nrows());
		assertEquals(nTrees, all.ncols);
		int[] yhat = et.getValues(et.input);
		// check if their mean is equal to extraTree predictions:
		//System.out.println(all);
		int errors = 0;
		for (int row=0; row<yhat.length; row++) {
			if (yhat[row] != et.output[row]) {
				errors++;
			}
		}
		System.out.println( String.format("Error rate: %1.3f", errors / (double) yhat.length) );
	}

	@Test
	public void testSelectTrees() {
		int ndata = 50;
		int ndim  = 5;
		FactorExtraTrees et = getSampleData(ndata, ndim);
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
		FactorExtraTrees et2 = et.selectTrees(selection);
		
		Matrix m  = et.getAllValues(et.input);
		Matrix m2 = et2.getAllValues(et.input);
		
		// checking that the chosen trees are the same:
		assertEquals( treeIds.length, et2.trees.size() );
		assertEquals( treeIds.length, m2.ncols );
		for (int i=0; i<treeIds.length; i++) {
			assertTrue( et2.trees.get(i)==et.trees.get(treeIds[i]) );
			// checking the getAll: (first)
			for (int row=0; row<m.nrows; row++) {
				assertEquals( m2.get(row, i), m.get(row, treeIds[i]), 1e-6);
			}
		}
	}
	
	public FactorExtraTrees getSampleData2(int ndata) {
		int[] output = new int[ndata];
		int ndim = 3;
		Random random = new Random();
		Matrix m = new Matrix(ndata, ndim);
		for (int i=0; i<m.nrows; i++) {
			output[i] = Math.random()>0.5 ?1 :0;
			if (output[i]==1) {
				m.set(i, 0, random.nextGaussian() + 1 );
				m.set(i, 1, random.nextGaussian() + 3 );
				m.set(i, 2, random.nextGaussian() + 5 );
			} else {
				m.set(i, 0, random.nextGaussian() );
				m.set(i, 1, random.nextGaussian() );
				m.set(i, 2, random.nextGaussian() );
			}
		}
		FactorExtraTrees et = new FactorExtraTrees(m, output);
		return et;
	}

	// data for testing weights:
	public static FactorExtraTrees getSampleDataW(boolean weights) {
		int ndata = 10;
		int ndim = 2;
		int[] output = new int[ndata];
		double[] v = new double[ndata*ndim];
		double[] w = new double[ndata];
		Matrix m = new Matrix(v, ndata, ndim);
		for (int i=0; i < m.nrows; i++) {
			m.set(i, 0, i);
			m.set(i, 1, (i < m.nrows/2) ?1.0 :0.0 );
			if (i < 5) {
				output[i] = 1;
			} else {
				output[i] = 0;
			}
			w[i] = (i < m.nrows/2) ?1.0 :0.2;
		}
		FactorExtraTrees et = new FactorExtraTrees(m, output);
		if (weights) {
			et.setWeights(w);
		}
		return et;
	}

	@Test
	public void testWeights() {
		FactorExtraTrees et = getSampleDataW(false);
		FactorExtraTrees etw = getSampleDataW(true);
		
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
		assertTrue(result.rightConst);
		assertTrue( ! result.leftConst);
		
		double giniLeft  = 1 - Math.pow(5.0/6.0, 2) - Math.pow(1.0/6, 2);
		double giniLeftw = 1 - Math.pow(5.0/5.2, 2) - Math.pow(0.2/5.2, 2);
		// score for case of no weights:
		double score  = giniLeft * 6.0 / 10;
		// score for case of weights:
		double scorew = giniLeftw * 5.2 / 6.0;
		assertEquals(score,  result.score, 1e-7);
		assertEquals(scorew, resultw.score, 1e-7);
	}

	@Test
	public void testWeightTraining() {
		FactorExtraTrees etw = getSampleDataW(true);
		etw.learnTrees(2, 1, 20);
		etw.getValues(etw.input);
	}

}
