package org.extratrees;
import static org.junit.Assert.*;

import java.util.Random;

import org.extratrees.FactorExtraTrees;
import org.extratrees.Matrix;
import org.junit.Test;


public class FactorTests {

	@Test
	public void testGini() {
		int[] counts;
		counts = new int[]{0,0,5,0,0};
		assertEquals(0, FactorExtraTrees.getGiniIndex(counts), 1e-6 );
		
		counts = new int[]{0,0,1,0,0};
		assertEquals(0, FactorExtraTrees.getGiniIndex(counts), 1e-6 );
		
		counts = new int[]{1,0,1,0,0};
		assertEquals(0.5, FactorExtraTrees.getGiniIndex(counts), 1e-6 );
		
		counts = new int[]{7,0,7,0,0};
		assertEquals(0.5, FactorExtraTrees.getGiniIndex(counts), 1e-6 );
		
		counts = new int[]{7,7,7,7,7};
		assertEquals(0.8, FactorExtraTrees.getGiniIndex(counts), 1e-6 );
		
		counts = new int[]{4,0,0,0,1};
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
		FactorExtraTrees et = getSampleData(ndata, 5);
		et.learnTrees(5, 4, 10);
		// get all predictions by trees:
		Matrix all = et.getAllValues(et.input);
		//int[] yhat = et.getValues(et.input);
		// check if their mean is equal to extraTree predictions:
		System.out.println(all);
		/*
		for (int row=0; row<yhat.length; row++) {
			double sum = 0;
			for (int j=0; j<all.ncols; j++) {
				sum += all.get(row, j);
			}
			assertEquals("row="+row, yhat[row], sum/all.ncols, 1e-6);
		}*/
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
		Matrix m2 = et2.getAllValues(et2.input);
		
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

	/*
	@Test
	public void testSSL() {
		int ndata = 50;
		int nUnlabeled = 500;
		FactorExtraTrees et = getSampleData2(ndata);
		et.setUnlabeled( getSampleData2(nUnlabeled).input );
		int ntrees = 10;
		et.learnTrees(5, 4, ntrees);
	}*/

}
