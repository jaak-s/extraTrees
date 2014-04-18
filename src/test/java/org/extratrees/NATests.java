package org.extratrees;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class NATests {
	/**
	 * @param ndata
	 * @param ndim
	 */
	public static FactorExtraTrees getFET(int ndata, int ndim) {
		int[] output = new int[ndata];
		Matrix m = new Matrix(ndata, ndim);
		// generate values for all outputs
		for (int row=0; row<output.length; row++) {
			m.set(row, 1, row/(double)output.length);
			m.set(row, 2, 0.5);
			output[row] = m.get(row, 1)>0.55 ?1 :0;
		}
		FactorExtraTrees et = new FactorExtraTrees(m, output);
		return et;
	}

	/**
	 * @param ndata
	 * @param ndim
	 */
	public static ExtraTrees getET(int ndata, int ndim) {
		double[] output = new double[ndata];
		Matrix m = new Matrix(ndata, ndim);
		// generate values for all outputs
		for (int row=0; row<output.length; row++) {
			m.set(row, 1, row/(double)output.length);
			m.set(row, 2, 0.5);
			output[row] = m.get(row, 1)>0.55 ?1 :0;
		}
		ExtraTrees et = new ExtraTrees(m, output);
		return et;
	}


	@Test
	public void testFactorNACalc() {
		FactorExtraTrees et = getFET(10, 5);
		double gini;
		gini = 1 - (0.4*0.4 + 0.6*0.6);
		assertEquals( gini, et.get1NaNScore( AbstractTrees.seq(10) ), 1e-6);
		
		gini = 1 - (Math.pow(3.0/9.0, 2) + Math.pow(6.0/9.0, 2));
		assertEquals( gini, et.get1NaNScore( AbstractTrees.seq(9) ), 1e-6);
		
	}

	@Test
	public void testRegressionNACalc() {
		ExtraTrees et = getET(10, 5);
		double var, mean;
		mean = 4.0 / 10.0;
		var  = 0.6 * Math.pow(mean, 2) + 0.4 * Math.pow(1-mean, 2);
		assertEquals( var, et.get1NaNScore( AbstractTrees.seq(10) ), 1e-6);
		
		mean = 3.0 / 9.0;
		var  = 6/9.0 * Math.pow(mean, 2) + 3/9.0 * Math.pow(1-mean, 2);
		assertEquals( var, et.get1NaNScore( AbstractTrees.seq(9) ), 1e-6);
	}

}
