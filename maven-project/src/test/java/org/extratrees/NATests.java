package org.extratrees;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.extratrees.AbstractTrees.CutResult;
import org.extratrees.data.Matrix;
import org.junit.Test;

public class NATests {
	/**
	 * @param ndata
	 * @param ndim
	 */
	public static FactorExtraTrees getFET(int ndata, int ndim, boolean useWeights) {
		int[] output = new int[ndata];
		Matrix m = new Matrix(ndata, ndim);
		// generate values for all outputs
		for (int row=0; row<output.length; row++) {
			m.set(row, 1, row/(double)output.length);
			m.set(row, 2, 0.5);
			if (row == 5 || row == 6 || row == 7) {
				m.set(row, 2, Double.NaN);
			}
			output[row] = m.get(row, 1)>0.55 ?1 :0;
		}
		FactorExtraTrees et = new FactorExtraTrees(m, output);
		et.setHasNaN(true);
		if (useWeights) {
			double[] w = new double[ndata];
			for (int i=0; i < w.length; i++) {
				w[i] = 0.5;
			}
			et.setWeights(w);
		}
		return et;
	}

	@Test
	public void testFactorNACalc() {
		FactorExtraTrees et = getFET(10, 5, false);
		FactorExtraTrees etw = getFET(10, 5, true);
		double gini;
		gini = 1 - (0.4*0.4 + 0.6*0.6);
		assertEquals( gini, et.get1NaNScore( AbstractTrees.seq(10) ), 1e-6);
		
		gini = 1 - (Math.pow(3.0/9.0, 2) + Math.pow(6.0/9.0, 2));
		assertEquals( gini, et.get1NaNScore( AbstractTrees.seq(9) ), 1e-6);

		// testing NaN counts:
		CutResult cr;
		cr = new CutResult();
		et.calculateCutScore(AbstractTrees.seq(9), 2, 0.5, cr);
		assertEquals( 3.0, cr.nanWeigth, 1e-6 );
		
		cr = new CutResult();
		et.calculateCutScore(AbstractTrees.seq(9), 1, 0.5, cr);
		assertEquals( 0.0, cr.nanWeigth, 1e-6 );
		
		// testing weights:
		cr = new CutResult();
		etw.calculateCutScore(AbstractTrees.seq(9), 2, 0.5, cr);
		assertEquals( 1.5, cr.nanWeigth, 1e-6 );
	}

	@Test
	public void testFactorNALearn() {
		int ndim = 5;
		FactorExtraTrees et = getFET(100, ndim, false);
		FactorExtraTrees etw = getFET(100, ndim, true);
		et.learnTrees(3, 3, 5);
		etw.learnTrees(3, 3, 5);
		
		double[] x = new double[ndim];
		for (int i=0; i < x.length; i++) { 
			x[i] = Double.NaN;
		}
		int[] val;
		val = et.getValues( new Matrix(x, 1, ndim) );
		assertEquals( -1, val[0] );
		val = etw.getValues( new Matrix(x, 1, ndim) );
		assertEquals( -1, val[0] );
		
		val = et.getValues( et.input );
		for (int i = 0; i < val.length; i++) {
			assertTrue( val[i] >= 0 );
		}
	}


	/**
	 * @param ndata
	 * @param ndim
	 */
	public static ExtraTrees getET(int ndata, int ndim, boolean useWeights) {
		double[] output = new double[ndata];
		Matrix m = new Matrix(ndata, ndim);
		// generate values for all outputs
		for (int row=0; row < output.length; row++) {
			m.set(row, 1, row/(double)output.length);
			m.set(row, 2, 0.5);
			if (row == 5 || row == 6 || row == 7) {
				m.set(row, 2, Double.NaN);
			}
			output[row] = m.get(row, 1)>0.55 ?1 :0;
		}
		ExtraTrees et = new ExtraTrees(m, output);
		et.setHasNaN(true);
		if (useWeights) {
			double[] w = new double[ndata];
			for (int i=0; i < w.length; i++) {
				w[i] = 0.5;
			}
			et.setWeights(w);
		}
		return et;
	}

	@Test
	public void testRegressionNACalc() {
		ExtraTrees et = getET(10, 5, false);
		ExtraTrees etw = getET(10, 5, true);
		double var, mean;
		mean = 4.0 / 10.0;
		var  = 0.6 * Math.pow(mean, 2) + 0.4 * Math.pow(1-mean, 2);
		assertEquals( var, et.get1NaNScore( AbstractTrees.seq(10) ), 1e-6);
		
		mean = 3.0 / 9.0;
		var  = 6/9.0 * Math.pow(mean, 2) + 3/9.0 * Math.pow(1-mean, 2);
		assertEquals( var, et.get1NaNScore( AbstractTrees.seq(9) ), 1e-6);
		
		// testing NaN counts:
		CutResult cr;
		cr = new CutResult();
		et.calculateCutScore(AbstractTrees.seq(9), 2, 0.5, cr);
		assertEquals( 3.0, cr.nanWeigth, 1e-6 );
		
		cr = new CutResult();
		et.calculateCutScore(AbstractTrees.seq(9), 1, 0.5, cr);
		assertEquals( 0.0, cr.nanWeigth, 1e-6 );
		
		// testing weights:
		cr = new CutResult();
		etw.calculateCutScore(AbstractTrees.seq(9), 2, 0.5, cr);
		assertEquals( 1.5, cr.nanWeigth, 1e-6 );
	}

	@Test
	public void testRegressionNALearn() {
		int ndim = 5;
		ExtraTrees et = getET(100,  ndim, false);
		ExtraTrees etw = getET(100, ndim, true);
		et.learnTrees(3, 3, 5);
		etw.learnTrees(3, 3, 5);
		
		double[] x = new double[ndim];
		for (int i=0; i < x.length; i++) { 
			x[i] = Double.NaN;
		}
		double[] val;
		val = et.getValues( new Matrix(x, 1, ndim) );
		assertTrue( Double.isNaN(val[0]) );
		val = etw.getValues( new Matrix(x, 1, ndim) );
		assertTrue( Double.isNaN(val[0]) );
		
		// checking if getRange works with NaN
		double[] col2   = ((Matrix)et.input).getCol(2);
		double[] range2 = AbstractTrees.getRange( col2 );
		assertEquals(0.5, range2[0], 1e-6);
		assertEquals(0.5, range2[1], 1e-6);
	}

	/**
	 * @param ndata
	 * @param ndim
	 */
	public static QuantileExtraTrees getQET(int ndata, int ndim, boolean useWeights) {
		double[] output = new double[ndata];
		Matrix m = new Matrix(ndata, ndim);
		// generate values for all outputs
		for (int row=0; row < output.length; row++) {
			m.set(row, 1, row/(double)output.length);
			m.set(row, 2, 0.5);
			if (row == 5 || row == 6 || row == 7) {
				m.set(row, 2, Double.NaN);
			}
			output[row] = m.get(row, 1)>0.55 ?1 :0;
		}
		QuantileExtraTrees et = new QuantileExtraTrees(m, output);
		et.setHasNaN(true);
		if (useWeights) {
			double[] w = new double[ndata];
			for (int i=0; i < w.length; i++) {
				w[i] = 0.5;
			}
			et.setWeights(w);
		}
		return et;
	}

	@Test
	public void testQuantileNALearn() {
		int ndim = 5;
		QuantileExtraTrees et = getQET(100,  ndim, false);
		et.learnTrees(3, 3, 5);
		
		double[] x = new double[ndim];
		for (int i=0; i < x.length; i++) { 
			x[i] = Double.NaN;
		}
		double[] val;
		val = et.getQuantiles( new Matrix(x, 1, ndim), 0.5 );
		assertTrue( Double.isNaN(val[0]) );
	}

}
