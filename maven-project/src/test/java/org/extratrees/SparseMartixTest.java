package org.extratrees;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Random;

import org.extratrees.data.CSparseMatrix;
import org.extratrees.data.Matrix;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class SparseMartixTest {

	@Test
	public void testCreate() {
		int[] r = {0,  2, 3, 3,  3, 5};
		int[] c = {10, 1, 2, 7, 13, 5};
		double[] v = {1, 0.5, 0, 1.5, -1.0, 1};
		int nrows = 10;
		int ncols = 14;
		
		Matrix m = new Matrix(nrows, ncols);
		for (int i = 0; i < r.length; i++) {
			m.set(r[i], c[i], v[i]);
		}
		
		CSparseMatrix sm = new CSparseMatrix(r, c, v, nrows, ncols);
		
		assertEquals( 0, sm.get(0, 0), 1e-9);
		for (int row = 0; row < nrows; row++) {
			for (int col = 0; col < ncols; col++) {
				assertEquals(m.get(row, col), sm.get(row, col), 1e-9);
			}
		}
	}
	
	@Test
	public void testNaNValues() {
		int[] r = {0,  2, 3, 3,  3, 5};
		int[] c = {10, 1, 2, 7, 13, 5};
		double[] v = {1, 0.5, 0, 1.5, -1.0, Double.NaN};
		int nrows = 10;
		int ncols = 14;
		
		Matrix m = new Matrix(nrows, ncols);
		for (int i = 0; i < r.length; i++) {
			m.set(r[i], c[i], v[i]);
		}
		
		CSparseMatrix sm = new CSparseMatrix(r, c, v, nrows, ncols);
		
		assertEquals( 0, sm.get(0, 0), 1e-9);
		for (int row = 0; row < nrows; row++) {
			for (int col = 0; col < ncols; col++) {
				assertEquals(m.get(row, col), sm.get(row, col), 1e-9);
			}
		}
	}

	
	@Rule
	public ExpectedException expectedEx = ExpectedException.none();

	@Test
	public void testTooSmallNrowException() {
		int[] r = {0,  2, 3, 3,  3, 5};
		int[] c = {10, 1, 2, 7, 13, 5};
		double[] v = {1, 0.5, 0, 1.5, -1.0, 1};
		
		expectedEx.expect(IllegalArgumentException.class);
		
		new CSparseMatrix(r, c, v, 5, 14);
		
	}

	@Test
	public void testTooSmallNcolException() {
		int[] r = {0,  2, 3, 3,  3, 5};
		int[] c = {10, 1, 2, 7, 13, 5};
		double[] v = {1, 0.5, 0, 1.5, -1.0, 1};

		expectedEx.expect(IllegalArgumentException.class);
		
		new CSparseMatrix(r, c, v, 6, 13);
	}
	
	public static FactorExtraTrees getSparseMatrix(int ndata) {
		int ncol      = 100;
		ArrayList<Integer> r = new ArrayList<Integer>();
		ArrayList<Integer> c = new ArrayList<Integer>();
		ArrayList<Double>  v = new ArrayList<Double> ();
		ArrayList<Integer> y = new ArrayList<Integer>();
		ShuffledIterator<Integer> cols = new ShuffledIterator<Integer>(
				AbstractTrees.getSequenceSet(ncol), new Random());
		Random rnd = new Random();
		
		for (int row = 0; row < ndata; row++) {
			cols.reset();
			int yval = row % 2;
			
			// number of random columns to add:
			int ncol_i = 10 + (int)(Math.random() * 10);
			
			for (int i = 0; i < ncol_i; i++) {
				int col = cols.next();
				double value = rnd.nextGaussian() + yval;
				r.add( row );
				c.add( col );
				v.add( value );
			}
			y.add( yval );
		}
		
		CSparseMatrix sm = new CSparseMatrix(
				FactorExtraTrees.listToArray(r),
				FactorExtraTrees.listToArray(c),
				ExtraTrees.listToDArray(v),
				ndata,
				ncol
		);
		
		assertEquals( r.size(), sm.getNumNonZero() );
		assertEquals( y.size(), sm.nrows() );
		int[] yi = AbstractTrees.listToArray(y);
		
		return new FactorExtraTrees(sm, yi );
	}
	
	@Test
	public void testLearningFromCSparseMatrix() {
		FactorExtraTrees fet = getSparseMatrix(200);
		FactorExtraTrees fet2 = getSparseMatrix(200);
		fet.learnTrees(1, 50, 100);
		
		int[] yhat = fet.getValues(fet2.input);
		
		int errors = 0;
		for (int i = 0; i < yhat.length; i++) {
			if (yhat[i] != fet2.output[i]) {
				errors++;
			}
		}
		double mean_error = errors / (double)yhat.length;
		
		System.out.println("Sparse Matrix (mean error): " + mean_error);
	}

}
