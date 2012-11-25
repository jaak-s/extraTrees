import static org.junit.Assert.*;

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
		fet.learnTrees(2, 4, 10);
	}
}
