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
		System.out.println(all);
		for (int row=0; row<yhat.length; row++) {
			double sum = 0;
			for (int j=0; j<all.ncols; j++) {
				sum += all.get(row, j);
			}
			assertEquals("row="+row, yhat[row], sum/all.ncols, 1e-6);
		}
	}

}
