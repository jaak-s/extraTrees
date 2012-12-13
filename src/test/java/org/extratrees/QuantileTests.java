package org.extratrees;
import static org.junit.Assert.*;

import java.util.Arrays;

import org.extratrees.Matrix;
import org.junit.Test;

public class QuantileTests {
	public static QuantileExtraTrees getSampleData(int ndata, int ndim) {
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
		QuantileExtraTrees et = new QuantileExtraTrees(m, output);
		return et;
	}
	
	@Test
	public void testQuantile() {
		QuantileExtraTrees qet = getSampleData(200, 5);
		qet.learnTrees(5, 3, 100);
		QuantileExtraTrees qet2 = getSampleData(5, 5);
		double[] values = qet.getValues(qet2.input);
		double[] qs0_5 = qet.getQuantile( qet2.input, 0.5 );
		double[] qs0_6 = qet.getQuantile( qet2.input, 0.6 );
		double[] qs0_9 = qet.getQuantile( qet2.input, 0.9 );
		System.out.println( Arrays.toString(values) );
		System.out.println( Arrays.toString(qs0_5) );
		System.out.println( Arrays.toString(qs0_6) );
		System.out.println( Arrays.toString(qs0_9) );
	}

}
