package org.extratrees;

import static org.junit.Assert.*;

import org.junit.Test;

public class MarginCriteriaTest {

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
	public void test() {
		MarginCriteria mc = new MarginCriteria();
		int ndata = 21;
		int[] ids = new int[ndata];
		for (int i=0; i<ids.length; i++) {
			ids[i] = i;
		}
		
		Matrix inputs = getSampleData(ndata);
		assertEquals( 0.5/(ndata-1), mc.getCriteria(inputs, ids, 0, 5.5), 1e-6 );
		assertEquals( 0.5, mc.getCriteria(inputs, ids, 1, 1.5), 1e-6 );
		assertEquals( 0.4, mc.getCriteria(inputs, ids, 1, 1.4), 1e-6 );
		// testing selection:
		assertEquals( 0.5, mc.getCriteria(inputs, new int[]{0, ndata-1}, 0, (ndata-1)/2.0), 1e-6 );

	}

}
