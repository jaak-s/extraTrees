package org.extratrees;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import org.junit.Test;

public class QuickSelectTests {

	@Test
	public void testSelect() {
		ArrayList<Double> a = new ArrayList<Double>(
				Arrays.asList(75.1, 40.0, 401.0, 10.0, 20.0, 500.0, -10.0, -10.0)
		);
		ArrayList<Double> b = new ArrayList<Double>(a);
		Collections.sort(b);
		
		for (int i=0; i<a.size(); i++) {
			assertEquals(b.get(i), QuickSelect.quickSelect(a, i+1), 1e-9);
		}
	}

}
