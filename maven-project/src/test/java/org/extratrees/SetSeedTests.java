package org.extratrees;

import static org.junit.Assert.assertArrayEquals;

import org.junit.Test;

public class SetSeedTests {

	@Test
	public void testSetSeed() {
		ExtraTrees test = ExtraTreeTests.getSampleData(50, 5);
		
		int nmin = 3;
		int K = 2;
		int ntrees = 10;
		final long SEED = 1000;
		
		ExtraTrees et1 = ExtraTreeTests.getSampleData(100, 5);
		ExtraTrees et2 = new ExtraTrees(et1.input, et1.output);
		
		et1.setSeed(SEED);
		et2.setSeed(SEED);
		
		et1.learnTrees(nmin, K, ntrees);
		et2.learnTrees(nmin, K, ntrees);
		
		double[] yhat1 = et1.getValues(test.input);
		double[] yhat2 = et2.getValues(test.input);
		
		assertArrayEquals(yhat1, yhat2, 1e-8);
	}

	@Test
	public void testSetSeedMultiThreaded() {
		ExtraTrees test = ExtraTreeTests.getSampleData(50, 5);
		
		int nmin = 3;
		int K = 2;
		int ntrees = 10;
		final long SEED = 1000;
		
		ExtraTrees et1 = ExtraTreeTests.getSampleData(100, 5);
		et1.setNumThreads(2);
		ExtraTrees et2 = new ExtraTrees(et1.input, et1.output);
		
		et1.setSeed(SEED);
		et2.setSeed(SEED);
		
		et1.learnTrees(nmin, K, ntrees);
		et2.learnTrees(nmin, K, ntrees);
		
		double[] yhat1 = et1.getValues(test.input);
		double[] yhat2 = et2.getValues(test.input);
		
		assertArrayEquals(yhat1, yhat2, 1e-8);
	}

	@Test
	public void testSetSeedWithSubset() {
		ExtraTrees test = ExtraTreeTests.getSampleData(50, 5);
		
		int nmin = 3;
		int K = 2;
		int ntrees = 10;
		final long SEED = 1000;
		
		ExtraTrees et1 = ExtraTreeTests.getSampleData(100, 5);
		ExtraTrees et2 = new ExtraTrees(et1.input, et1.output);
		
		et1.setSubsetting(80);
		et2.setSubsetting(80);
		
		et1.setSeed(SEED);
		et2.setSeed(SEED);
		
		et1.learnTrees(nmin, K, ntrees);
		et2.learnTrees(nmin, K, ntrees);
		
		double[] yhat1 = et1.getValues(test.input);
		double[] yhat2 = et2.getValues(test.input);
		
		assertArrayEquals(yhat1, yhat2, 1e-8);
	}

}
