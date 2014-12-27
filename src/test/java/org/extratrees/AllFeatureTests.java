package org.extratrees;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ ExtraTreeTests.class, FactorTests.class,
		MultitaskTests.class, NATests.class,
		QuantileTests.class, QuickSelectTests.class,
		SetSeedTests.class,
		ShuffleTests.class,
		SparseMartixTest.class,
		SubsetTests.class})
public class AllFeatureTests {

}
