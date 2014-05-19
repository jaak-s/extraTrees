package org.extratrees;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ BagTests.class, ExtraTreeTests.class, FactorTests.class,
		MultitaskTests.class, NATests.class,
		QuantileTests.class, QuickSelectTests.class, ShuffleTests.class })
public class AllFeatureTests {

}
