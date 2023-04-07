package org.extratrees;

import java.util.ArrayList;

/**
 * Only used for leaves of QuantileExtraTrees.
 * Adds list of double values to the leaves.
 * 
 * @author Jaak Simm
 */
public class QuantileBinaryTree extends BinaryTree {
	ArrayList<Double> values;
}
