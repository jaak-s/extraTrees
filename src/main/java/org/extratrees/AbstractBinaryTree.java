package org.extratrees;

import java.util.Set;

public abstract class AbstractBinaryTree {

	/** number of elements in the tree */
	public int    nSuccessors;
	/** feature ID used for cutting */
	public int    column=-1;
	/** threshold of cutting */
	public double threshold; 

	/** tasks that are active in this thread */
	Set<Integer> tasks;
}
