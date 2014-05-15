package org.extratrees;

import java.util.Set;

public abstract class AbstractBinaryTree <T extends AbstractBinaryTree<T>> {
	/** tree for elements below threshold.
	 * if left==null, it is a leaf node
     * if left!=null, not a leaf
	 *  */
	public T left;
	/** tree for elements equal or above threshold. */
	public T right;

	/** number of elements in the tree */
	public int    nSuccessors;
	/** feature ID used for cutting */
	public int    column=-1;
	/** threshold of cutting */
	public double threshold; 

	/** tasks that are active in this thread */
	Set<Integer> tasks;
}
