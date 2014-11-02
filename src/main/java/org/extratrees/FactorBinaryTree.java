package org.extratrees;

import java.io.Serializable;

public class FactorBinaryTree extends AbstractBinaryTree<FactorBinaryTree, Integer> implements Serializable {
	private static final long serialVersionUID = 6177152082032351501L;

	public FactorBinaryTree() {
	}
	
	@Override
	public Integer getNA() {
		return -1;
	}

	/** count the number of times each column is used */
	public int[] countColumns(int ncols) {
		int[] counts = new int[ncols];
		this.countColumns( counts );
		return(counts);
	}

	/** counts the number of times each column is used recursively */
	private void countColumns(int[] counts) {
		if (left==null) {
			return;
		}
		counts[this.column]++;
		this.left.countColumns(counts);
		this.right.countColumns(counts);
	}

}
