package org.extratrees;

public class BinaryTree extends AbstractBinaryTree<BinaryTree, Double> {
	public BinaryTree() {
	}
	
	@Override
	public Double getNA() {
		return AbstractTrees.NA;
	}
	
	/** count the number of times each column is used */
	public int[] countColumns(int ncols) {
		int[] counts = new int[ncols];
		this.countColumns( counts );
		return(counts);
	}
	
	/** counts the number of times each column is used recursively */
	public void countColumns(int[] counts) {
		if (left==null) {
			return;
		}
		counts[this.column]++;
		this.left.countColumns(counts);
		this.right.countColumns(counts);
	}
}
