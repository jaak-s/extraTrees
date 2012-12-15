package org.extratrees;

import java.util.ArrayList;

public class QuantileExtraTrees extends ExtraTrees {

	public QuantileExtraTrees(Matrix input, double[] output) {
		super(input, output);
	}
	
	/**
	 * 
	 * @param input
	 * @param quantile a value between 0.0 and 1.0. For median use 0.5
	 * @return return quantiles for each input row. 
	 */
	public double[] getQuantiles(Matrix input, double k) {
		double[] quantileValues = new double[input.nrows];
		ArrayList<Double> leafValues = new ArrayList<Double>(this.trees.size());
		double[] temp = new double[input.ncols];
		for (int row=0; row<input.nrows; row++) {
			// copy row to temp:
			input.copyRow(row, temp);
			// get all leaf values for temp:
			getLeafValues(temp, leafValues);
			// doing quickselect:
			quantileValues[row] = QuickSelect.quickSelect(leafValues, k);
		}
		return quantileValues;
	}
	
	/**
	 * Clears ArrayList {@code values} and adds leaf values of {@code input} to it.
	 * @param input
	 * @param values 
	 */
	public void getLeafValues(double[] input, ArrayList<Double> values) {
		values.clear();
		for(BinaryTree t : trees) {
			QuantileBinaryTree leaf = (QuantileBinaryTree)t.getLeaf(input);
			values.addAll( leaf.values );
		}
	}
	

	/**
	 *  @param ids a list of ids in training data to make the leaf, 
	 *  stores double values in the node.  
	 *  @return leaf node
	 */
	@Override
	public BinaryTree makeLeaf(int[] ids) {
		// terminal node:
		QuantileBinaryTree bt = new QuantileBinaryTree();
		bt.value = 0;
		bt.nSuccessors = ids.length;
		bt.values = new ArrayList<Double>(ids.length);
		for (int n=0; n<ids.length; n++) {
			bt.value += output[ids[n]];
			// adding values for the leaf node:
			bt.values.add( output[ids[n]] );
		}
		bt.value /= ids.length;
		return(bt);
	}
	

}
