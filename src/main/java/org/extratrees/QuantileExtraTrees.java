package org.extratrees;

import java.util.ArrayList;
import java.util.Set;

import org.extratrees.data.Array2D;
import org.extratrees.data.Row;

public class QuantileExtraTrees extends ExtraTrees {

	public QuantileExtraTrees(Array2D input, double[] output) {
		super(input, output);
	}
	
	/**
	 * 
	 * @param input
	 * @param quantile a value between 0.0 and 1.0. For median use 0.5
	 * @return return quantiles for each input row. 
	 */
	public double[] getQuantiles(Array2D input, double k) {
		double[] quantileValues = new double[input.nrows()];
		ArrayList<Double> leafValues = new ArrayList<Double>(this.trees.size());
		for (int row=0; row<input.nrows(); row++) {
			getLeafValues(input.getRow(row), leafValues);
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
	public void getLeafValues(Row input, ArrayList<Double> values) {
		values.clear();
		for(BinaryTree t : trees) {
			QuantileBinaryTree leaf = (QuantileBinaryTree)t.getLeaf(input);
			if (leaf != null) {
				values.addAll( leaf.values );
			}
		}
	}
	

	/**
	 *  @param ids a list of ids in training data to make the leaf, 
	 *  stores double values in the node.  
	 *  @return leaf node
	 */
	@Override
	public BinaryTree makeLeaf(int[] ids, Set<Integer> tasks) {
		// terminal node:
		QuantileBinaryTree bt = new QuantileBinaryTree();
		bt.value = 0d;
		bt.nSuccessors = ids.length;
		bt.tasks = tasks;
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
