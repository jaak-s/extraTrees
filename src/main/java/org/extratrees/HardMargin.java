package org.extratrees;

public class HardMargin {
	
	/**
	 * @param inputs NxD matrix of inputs (unlabeled)
	 * @param ids    array of IDs to use for calculation
	 * @param dim    dimension to cut
	 * @param cut    cutting point
	 * @return margin relative to the edges: margin / (max - min). 
	 * It is a value between 0.0 (no margin) and 0.5 (maximum margin).
	 * If max==min, returns 0.0.
	 */
	public double getCriteria(Matrix inputs, int[] ids, int dim, double cut) {
		double margin = Double.POSITIVE_INFINITY;
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		
		for (int i=0; i<ids.length; i++) {
			double val = inputs.get(ids[i], dim);
			double dist = Math.abs(cut-val);
			if (dist<margin) {
				margin = dist;
			}
			if (val < min) {
				min = val;
			}
			if (val > max) {
				max = val;
			}
		}
		if (min >= max) {
			// no points or all same:
			return 0;
		}
		return margin / (max-min);
	}
}
