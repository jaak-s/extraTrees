package org.extratrees;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.extratrees.data.Array2D;


public class FactorExtraTrees extends AbstractTrees<FactorBinaryTree, Integer> implements Serializable {
	private static final long serialVersionUID = 3625952360819832098L;
	
	transient int[] output;
	/** number of factors: */
	int nFactors;
	
	public FactorExtraTrees(int nFactors) {
		this.nFactors = nFactors;
	}

	public FactorExtraTrees(Array2D input, int[] output) {
		this(input, output, null);
	}
	
	/**
	 * @param input    - matrix of inputs, each row is an input vector
	 * @param output   - array of ints from 0 to nFactors-1 (class label)
	 * @param tasks    - array of task indeces from 0 nTasks-1, null if no multi-task learning
	 */
	public FactorExtraTrees(Array2D input, int[] output, int[] tasks) {
		if (input.nrows() != output.length) {
			throw(new IllegalArgumentException("Input and output do not have same length."));
		}
		if (tasks!=null && input.nrows() != tasks.length) {
			throw(new IllegalArgumentException("Input and tasks do not have the same number of data points."));
		}
		setInput(input);
		this.output = output;
		
		this.nFactors = 1;
		for (int i=0; i<output.length; i++) {
			if (output[i]<0) { 
				throw new RuntimeException("Bug: negative output (factor) values.");
			}
			if (nFactors<=output[i]) {
				nFactors = output[i]+1;
			}
		}
		setTasks(tasks);
	}
	
	public int getnFactors() {
		return nFactors;
	}
	
	public void setnFactors(int nFactors) {
		this.nFactors = nFactors;
	}
	
	/**
	 * @param selection
	 * @return new ExtraTrees object with the same input and output data with
	 * only the selected trees specified by {@code selection}.
	 */
	public FactorExtraTrees selectTrees(boolean[] selection) {
		FactorExtraTrees newET = new FactorExtraTrees(this.nFactors);
		newET.trees = new ArrayList<FactorBinaryTree>();
		for (int i=0; i<selection.length; i++) {
			if (!selection[i]) {
				continue;
			}
			newET.trees.add(this.trees.get(i));
		}
		return newET;
	}
	
	public class MajorityVote implements Aggregator<Integer> {
		int[] counts = new int[nFactors];
		
		@Override
		public void processLeaf(Integer leafValue) {
			counts[ leafValue ]++;
		}

		@Override
		public Integer getResult() {
			return getMaxIndex(counts);
		}
	}
	

	@Override
	Aggregator<Integer> getNewAggregator() {
		return new MajorityVote();
	}
	
	@Override
	double convertToDouble(Integer value) {
		return value >= 0 ? value : NA;
	}
	
	/**
	 * @param values
	 * @return index of the max positive value (first one if there are many) 
	 * or -1 if all values are non-positive.
	 */
	public static int getMaxIndex(int[] values) {
		int maxIndex = -1;
		int maxValue = 0;
		// adding counts:
		for (int i=0; i<values.length; i++) {
			if (values[i] > maxValue) {
				maxValue = values[i];
				maxIndex = i;
			}
		}
		return maxIndex;
	}

	/**
	 * @param values
	 * @return index of the max value (first one if there are many)
	 */
	public static int getMaxIndex(double[] values) {
		int maxIndex = -1;
		double maxValue = Double.NEGATIVE_INFINITY;
		// adding counts:
		for (int i=0; i<values.length; i++) {
			if (values[i]>maxValue) {
				maxValue = values[i];
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	
	private static int[] list2array(ArrayList<Integer> list) {
		int[] out = new int[ list.size() ];
		for (int i = 0; i < out.length; i++) {
			out[i] = list.get(i);
		}
		return out;
	}

	/**
	 * Object method, using the trees stored by learnTrees(...) method.
	 * @param input
	 * @return
	 */
	public int[] getValues(Array2D input) {
		return list2array( getValuesD(input) );
	}
	
	public int[] getValuesMT(Array2D newInput, int[] tasks) {
		return list2array( getValuesMTD(newInput, tasks) );
	}
	
	/**
	 * @param counts
	 * @return Gini index ( 1 - sum(f_i ^ 2) ), where f_i is 
	 *         the proportion of label i.<br>
	 *         Value 0 implies pure, high values mean noisy.  
	 */
	public static double getGiniIndex(double[] counts) {
		double sum = 0;
		double total = 0;
		for (int i=0; i<counts.length; i++) {
			sum   += counts[i]*counts[i];
			total += counts[i];
		}
		return 1 - sum / (double)( total*total );
	}
	

	/**
	 * @param counts
	 * @return Gini index ( 1 - sum(f_i ^ 2) ), where f_i is 
	 *         the proportion of label i.<br>
	 *         Value 0 implies pure, high values mean noisy.  
	 */
	public static double getGiniIndex(int[] counts) {
		long sum = 0;
		long total = 0;
		for (int i=0; i<counts.length; i++) {
			sum   += counts[i]*counts[i];
			total += counts[i];
		}
		return 1 - sum / (double)( total*total );
	}

	/**
	 * Makes tree that is filled with data.
	 * 
	 * @param leftTree
	 * @param rightTree
	 * @param col_best
	 * @param t_best
	 * @param nSuccessors
	 * @return
	 */
	@Override
	protected FactorBinaryTree makeFilledTree(FactorBinaryTree leftTree,
			FactorBinaryTree rightTree, 
			int col_best, double t_best, int nSuccessors) {
		FactorBinaryTree bt = new FactorBinaryTree();
		bt.column      = col_best;
		bt.threshold   = t_best;
		bt.nSuccessors = nSuccessors;
		bt.left  = leftTree;
		bt.right = rightTree;

		// TODO: add code for calculating value for intermediate nodes (for CV):
		// this feature needs refactoring
		//bt.value  = bt.left.value*bt.left.nSuccessors + bt.right.value*bt.right.nSuccessors;
		//bt.value /= bt.nSuccessors;

		return bt;
	}
	
	@Override
	protected TaskCutResult getTaskCut(int[] ids, Set<Integer> nodeTasks, double bestScore, int tree) {
		if (nFactors>2) {
			throw new RuntimeException("Multitask learning is not implemented 3 or more factors (classes).");
		}
		if (nodeTasks.size() <= 1) {
			return null;
		}

		int[][] factorTaskTable = getFactorTaskTable(ids);
		double[] p = getTaskScores(factorTaskTable);
		
		// counting if there are at least 2 tasks with samples
		if (! hasAtLeast2Tasks(nodeTasks, factorTaskTable) ) {
			return null;
		}
		
		double[] range = getRange(p);
		TaskCutResult bestResult = null;

		for (int repeat=0; repeat<this.numRandomTaskCuts; repeat++) {
			// get random cut:
			double t = getRandom(range[0], range[1], tree);
			TaskCutResult result = new TaskCutResult();
			calculateTaskCutScore(p, factorTaskTable, t, result);
			if (result.score < bestScore) {
				bestResult = result;
				bestScore = result.score;
			}
		}
		return bestResult;
	}

	protected boolean hasAtLeast2Tasks(Set<Integer> nodeTasks,
			int[][] factorTaskTable) {
		boolean has1Task = false;
		for (int task : nodeTasks) {
			if (factorTaskTable[0][task]>0 || factorTaskTable[1][task]>0) {
				if (has1Task) {
					return true;
				}
				has1Task = true;
			}
		}
		return false;
	}
	

	private int[][] getFactorTaskTable(int[] ids) {
		int[][] counts = new int[nFactors][nTasks];
		for (int i=0; i<nFactors; i++) {
			counts[i] = new int[nTasks];
		}
		for (int i=0; i<ids.length; i++) {
			int n = ids[i];
			counts[ output[n] ][ tasks[n] ]++;
		}
		return counts;
	}
	
	/**
	 * Assumes factors with 2 levels (binary classification).
	 * @param factorTaskTable
	 * @return
	 */
	private double[] getTaskScores(int[][] factorTaskTable) {
		int[][] counts = factorTaskTable;

		// calculate prior for regularization:
		double alpha  = 1;
		double[] sums = sumAlong2nd(counts);
		double prior  = (sums[0]+1)/(sums[0] + sums[1] + 2) * alpha;
		
		double[] scores = new double[nTasks];
		for (int task=0; task<nTasks; task++) {
			// regularized estimate for each task probability
			scores[task] = (counts[0][task] + prior) 
					     / (double)(counts[0][task] + counts[1][task] + alpha);
		}
		return scores;
	}
	
	/**
	 * @param factorTaskTable  array of size 2 (of arrays)
	 * @return
	 */
	public static double[] sumAlong2nd(int[][] factorTaskTable) {
		double[] sum = new double[2];
		for (int task=0; task<factorTaskTable[0].length; task++) {
			sum[0] += factorTaskTable[0][task];
			sum[1] += factorTaskTable[1][task];
		}
		return sum;
	}
	
	/**
	 * @return Gini index value given ids.
	 */
	@Override
	protected double get1NaNScore(int[] ids) {
		double[] factorCounts = new double[nFactors];
		for (int n=0; n<ids.length; n++) {
			int id = ids[n];
			factorCounts[ output[id] ] += (useWeights ?weights[id] :1.0);
		}
		return getGiniIndex(factorCounts);
	}
	
	/**
	 * Calculates the score for the cut. The smaller the better.
	 * @param ids
	 * @param col
	 * @param t
	 * @param result
	 */
	@Override
	protected void calculateCutScore(int[] ids, int col, double t,
			CutResult result) {
		if ( ! useWeights && ! hasNaN ) {
			int[][] factorCounts = new int[2][nFactors];
			for (int n=0; n<ids.length; n++) {
				factorCounts[input.get(ids[n], col) < t ?0 :1][ output[ids[n]] ]++;
			}
			result.countLeft  = sum(factorCounts[0]);
			result.countRight = sum(factorCounts[1]);
			double giniLeft  = getGiniIndex(factorCounts[0]);
			double giniRight = getGiniIndex(factorCounts[1]);
			result.score = (giniLeft*result.countLeft + giniRight*result.countRight) / (result.countLeft + result.countRight);
			result.leftConst  = giniLeft  < zero*zero;
			result.rightConst = giniRight < zero*zero;
		} else {
			// using weights, thus instead of counts we have weights, which are real numbers
			double[][] factorWeights = new double[2][nFactors];
			int[] branchCounts = new int[2];
			for (int n=0; n<ids.length; n++) {
				int id = ids[n];
				double value = input.get(id, col);
				double w = useWeights ?weights[id] :1.0;
				if (hasNaN) {
					if (Double.isNaN(value)) {
						result.nanWeigth += w;
						continue;
					}
				}
				int branch = value < t ?0 :1;
				factorWeights[ branch ][ output[id] ] += w;
				branchCounts[  branch ] ++;
			}
			result.countLeft  = branchCounts[0];
			result.countRight = branchCounts[1];
			double giniLeft  = getGiniIndex( factorWeights[0] );
			double giniRight = getGiniIndex( factorWeights[1] );
			double weightLeft  = sum(factorWeights[0]);
			double weightRight = sum(factorWeights[1]);
			result.score = (giniLeft*weightLeft + giniRight*weightRight) / (weightLeft + weightRight);
			result.leftConst  = giniLeft  < zero*zero;
			result.rightConst = giniRight < zero*zero;
		}
		
	}
	
	/**
	 * @param taskScores      scores
	 * @param factorTaskTable counts of factors and tasks
	 * @param t               cut
	 * @return GINI index for task cut: 1 - sum( (f_i)^2 )
	 */
	private void calculateTaskCutScore(double[] taskScores, int[][] factorTaskTable, double t, TaskCutResult result) {
		double[] leftCounts  = new double[nFactors];
		double[] rightCounts = new double[nFactors];
		result.leftTasks  = new HashSet<Integer>();
		result.rightTasks = new HashSet<Integer>();
		// TODO: make it use weights (or make an alternative function for that)		
		for (int task=0; task<factorTaskTable[0].length; task++) {
			if (taskScores[task] < t) {
				// task is going to the left branch
				for (int factor=0; factor<nFactors; factor++) {
					leftCounts[factor] += factorTaskTable[factor][task];
				}
				result.leftTasks.add(task);
			} else {
				// task is going to the right branch
				for (int factor=0; factor<nFactors; factor++) {
					rightCounts[factor] += factorTaskTable[factor][task];
				}
				result.rightTasks.add(task);
			}
		}
		// workaround for weights:
		result.countLeft  = (int)sum(leftCounts);
		result.countRight = (int)sum(rightCounts);
		cutResultFromCounts(result, leftCounts, rightCounts);
	}

	private void cutResultFromCounts(CutResult result, double[] leftCounts,
			double[] rightCounts) {
		double giniLeft  = getGiniIndex(leftCounts);
		double giniRight = getGiniIndex(rightCounts);
		//result.countLeft  = sum(leftCounts);
		//result.countRight = sum(rightCounts);

		result.score = (giniLeft*result.countLeft + giniRight*result.countRight) / (result.countLeft + result.countRight);
		result.leftConst  = giniLeft  < zero*zero;
		result.rightConst = giniRight < zero*zero;
	}


	/**
	 * @param ids
	 * @return builds a leaf node and returns it with the given ids.
	 */
	@Override
	public FactorBinaryTree makeLeaf(int[] ids, Set<Integer> tasks) {
		// terminal node:
		FactorBinaryTree bt = new FactorBinaryTree();
		bt.value = 0;
		bt.nSuccessors = ids.length;
		bt.tasks = tasks;
		// counting the factors:
		if ( ! useWeights) {
			int[] counts = new int[nFactors];
			for (int n=0; n<ids.length; n++) {
				counts[ output[ids[n]] ]++;
			}
			bt.value = getMaxIndex(counts);
		} else {
			// using weights to output the answer
			double[] counts = new double[nFactors];
			for (int n=0; n<ids.length; n++) {
				counts[ output[ids[n]] ] += weights[ ids[n] ];
			}
			bt.value = getMaxIndex(counts);
		}
		return(bt);
	}

}
