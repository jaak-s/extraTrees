package org.extratrees;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.extratrees.AbstractTrees.CutResult;


public class FactorExtraTrees extends AbstractTrees<FactorBinaryTree> {
	int[] output;
	/** number of factors: */
	int nFactors;
	//String[] factorNames;
	
	//ArrayList<FactorBinaryTree> trees;

	public FactorExtraTrees(Matrix input, int[] output) {
		this(input, output, null);
	}
	
	/**
	 * @param input    - matrix of inputs, each row is an input vector
	 * @param output   - array of ints from 0 to nFactors-1 (class label)
	 * @param tasks    - array of task indeces from 0 nTasks-1, null if no multi-task learning
	 */
	public FactorExtraTrees(Matrix input, int[] output, int[] tasks) {
		if (input.nrows!=output.length) {
			throw(new IllegalArgumentException("Input and output do not have same length."));
		}
		if (tasks!=null && input.nrows!=tasks.length) {
			throw(new IllegalArgumentException("Input and tasks do not have the same number of data points."));
		}
		this.input  = input;
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
		// making a list of tasks:
		this.tasks  = tasks; 
		this.nTasks = 1;
		if (this.tasks!=null) {
			for (int i=0; i<tasks.length; i++) {
				//taskNames.add(tasks[i]);
				nTasks = tasks[i] + 1;
			}
		}
		
		// making cols list for later use:
		this.cols = new ArrayList<Integer>(input.ncols);
		for (int i=0; i<input.ncols; i++) {
			cols.add(i);
		}
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
		FactorExtraTrees newET = new FactorExtraTrees(input, output, tasks);
		newET.trees = new ArrayList<FactorBinaryTree>();
		for (int i=0; i<selection.length; i++) {
			if (!selection[i]) {
				continue;
			}
			newET.trees.add(this.trees.get(i));
		}
		return newET;
	}

	
	/** Builds trees with ids */
	/*
	public ArrayList<FactorBinaryTree> buildTrees(int nmin, int K, int nTrees, int[] ids) {
		ArrayList<FactorBinaryTree> trees = new ArrayList<FactorBinaryTree>(nTrees);
		ShuffledIterator<Integer> cols = new ShuffledIterator<Integer>(this.cols);
		for (int t=0; t<nTrees; t++) {
			trees.add( this.buildTree(nmin, K, ids, null, cols) );
		}
		return trees;		
	}*/

	/** Average of several trees: */
	public static int getValue(ArrayList<FactorBinaryTree> trees, double[] input, int nFactors) {
		int[] counts = new int[nFactors];
		for(FactorBinaryTree t : trees) {
			counts[ t.getValue(input) ]++;
		}
		return getMaxIndex(counts);
	}
	
	/** Average of several trees, using nmin as depth */
	/*
	public static double getValue(ArrayList<FactorBinaryTree> trees, double[] input, int nmin, int nFactors) {
		int[] counts = new int[nFactors];
		for(FactorBinaryTree t : trees) {
			counts[ t.getValue(input, nmin) ]++;
		}
		return getMaxIndex(counts);
	}*/

	/**
	 * @param values
	 * @return index of the max value (first one if there are many)
	 */
	public static int getMaxIndex(int[] values) {
		int maxIndex = -1;
		int maxValue = Integer.MIN_VALUE;
		// adding counts:
		for (int i=0; i<values.length; i++) {
			if (values[i]>maxValue) {
				maxValue = values[i];
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	
	/**
	 * @param input
	 * @return matrix of predictions where
	 * output[i, j] gives prediction made for i-th row of input by j-th tree.
	 * All values are integers. 
	 */
	public Matrix getAllValues(Matrix input) {
		Matrix out = new Matrix( input.nrows, trees.size() );
		// temporary vector:
		double[] temp = new double[input.ncols];
		for (int row=0; row<input.nrows; row++) {
			input.copyRow(row, temp);
			for (int j=0; j<trees.size(); j++) {
				out.set( row, j, trees.get(j).getValue(temp) );
			}
		}
		return out;
	}
	
	/**
	 * Object method, using the trees stored by learnTrees(...) method.
	 * @param input
	 * @return
	 */
	public int[] getValues(Matrix input) {
		return getValues(this.trees, input, this.nFactors);
	}
	
	/** Average of several trees for many samples */
	public static int[] getValues(ArrayList<FactorBinaryTree> trees, Matrix input, int nFactors) {
		int[]  values = new int[input.nrows];
		double[] temp = new double[input.ncols];
		for (int row=0; row<input.nrows; row++) {
			// copying matrix row to temp:
			for (int col=0; col<input.ncols; col++) {
				temp[col] = input.get(row, col);
			}
			values[row] = getValue(trees, temp, nFactors);
		}
		return values;
	}
	
	public int[] getValuesMT(Matrix newInput, int[] tasks) {
		int[] values = new int[newInput.nrows];
		double[] temp = new double[newInput.ncols];
		for (int row=0; row<newInput.nrows; row++) {
			// copying matrix row to temp:
			for (int col=0; col<newInput.ncols; col++) {
				temp[col] = newInput.get(row, col);
			}
			values[row] = this.getValueMT(temp, tasks[row]);
		}
		return values;
	}

	public int getValueMT(double[] x, int task) {
		int[] counts = new int[nFactors];
		for(FactorBinaryTree t : trees) {
			counts[ t.getValueMT(x, task) ]++;
		}
		return getMaxIndex(counts);
	}

	/**
	 * @param nmin - number of elements in leaf node
	 * @param K    - number of choices
	 *//*
	public FactorBinaryTree buildTree(int nmin, int K) {
		// generating full list of ids:
		int[]    ids = new int[input.nrows];
		for (int i=0; i<ids.length; i++) {
			ids[i] = i;
		}
		ShuffledIterator<Integer> cols = new ShuffledIterator<Integer>(this.cols);
		HashSet<Integer> taskSet = new HashSet<Integer>(tasks.length);
		for (int n : tasks) {
			taskSet.add(n);
		}
		return buildTree(nmin, K, ids, cols, new HashSet<Integer>() );
		//return buildTree(nmin, K, ids, unlabeledIds, cols);
	}*/
	
	/**
	 * @param counts
	 * @return Gini index ( 1 - sum(f_i ^ 2) ), where f_i is 
	 *         the proportion of label i.<br>
	 *         Value 0 implies pure, high values mean noisy.  
	 */
	public static double getGiniIndex(int[] counts) {
		int sum = 0;
		int total = 0;
		for (int i=0; i<counts.length; i++) {
			sum   += counts[i]*counts[i];
			total += counts[i];
		}
		return 1 - sum / (double)( total*total );
	}
	
	/*
	public FactorBinaryTree buildTree(int nmin, int K, int[] ids, 
			ShuffledIterator<Integer> randomCols) 
	{
		if (ids.length<nmin) {
			return makeLeaf(ids);
		}
		// doing a shuffle of cols:
		randomCols.reset();
		
		// trying K trees or the number of non-constant columns,
		// whichever is smaller:
		int k = 0, col_best=-1;
		double t_best=Double.NaN;
		CutResult bestResult = new CutResult();
		bestResult.score = Double.POSITIVE_INFINITY;
		
		while( randomCols.hasNext() ) {
			int col = randomCols.next();
			// calculating columns min and max:
			double[] range = getRange(ids, col, input);
			if (range[1]-range[0] < zero) {
				// skipping, because column is constant
				continue;
			}
			double diff = (range[1]-range[0]);
			for (int repeat=0; repeat<this.numRandomCuts; repeat++) {
				// picking random test point:
				double t;
				t = getRandomCut(range[0], diff, repeat);
			
				// calculating GINI impurity index (0 - pure, 1 - noisy):
				CutResult result = new CutResult();
				calculateCutScore(ids, col, t, result);
				
				if (result.score < bestResult.score) {
					col_best   = col;
					t_best     = t;
					
					bestResult.score = result.score;
					bestResult.leftConst  = result.leftConst;
					bestResult.rightConst = result.rightConst;
					bestResult.countLeft  = result.countLeft;
					bestResult.countRight = result.countRight;
				}
			}

			k++;
			if (k>=K) {
				// checked enough columns, stopping:
				break;
			}
		}
		// no score has been found, all inputs are constant:
		if (col_best<0) {
			return makeLeaf(ids);
		}

		// outputting the tree using the best score cut:
		int[] idsLeft  = new int[bestResult.countLeft];
		int[] idsRight = new int[bestResult.countRight];
		int nLeft=0, nRight=0;
		for (int n=0; n<ids.length; n++) {
			if (input.get(ids[n], col_best) < t_best) {
				// element goes to the left tree:
				idsLeft[nLeft] = ids[n];
				nLeft++;
			} else {
				// element goes to the right tree:
				idsRight[nRight] = ids[n];
				nRight++;
			}
		}

		// calculating sub trees:
		FactorBinaryTree leftTree, rightTree;
		if (bestResult.leftConst) { 
			// left child's output is constant
			leftTree = makeLeaf(idsLeft); 
		} else {
			leftTree  = this.buildTree(nmin, K, idsLeft, randomCols);
		}
		if (bestResult.rightConst) {
			// right child's output is constant
			rightTree = makeLeaf(idsRight);
		} else {
			rightTree = this.buildTree(nmin, K, idsRight, randomCols);
		}

		FactorBinaryTree bt = makeFilledTree(leftTree, rightTree,
				col_best, t_best, ids.length);
		return bt;
	}*/

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
	protected TaskCutResult getTaskCut(int[] ids, HashSet<Integer> tasks, double bestScore) {
		if (nFactors>2) {
			throw new RuntimeException("Multitask learning is not implemented 3 or more factors (classes).");
		}

		int[][] factorTaskTable = getFactorTaskTable(ids);
		double[] p = getTaskScores(factorTaskTable);
		
		// counting if there are at least 2 tasks with samples
		int numNonEmptyTasks = 0;
		for (int task=0; task<nTasks; task++) {
			if (factorTaskTable[0][task]>0 || factorTaskTable[1][task]>0) {
				numNonEmptyTasks++;
				if (numNonEmptyTasks>1) {
					break;
				}
			}
		}
		if (numNonEmptyTasks<=1) {
			// not enough tasks to do splitting
			return null;
		}
		
		double[] range = getRange(p);
		TaskCutResult bestResult = null;

		for (int repeat=0; repeat<this.numRandomTaskCuts; repeat++) {
			// get random cut:
			double t = getRandom(range[0], range[1]);
			TaskCutResult result = new TaskCutResult();
			calculateTaskCutScore(p, factorTaskTable, t, result);
			if (result.score < bestScore) {
				bestResult = result;
				bestScore = result.score;
			}
		}
		return bestResult;
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
	 * Assumes only 2 factors.
	 * @param factorTaskTable
	 * @return
	 */
	private double[] getTaskScores(int[][] factorTaskTable) {
		int[][] counts = factorTaskTable;
		double[] scores = new double[nTasks];
		for (int task=0; task<nTasks; task++) {
			// regularized estimate for each task probability
			scores[task] = (counts[0][task] + 1) 
					     / (double)(counts[0][task] + counts[1][task] + 2);
		}
		return scores;
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
		int[] factorCountLeft  = new int[nFactors];
		int[] factorCountRight = new int[nFactors];
		for (int n=0; n<ids.length; n++) {
			if (input.get(ids[n], col) < t) {
				//result.countLeft++;
				factorCountLeft[ output[ids[n]] ]++;
			} else {
				//result.countRight++;
				factorCountRight[ output[ids[n]] ]++;
			}
		}
		/* OLD CODE
		// calculating score:
		double giniLeft  = getGiniIndex(factorCountLeft);
		double giniRight = getGiniIndex(factorCountRight);
		
		result.score = (giniLeft*result.countLeft + giniRight*result.countRight) / ids.length;
		result.leftConst  = giniLeft  < zero*zero;
		result.rightConst = giniRight < zero*zero;*/
		cutResultFromCounts( result, factorCountLeft, factorCountRight );
	}
	
	/**
	 * @param taskScores      scores
	 * @param factorTaskTable counts of factors and tasks
	 * @param t               cut
	 * @return GINI index for task cut: 1 - sum( (f_i)^2 )
	 */
	private void calculateTaskCutScore(double[] taskScores, int[][] factorTaskTable, double t, TaskCutResult result) {
		int[] leftCounts  = new int[nFactors];
		int[] rightCounts = new int[nFactors];
		result.leftTasks  = new HashSet<Integer>();
		result.rightTasks = new HashSet<Integer>();
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
		cutResultFromCounts(result, leftCounts, rightCounts);
	}

	private void cutResultFromCounts(CutResult result, int[] leftCounts,
			int[] rightCounts) {
		double giniLeft  = getGiniIndex(leftCounts);
		double giniRight = getGiniIndex(rightCounts);
		result.countLeft  = sum(leftCounts);
		result.countRight = sum(rightCounts);

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
		int[] counts = new int[nFactors];
		for (int n=0; n<ids.length; n++) {
			counts[ output[ids[n]] ]++;
		}
		bt.value = getMaxIndex(counts);
		return(bt);
	}

}
