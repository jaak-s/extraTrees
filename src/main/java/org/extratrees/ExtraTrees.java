package org.extratrees;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class ExtraTrees extends AbstractTrees<BinaryTree, Double> {
	double[] output;
	double[] outputSq;
	
	// defined in AbstractTrees:
	//ArrayList<BinaryTree> trees;
	public ExtraTrees(Matrix input, double[] output) {
		this(input, output, null);
	}
	

	/**
	 * @param input    - matrix of inputs, each row is an input vector
	 * @param output   - array of output values (doubles)
	 * @param tasks    - array of task indices from 0 nTasks-1, null if no multi-task learning
	 */
	public ExtraTrees(Matrix input, double[] output, int[] tasks) {
		if (input.nrows() != output.length) {
			throw(new IllegalArgumentException("Input and output do not have same length."));
		}
		if (tasks!=null && input.nrows() != tasks.length) {
			throw(new IllegalArgumentException("Input and tasks do not have the same number of data points."));
		}
		setInput(input);
		this.output = output;
		this.outputSq = new double[output.length];
		for (int i=0; i<output.length; i++) {
			this.outputSq[i] = this.output[i]*this.output[i];
		}
		setTasks(tasks);
	}
	
	/**
	 * @param selection
	 * @return new ExtraTrees object with the same input and output data with
	 * only the selected trees specified by {@code selection}.
	 */
	public ExtraTrees selectTrees(boolean[] selection) {
		ExtraTrees newET = new ExtraTrees(input, output);
		newET.trees = new ArrayList<BinaryTree>();
		for (int i=0; i<selection.length; i++) {
			if (!selection[i]) {
				continue;
			}
			newET.trees.add(this.trees.get(i));
		}
		return newET;
	}
	
	/** Builds trees with ids */
	public ArrayList<BinaryTree> buildTrees(int nmin, int K, int nTrees, int[] ids) {
		ArrayList<BinaryTree> trees = new ArrayList<BinaryTree>(nTrees);
		ShuffledIterator<Integer> cols = new ShuffledIterator<Integer>(this.cols);
		for (int t=0; t<nTrees; t++) {
			trees.add( this.buildTree(nmin, K, ids, cols, getSequenceSet(nTasks)) );
		}
		return trees;
	}
	
	public class ArithmeticMean implements Aggregator<Double> {
		double sum = 0;
		int count;
		
		@Override
		public void processLeaf(Double leafValue) {
			sum += leafValue;
			count++;
		}

		@Override
		public Double getResult() {
			if (count == 0) {
				return NA;
			}
			return sum/count;
		}
	}
	
	@Override
	Aggregator<Double> getNewAggregator() {
		return new ArithmeticMean();
	}
	
	/**
	 * @param input
	 * @return matrix of predictions where
	 * output[i, j] gives prediction made for i-th row of input by j-th tree. 
	 */
	public Matrix getAllValues(Matrix input) {
		Matrix out = new Matrix( input.nrows(), trees.size() );
		// temporary vector:
		double[] temp = new double[input.ncols()];
		for (int row=0; row<input.nrows(); row++) {
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
	public double[] getValues(Matrix input) {
		double[] values = new double[input.nrows()];
		double[] temp = new double[input.ncols()];
		for (int row=0; row<input.nrows(); row++) {
			// copying matrix row to temp:
			input.copyRow(row, temp);
			values[row] = getValue(temp);
		}
		return values;
	}

	public double[] getValuesMT(Matrix newInput, int[] tasks) {
		double[] values = new double[newInput.nrows()];
		double[] temp = new double[newInput.ncols()];
		for (int row=0; row<newInput.nrows(); row++) {
			// copying matrix row to temp:
			for (int col=0; col<newInput.ncols(); col++) {
				temp[col] = newInput.get(row, col);
			}
			values[row] = this.getValueMT(temp, tasks[row]);
		}
		return values;
	}

	/**
	 * @param input
	 * @return matrix of predictions where
	 * output[i, j] gives prediction made for i-th row of input by j-th tree.
	 */
	public Matrix getAllValuesMT(Matrix input, int[] tasks) {
		if (input.nrows() != tasks.length) {
			throw new IllegalArgumentException("Inputs and tasks do not have the same length.");
		}
		Matrix out = new Matrix( input.nrows(), trees.size() );
		// temporary vector:
		double[] temp = new double[input.ncols()];
		for (int row=0; row < input.nrows(); row++) {
			input.copyRow(row, temp);
			for (int j=0; j<trees.size(); j++) {
				out.set( row, j, trees.get(j).getValueMT(temp, tasks[row]) );
			}
		}
		return out;
	}

	

	@Override
	protected BinaryTree makeFilledTree(BinaryTree leftTree, BinaryTree rightTree,
			int col_best, double t_best, int nSuccessors) {
		BinaryTree bt = new BinaryTree();
		bt.column    = col_best;
		bt.threshold = t_best;
		bt.nSuccessors = nSuccessors;
		bt.left   = leftTree;
		bt.right  = rightTree;
		// value in intermediate nodes (used for CV):
		bt.value  = bt.left.value*bt.left.nSuccessors + bt.right.value*bt.right.nSuccessors;
		bt.value /= bt.nSuccessors;
		return bt;
	}
	
	@Override
	protected double get1NaNScore(int[] ids) {
		double sum=0, sumSq=0, weight=0;
		
		for (int n=0; n < ids.length; n++) {
			int id = ids[n];
			double w = useWeights ?weights[id] :1.0;
			weight += w;
			sum    += output[  id] * w;
			sumSq  += outputSq[id] * w;
		}

		double var = sumSq/weight - (sum/weight)*(sum/weight);
		return var;
	}

	/**
	 * result.score = weight(left) * var(left) + 
	 *                weight(right) * var(right)
	 */
	@Override
	protected void calculateCutScore(int[] ids, int col, double t,
			CutResult result) {
		// calculating score:
		double sumLeft=0, sumRight=0;
		double sumSqLeft=0, sumSqRight=0;
		double weightLeft=0, weightRight=0;
		for (int n=0; n<ids.length; n++) {
			int id = ids[n];
			double w = useWeights ?weights[id] :1.0;
			double value = input.get(id, col);
			if (hasNaN) {
				if (Double.isNaN(value)) {
					result.nanWeigth += w;
					continue;
				}
			}
			if (value < t) {
				result.countLeft++;
				weightLeft += w;
				sumLeft    += output[  id] * w;
				sumSqLeft  += outputSq[id] * w;
			} else {
				result.countRight++;
				weightRight += w;
				sumRight    += output[  id] * w;
				sumSqRight  += outputSq[id] * w;
			}
		}
		// calculating score:
		cutResultFromSums(result, sumLeft, sumRight, 
				sumSqLeft, sumSqRight, 
				weightLeft, weightRight);
		// value in intermediate nodes (used for CV):
	}

	/**
	 * 
	 * @param result
	 * @param sumLeft
	 * @param sumRight
	 * @param sumSqLeft
	 * @param sumSqRight
	 * @param countLeft   separate left  count (regularized in the case of task cut)
	 * @param countRight  separate right count (regularized in the case of task cut)
	 * @param countNaN    NaN count
	 * @param nanPenalty  penalty per NaN
	 */
	private void cutResultFromSums(CutResult result, double sumLeft,
			double sumRight, double sumSqLeft, double sumSqRight, 
			double countLeft, double countRight) {
		double varLeft  = sumSqLeft/countLeft  - 
				(sumLeft/countLeft)*(sumLeft/countLeft);
		double varRight = sumSqRight/countRight- 
				(sumRight/countRight)*(sumRight/countRight);
		result.score = (countLeft*varLeft + countRight*varRight);// / ids.length / var;
		result.leftConst  = (varLeft<zero*zero);
		result.rightConst = (varRight<zero*zero);
	}

	@Override
	protected TaskCutResult getTaskCut(int[] ids, Set<Integer> nodeTasks,
			double bestScore) {
		// return null if not at least 2 tasks
		if (nodeTasks.size() <= 1) {
			return null;
		}
		
		double mean  = getOutputMean(ids);
		int[] counts   = new int[nTasks];
		double[] regcounts= new double[nTasks];
		double[] sums  = new double[nTasks];
		double[] sumSq = new double[nTasks];
		double[] p = getTaskScores(ids, mean, nodeTasks, counts, regcounts, sums, sumSq);

		// check if there are at least two tasks
		if (! hasAtLeast2Tasks(ids) ) {
			return null;
		}
		
		double[] range = getRange(p);
		TaskCutResult bestResult = null;
		
		for (int repeat=0; repeat<this.numRandomTaskCuts; repeat++) {
			// get random cut:
			double t = getRandom(range[0], range[1]);
			TaskCutResult result = new TaskCutResult();
			calculateTaskCutScore(p, counts, regcounts, sums, sumSq, mean, t, result, nodeTasks);
			if (result.score < bestScore) {
				bestResult = result;
				bestScore  = result.score;
			}
		}
		return bestResult;
	}

	private void calculateTaskCutScore(double[] taskScores,
			int[] counts,
			double[] regcounts, 
			double[] sums,
			double[] sumSq,
			double mean, 
			double t,
			TaskCutResult result, 
			Set<Integer> nodeTasks) 
	{
		// TODO: support weights
		double sumLeft    = 0;
		double sumRight   = 0;
		double sumSqLeft  = 0;
		double sumSqRight = 0;
		double regcountLeft  = 0;
		double regcountRight = 0;
		result.leftTasks  = new HashSet<Integer>();
		result.rightTasks = new HashSet<Integer>();
		result.countLeft  = 0;
		result.countRight = 0;
		for (int task : nodeTasks) {
			if (taskScores[task] < t) {
				// left branch
				result.leftTasks.add(task);
				result.countLeft += counts[task];
				regcountLeft += regcounts[task];
				sumLeft   += sums[task];
				sumSqLeft += sumSq[task];
			} else {
				// right branch
				result.rightTasks.add(task);
				result.countRight += counts[task];
				regcountRight += regcounts[task];
				sumRight   += sums[task];
				sumSqRight += sumSq[task];
			}
		}
		cutResultFromSums(result, sumLeft, sumRight, sumSqLeft, sumSqRight, 
				regcountLeft, regcountRight);
	}

	private boolean hasAtLeast2Tasks(int[] ids) {
		int task0 = this.tasks[ids[0]];
		for (int i=1; i<ids.length; i++) {
			if (task0 != this.tasks[ids[i]]) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * 
	 * @param ids
	 * @param priorMean
	 * @param nodeTasks
	 * @param counts    filled by this method (NOT adjusted for regularization)
	 * @param regcounts filled by this method (adjusted for regularization)
	 * @param sums      filled by this method (adjusted for regularization)
	 * @param sumSq     filled by this method (adjusted for regularization)
	 * @return
	 */
	private double[] getTaskScores(int[] ids, double priorMean, Set<Integer> nodeTasks, 
			int[] counts,
			double[] regcounts,
			double[] sums, 
			double[] sumSq )
	{
		// calculate prior for regularization:
		double alpha = 1;
		
		double[] scores = new double[nTasks];
		for (int i=0; i<ids.length; i++) {
			int n = ids[i];
			counts[tasks[n]] += 1;
			sums[tasks[n]]   += output[n];
			sumSq[tasks[n]]  += outputSq[n];
		}
		for (int task : nodeTasks) {
			// only regularization to scores (sumsq, sums, regcounts are unaffected):
			regcounts[task] = counts[task];
			scores[task]  = (sums[task]+priorMean*alpha) / (counts[task] + alpha);
			
			// regularization:
			//sums[task]   += priorMean*alpha;
			//sumSq[task]  += priorMean*priorMean*alpha;
			// when regularization is used for comparison:
			//regcounts[task] = counts[task] + alpha;
			// calculating regularized score:
			//scores[task]  = sums[task] / (counts[task] + alpha);
		}
		
		return scores;
	}

	private double getOutputMean(int[] ids) {
		double mean = 0;
		
		for (int i=0; i<ids.length; i++) {
			mean += this.output[ids[i]];
		}
		mean /= ids.length;
		return mean;
	}
	
	/**
	 * @param ids
	 * @return builds a leaf node and returns it with the given ids.
	 */
	@Override
	public BinaryTree makeLeaf(int[] ids, Set<Integer> tasks) {
		// terminal node:
		BinaryTree bt = new BinaryTree();
		bt.value = 0d;
		bt.nSuccessors = ids.length;
		bt.tasks = tasks;
		double sumWeights = 0;
		for (int n=0; n<ids.length; n++) {
			double w = useWeights ?weights[ids[n]] :1.0;
			bt.value += output[ids[n]] * w;
			sumWeights += w;
		}
		bt.value /= sumWeights;
		return(bt);
	}

}
