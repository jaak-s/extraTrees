package org.extratrees;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;

public class ExtraTrees extends AbstractTrees<BinaryTree> {
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
	 * @param tasks    - array of task indeces from 0 nTasks-1, null if no multi-task learning
	 */
	public ExtraTrees(Matrix input, double[] output, int[] tasks) {
		if (input.nrows!=output.length) {
			throw(new IllegalArgumentException("Input and output do not have same length."));
		}
		if (tasks!=null && input.nrows!=tasks.length) {
			throw(new IllegalArgumentException("Input and tasks do not have the same number of data points."));
		}
		this.input = input;
		this.output = output;
		this.outputSq = new double[output.length];
		for (int i=0; i<output.length; i++) {
			this.outputSq[i] = this.output[i]*this.output[i]; 
		}
		setTasks(tasks);
		
		// making cols list for later use:
		this.cols = new ArrayList<Integer>(input.ncols);
		for (int i=0; i<input.ncols; i++) {
			cols.add(i);
		}
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
	
	/** Average of several trees: */
	public static double getValue(ArrayList<BinaryTree> trees, double[] input) {
		double output = 0;
		for(BinaryTree t : trees) {
			output += t.getValue(input);
		}
		return output/trees.size();
	}

	/** Average of several trees, using nmin as depth */
	public static double getValue(ArrayList<BinaryTree> trees, double[] input, int nmin) {
		double output = 0;
		for(BinaryTree t : trees) {
			output += t.getValue(input, nmin);
		}
		return output/trees.size();
	}
	
	/**
	 * @param input
	 * @return matrix of predictions where
	 * output[i, j] gives prediction made for i-th row of input by j-th tree. 
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
	public double[] getValues(Matrix input) {
		return getValues(this.trees, input);
	}

	/** Average of several trees for many samples */
	public static double[] getValues(ArrayList<BinaryTree> trees, Matrix input) {
		double[] values = new double[input.nrows];
		double[] temp = new double[input.ncols];
		for (int row=0; row<input.nrows; row++) {
			// copying matrix row to temp:
			input.copyRow(row, temp);
			values[row] = getValue(trees, temp);
		}
		return values;
	}
	
	public double[] getValuesMT(Matrix newInput, int[] tasks) {
		double[] values = new double[newInput.nrows];
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

	public double getValueMT(double[] x, int task) {
		double mean = 0;
		for(BinaryTree t : trees) {
			mean += t.getValueMT(x, task);
		}
		mean /= trees.size();
		return mean;
	}

	/**
	 * @param input
	 * @return matrix of predictions where
	 * output[i, j] gives prediction made for i-th row of input by j-th tree.
	 */
	public Matrix getAllValuesMT(Matrix input, int[] tasks) {
		if (input.nrows!=tasks.length) {
			throw new IllegalArgumentException("Inputs and tasks do not have the same length.");
		}
		Matrix out = new Matrix( input.nrows, trees.size() );
		// temporary vector:
		double[] temp = new double[input.ncols];
		for (int row=0; row<input.nrows; row++) {
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
	protected void calculateCutScore(int[] ids, int col, double t,
			CutResult result) {
		// calculating score:
		double sumLeft=0, sumRight=0;
		double sumSqLeft=0, sumSqRight=0;
		for (int n=0; n<ids.length; n++) {
			if (input.get(ids[n], col) < t) {
				result.countLeft++;
				sumLeft   += output[  ids[n]];
				sumSqLeft += outputSq[ids[n]];
			} else {
				result.countRight++;
				sumRight   += output[  ids[n]];
				sumSqRight += outputSq[ids[n]];
			}
		}
		// calculating score:
		cutResultFromSums(result, sumLeft, sumRight, sumSqLeft, sumSqRight, result.countLeft, result.countRight);
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
	 */
	private void cutResultFromSums(CutResult result, double sumLeft,
			double sumRight, double sumSqLeft, double sumSqRight, 
			double countLeft, double countRight) {
		double varLeft  = sumSqLeft/countLeft  - 
				(sumLeft/countLeft)*(sumLeft/countLeft);
		double varRight = sumSqRight/countRight- 
				(sumRight/countRight)*(sumRight/countRight);
		// TODO: move var and var<zero*zero outside this loop:
		//double var = (sumSqLeft+sumSqRight)/ids.length - Math.pow((sumLeft+sumRight)/ids.length, 2.0);
		// the smaller the score the better:
		result.score = (result.countLeft*varLeft + result.countRight*varRight);// / ids.length / var;
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
		double alpha = 0.5;
		
		double[] scores = new double[nTasks];
		for (int i=0; i<ids.length; i++) {
			int n = ids[i];
			counts[tasks[n]] += 1;
			sums[tasks[n]]   += output[n];
			sumSq[tasks[n]]  += outputSq[n];
		}
		for (int task : nodeTasks) {
			// regularization:
			sums[task]   += priorMean*alpha;
			sumSq[task]  += priorMean*priorMean*alpha;
			regcounts[task] = counts[task] + alpha;
			// calculating regularized score:
			scores[task]  = sums[task] / (counts[task] + alpha);
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
		bt.value = 0;
		bt.nSuccessors = ids.length;
		bt.tasks = tasks;
		for (int n=0; n<ids.length; n++) {
			bt.value += output[ids[n]];
		}
		bt.value /= ids.length;
		return(bt);
	}
	
	
	
	/**
	 * @param ndata
	 * @param ndim
	 * @return
	 */
	
	
	public static ExtraTrees getSampleData(int ndata, int ndim) {
		double[] output = new double[ndata];
		double[] v = new double[ndata*ndim];
		for (int i=0; i<v.length; i++) {
			v[i] = Math.random();
		}
		Matrix m = new Matrix(v, ndata, ndim);
		// generate values for all outputs
		for (int row=0; row<output.length; row++) {
			m.set(row, 2, 0.5);
			output[row] = m.get(row, 1)+0.2*m.get(row, 3);
		}
		ExtraTrees et = new ExtraTrees(m, output);
		return et;
	}
	
	/**
	 * 
	 * @param trees
	 * @param testInput
	 * @param testOutput
	 * @return mean squared error on the test input and output.
	 */
	public static double getMeanSqError(ArrayList<BinaryTree> trees, Matrix testInput, double[] testOutput) {
		double error = 0;
		//double[] output_hat = getValues(trees, testInput);
		double[] temp = new double[testInput.ncols];
		for (int row=0; row<testInput.nrows; row++) {
			// copying matrix row to temp:
			for (int col=0; col<testInput.ncols; col++) {
				temp[col] = testInput.get(row, col);
			}
			error += Math.pow(getValue(trees, temp) - testOutput[row], 2);
		}

		return error/testOutput.length;
	}

	/** 
	 * mean squared error for CV
	 * @return mean squared error calculated only 
	 *         using data rows in testIds as test input-output.
	 **/
	public static double getMeanSqError(ArrayList<BinaryTree> trees, 
			Matrix testInput, double[] testOutput, 
			int nmin, int[] testIds)
	{
		double error = 0;
		//double[] output_hat = getValues(trees, testInput);
		double[] temp = new double[testInput.ncols];
		for (int n=0; n<testIds.length; n++) {
			int row = testIds[n];
			// copying matrix row to temp:
			for (int col=0; col<testInput.ncols; col++) {
				temp[col] = testInput.get(row, col);
			}
			error += Math.pow(getValue(trees, temp, nmin) - testOutput[row], 2);
		}

		return error/testIds.length;
	}

	public static double getMeanAbsError(ArrayList<BinaryTree> trees, Matrix testInput, double[] testOutput) {
		double error = 0;
		double[] output_hat = getValues(trees, testInput);
		for (int n=0; n<testOutput.length; n++) {
			error += Math.abs(output_hat[n] - testOutput[n]);
		}
		return error/testOutput.length;
	}
	
	public ArrayList<BinaryTree> buildTreeCV(int K, int nTrees) {
		int[] nmins = {2, 3, 5, 9, 14};
		int trainSize = (int)(2d/3d*output.length);
		
		// randomizing data:
		Integer[] ids = new Integer[this.output.length];
		for (int n=0; n<ids.length; n++) { ids[n]=n; }
		Collections.shuffle( Arrays.asList(ids) );
		
		// choosing trainSize of them as training data, others as test:
		int[] idsTrain = new int[trainSize];
		int[] idsTest  = new int[output.length - trainSize];
		for (int n=0; n<idsTrain.length; n++) {
			idsTrain[n] = ids[n];
		}
		for (int n=0; n<idsTest.length; n++) {
			idsTest[n]  = ids[n+idsTrain.length];
		}
		
		// building model:
		ArrayList<BinaryTree> t = this.buildTrees(2, K, nTrees, idsTrain);
		double[] errors = new double[nmins.length];
		double error_best = Double.POSITIVE_INFINITY;
		int nmin_best = nmins[0];
		for (int i=0; i<nmins.length; i++) {
			errors[i] = getMeanSqError(t, this.input, this.output, nmins[i], idsTest);
			if (errors[i]<error_best) {
				nmin_best = nmins[i];
				error_best = errors[i];
			}
		}
		//System.out.println( Arrays.toString(errors) );
		//System.out.println( String.format("Using nmin=%d based on CV.", nmin_best) );
		// building full model:
		ArrayList<BinaryTree> trees = this.buildTrees(nmin_best, K, nTrees);
		return trees;
	}
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int ndata  = 10000, ndim=7;
		int nTrees = 15;
		
		ExtraTrees et = getSampleData(ndata, ndim);
		Date t3 = new Date();
//		ArrayList<BinaryTree> m = et.buildTreeCV(ndim, nTrees);
		Date t4 = new Date();
		System.out.println( "Took: " + (t4.getTime()-t3.getTime())/1000.0 + "s");
		//int x=1;
		//if (x==1) return;
		Date t1 = new Date();
		et.learnTrees(2, 6, nTrees);
		ArrayList<BinaryTree> trees = et.trees;
		Date t2 = new Date();
		
		ExtraTrees et2 = getSampleData(1000, ndim);
		double[] output_hat = getValues(trees, et2.input);
		for (int row=0; row<et2.output.length; row++) {
			System.out.print( String.format("%d\t%1.3f %1.3f", row, et2.output[row], output_hat[row]) );
			System.out.println();
		}
		System.out.println( "Took: " + (t2.getTime()-t1.getTime())/1000.0 + "s");

		int[] ids = new int[et2.output.length];
		for (int n=0; n<ids.length; n++) { ids[n]=n; }
		double e1 = getMeanSqError(trees, et2.input, et2.output);
		double e2 = getMeanSqError(trees, et2.input, et2.output, 5, ids );
		System.out.println( "Error: " + e1 );
		System.out.println( "Error: " + e2 );
		

//		System.out.println( Arrays.toString( bt.countColumns(m.ncols)) );
	}


}
