package org.extratrees;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public abstract class AbstractTrees<E> {
	protected final static double zero=1e-7;

	ArrayList<E> trees;

	/** number of threads */
	int numThreads = 1;
	
	/** number of times task cut is tried */
	int numRandomTaskCuts = 1;
	
	/** probability of trying task cutting (lambda) */
	double probOfTaskCuts = 1.0;

	/** number of random cuts tried for each feature */
	int numRandomCuts = 1;
	/** whether random cuts are totally uniform or evenly uniform */
	boolean evenCuts = false;


	// for semi-supervised learning:
	/** unsupervised points, used for semi-supervised learning */
	//Matrix unlabeled;
	/** weight of unsupervised learning (used to make semi-supervised signal) */
	//double weightOfUSL = 1.0;
	
//	HardMargin margin = new HardMargin();
	
	/*
	public void setUnlabeled(Matrix unlabeled) {
		this.unlabeled = unlabeled;
	}
	
	public Matrix getUnlabeled() {
		return unlabeled;
	}*/
	
//	public void setWeightOfUSL(double weightOfUSL) {
//		this.weightOfUSL = weightOfUSL;
//	}
//	
//	public double getWeightOfUSL() {
//		return weightOfUSL;
//	}

	public int getNumThreads() {
		return numThreads;
	}
	
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	
	public int getNumRandomTaskCuts() {
		return numRandomTaskCuts;
	}
	/** sets the number of times random cut for tasks is tried (assuming multitask setting).  
	 *  default is 1.
	 * */
	public void setNumRandomTaskCuts(int numRandomTaskCuts) {
		this.numRandomTaskCuts = numRandomTaskCuts;
	}
	
	public double getProbOfTaskCuts() {
		return probOfTaskCuts;
	}
	/**
	 * sets the probability of trying task cutting.
	 * @param probOfTaskCuts
	 */
	public void setProbOfTaskCuts(double probOfTaskCuts) {
		this.probOfTaskCuts = probOfTaskCuts;
	}
	
	/**
	 * @return number of trees used in ExtraTrees
	 */
	public int getNumTrees() {
		return trees.size();
	}
	
	public boolean isEvenCuts() {
		return evenCuts;
	}
	
	/**
	 * @param evenCuts - whether the random cuts (if more than 1) are
	 * sampled from fixed even intervals (true) 
	 * or just sampled ordinary uniform way (false)
	 */
	public void setEvenCuts(boolean evenCuts) {
		this.evenCuts = evenCuts;
	}
	
	public int getNumRandomCuts() {
		return numRandomCuts;
	}
	
	public void setNumRandomCuts(int numRandomCuts) {
		this.numRandomCuts = numRandomCuts;
	}

	/**
	 * @param col_min
	 * @param diff
	 * @param repeat  only used when evenCuts==true.
	 * @return random cut from col_min to col_min+diff.
	 */
	protected double getRandomCut(double col_min, double diff, int repeat) {
		double t;
		if (evenCuts) {
			double iStart = col_min + repeat*diff/numRandomCuts;
			double iStop  = col_min + (repeat+1)*diff/numRandomCuts;
			t = Math.random()*(iStop-iStart) + iStart;
		} else {
			t = Math.random()*diff + col_min;
		}
		return t;
	}
	
	/**
	 * @param ids
	 * @param col
	 * @param input
	 * @return array of min and max values.
	 */
	protected double[] getRange(int[] ids, int col, Matrix input) {
		double[] range = new double[2];
		range[0] = Double.POSITIVE_INFINITY;
		range[1] = Double.NEGATIVE_INFINITY;
		for (int n=0; n<ids.length; n++) {
			double v = input.get(ids[n], col);
			if ( v<range[0] ) { range[0] = v; }
			if ( v>range[1] ) { range[1] = v; }
		}
		return range;
	}

	
	/**
	 * @param inputs
	 * @param ids
	 * @param dim
	 * @param cut
	 * @return score of unsupervised error (the bigger the worse), i.e., 
	 * how far the margin is from the best margin (0.5).
	 * Gives 0 if unlabeled is null. 
	 * Uses <b>weightOfSSL</b> to adjust the result.
	 */
	/*
	public double calculateUSL(Matrix inputs, int[] ids, int dim, double cut) {
		if (unlabeled==null) {
			return 0;
		}
		// best margin is 0.5
		return weightOfUSL * ( 0.5 - margin.getCriteria(inputs, ids, dim, cut) );
	}
	*/
	
	/**
	 * Same as buildTrees() except computes in parallel.
	 * @param nmin
	 * @param K
	 * @param nTrees
	 * @return
	 */
	public ArrayList<E> buildTreesParallel(int nmin, int K, int nTrees) {
		// creating a thread pool and using it to compute nTrees:
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		List<TreeCallable> callables = new ArrayList<TreeCallable>();
		
		// adding tasks: using the same task for each tree
		TreeCallable task = new TreeCallable(nmin, K);
		for (int i=0; i<nTrees; i++) {
			callables.add( task );
		}
		// computing and fetching results:
		List<Future<E>> results;
		try {
			results = executor.invokeAll(callables);
		} catch (InterruptedException e) {
			// not solving this error here:
			throw new RuntimeException(e);
		}
		// fetching all BinaryTrees and storing them:
		ArrayList<E> trees = new ArrayList<E>(nTrees);
		for (Future<E> f : results) {
			try {
				trees.add( f.get() );
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			} catch (ExecutionException e) {
				throw new RuntimeException(e);
			}
		}
		executor.shutdown();
		return trees;
	}
	
	
	/** Nested class for making BinaryTrees */
	public class TreeCallable implements Callable<E> {
		int nmin, K;
		public TreeCallable(int nmin, int K) {
			this.nmin = nmin;
			this.K    = K;
		}
		@Override
		public E call() throws Exception {
			return AbstractTrees.this.buildTree(nmin, K);
		}
	}

	/**
 	 * good values:
	 * n_min = 2 (size of tree element)
	 * K = 5     (# of random choices)
	 * M = 50    (# of trees)
	 * if n_min is chosen by CV, then we have pruned version

	 * @param nmin   - size of tree element
	 * @param K      - # of random choices
	 * @param nTrees - # of trees
	 * Single threaded computation.
	 * @return learned trees
	 */
	public ArrayList<E> buildTrees(int nmin, int K, int nTrees) {
		ArrayList<E> trees = new ArrayList<E>(nTrees);
		// single-threading:
		for (int t=0; t<nTrees; t++) {
			trees.add( this.buildTree(nmin, K) );
		}
		return trees;
	}

	/**
	 * stores trees with the AbstractTrees object. 
	 * Uses multiple threads if set.
	 * @param nmin
	 * @param K
	 * @param nTrees
	 */
	public void learnTrees(int nmin, int K, int nTrees) {
		if (numThreads<=1) {	//if (numThreads<=1) {
			this.trees = buildTrees(nmin, K, nTrees);
		} else {
			this.trees = buildTreesParallel(nmin, K, nTrees);
		}
	}
	

	public abstract E buildTree(int nmin, int k);


	/**
	 * @param m    data matrix
	 * @param ids  row numbers to be split
	 * @param dim  which dimension to use for splitting
	 * @param cut  splitting value
	 * @return filters ids <b>ids</b> into two arrays, one whose values in data values
	 * {@code m[id, dim]} are below cut and others whose values are higher. 
	 */
	public static int[][] splitIds(Matrix m, int[] ids, int dim, double cut) {
		int[][] out = new int[2][];
		int lenLower = 0;
		for (int i=0; i<ids.length; i++) {
			if (m.get(ids[i], dim)<cut) {
				lenLower++;
			}
		}
		// two vectors: lower and higher
		out[0] = new int[lenLower];
		out[1] = new int[ids.length - lenLower];
		
		int i0 = 0;
		int i1 = 0;
		for (int i=0; i<ids.length; i++) {
			if (m.get(ids[i], dim)<cut) {
				out[0][i0] = ids[i];
				i0++;
			} else {
				out[1][i1] = ids[i];
				i1++;
			}
		}
		
		return out;
	}
}
