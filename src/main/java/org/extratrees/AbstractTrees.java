package org.extratrees;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public abstract class AbstractTrees<E> {
	ArrayList<E> trees;

	/** number of threads */
	int numThreads = 1;

	public int getNumThreads() {
		return numThreads;
	}
	
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	
	/**
	 * @return number of trees used in ExtraTrees
	 */
	public int getNumTrees() {
		return trees.size();
	}

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


}
