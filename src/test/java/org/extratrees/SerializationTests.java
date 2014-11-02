package org.extratrees;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import org.extratrees.data.Matrix;
import org.junit.Test;

public class SerializationTests {

	public static ExtraTrees getSampleData(int ndata, int ndim) {
		double[] output = new double[ndata];
		double[] v = new double[ndata * ndim];
		for (int i = 0; i < v.length; i++) {
			v[i] = Math.random();
		}
		Matrix m = new Matrix(v, ndata, ndim);
		// generate values for all outputs
		for (int row = 0; row < output.length; row++) {
			m.set(row, 2, 0.5);
			output[row] = m.get(row, 1) + 0.2 * m.get(row, 3);
		}
		ExtraTrees et = new ExtraTrees(m, output);
		return et;
	}

	@Test
	public void testEtSerialization() {
		try {
			File temp = File.createTempFile("et.ser", ".tmp");
			temp.deleteOnExit();

			FileOutputStream fileOut = new FileOutputStream(temp);
			ObjectOutputStream streamOut = new ObjectOutputStream(fileOut);

			ExtraTrees et = getSampleData(50, 5);
			et.learnTrees(5, 4, 10);

			streamOut.writeObject(et);
			streamOut.close();
			fileOut.close();
			System.out.printf("Serialized ET to temp file.");

			// deserializing
			ExtraTrees et2;
			FileInputStream fileIn = new FileInputStream(temp);
			ObjectInputStream in = new ObjectInputStream(fileIn);
			et2 = (ExtraTrees) in.readObject();
			in.close();
			fileIn.close();
			
			assertTrue(et2 != null);
			assertTrue(et2.trees != null);
			assertEquals(et.trees.size(), et2.trees.size());
			
			double[] yhat  = et.getValues(et.input);
			double[] yhat2 = et2.getValues(et.input);
			for (int i = 0; i < yhat.length; i++) {
				assertEquals(yhat[i], yhat2[i], 1e-6);
			}
			
		} catch (IOException e) {
			e.printStackTrace();
			return;
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			return;
		}
	}
	
	/**
	 * @param ndata
	 * @param ndim
	 * @return FactorExtraTrees, with factors of 0, 1, 2.
	 */
	public static FactorExtraTrees getFactorData(int ndata, int ndim) {
		int[] output = new int[ndata];
		double[] v   = new double[ndata*ndim];
		for (int i=0; i<v.length; i++) {
			v[i] = Math.random();
		}
		Matrix m = new Matrix(v, ndata, ndim);
		// generate values for all outputs
		for (int row=0; row<output.length; row++) {
			m.set(row, 2, 0.5);
			output[row] = (int) Math.floor(m.get(row, 1)+2*m.get(row, 3));
		}
		FactorExtraTrees et = new FactorExtraTrees(m, output);
		return et;
	}

	@Test
	public void testFactorSerialization() {
		try {
			File temp = File.createTempFile("fet.ser", ".tmp");
			temp.deleteOnExit();

			FileOutputStream fileOut = new FileOutputStream(temp);
			ObjectOutputStream streamOut = new ObjectOutputStream(fileOut);

			FactorExtraTrees fet = getFactorData(100, 5);
			fet.learnTrees(2, 4, 10);

			streamOut.writeObject(fet);
			streamOut.close();
			fileOut.close();
			System.out.printf("Serialized FactorET to temp file.");

			// deserializing
			FactorExtraTrees fet2;
			FileInputStream fileIn = new FileInputStream(temp);
			ObjectInputStream in = new ObjectInputStream(fileIn);
			fet2 = (FactorExtraTrees) in.readObject();
			in.close();
			fileIn.close();
			
			assertTrue(fet2 != null);
			assertTrue(fet2.trees != null);
			assertEquals(fet.trees.size(), fet2.trees.size());
			
			int[] yhat  = fet.getValues(fet.input);
			int[] yhat2 = fet2.getValues(fet.input);
			for (int i = 0; i < yhat.length; i++) {
				assertEquals(yhat[i], yhat2[i]);
			}
			
		} catch (IOException e) {
			e.printStackTrace();
			return;
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			return;
		}
	}

}
