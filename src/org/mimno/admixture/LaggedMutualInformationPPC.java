package org.mimno.admixture;

import java.util.*;
import java.util.zip.*;
import java.io.*;

public class LaggedMutualInformationPPC extends Replicator {

	public LaggedMutualInformationPPC(File dataFile, File genomePopFile, File posPopFile) throws IOException {
		super(dataFile, genomePopFile, posPopFile);
	}

	public int sum(int[] counts) {
		int n = 0;
		for (int i = 0; i < counts.length; i++) {
			n += counts[i];
		}
		return n;
	}

	public double mutualInformation(int[] counts) {
		double normalizer = 1.0 / (counts[0] + counts[1] + counts[2] + counts[3]);

		double p00 = counts[0] * normalizer;
		double p01 = counts[1] * normalizer;
		double p10 = counts[2] * normalizer;
		double p11 = counts[3] * normalizer;

		double p = p10 + p11;
		double q = p01 + p11;

		if (p == 1.0 || q == 1.0) { return 0.0; }

		double mi = 0.0;
		if (p00 > 0.0) {
			mi += p00 * Math.log(p00 / ((1.0 - p) * (1.0 - q)));
		}
		if (p01 > 0.0) {
			mi += p01 * Math.log(p01 / ((1.0 - p) * q));
		}
		if (p10 > 0.0) {
			mi += p10 * Math.log(p10 / (p * (1.0 - q)));
		}
		if (p11 > 0.0) {
			mi += p11 * Math.log(p11 / (p * q));
		}
		
		if (mi < 0.0) {
			mi = 0.0;
		}
		
		return mi;
	}

	public void countAllelePair(byte[][] hidden, byte[] alleles, int[][] distributions, int locus1, int locus2) {
		byte a0 = getCopy(alleles[locus1], 0);
		byte a1 = getCopy(alleles[locus1], 1);
		byte b0 = getCopy(alleles[locus2], 0);
		byte b1 = getCopy(alleles[locus2], 1);
		
		int popA0 = hidden[locus1][0];
		int popA1 = hidden[locus1][1];
		int popB0 = hidden[locus2][0];
		int popB1 = hidden[locus2][1];
		
		// I'm going to be adding up fractions of units (halves, quarters),
		//  but I want still to use int's instead of doubles to save memory.
		// Therefore I'm going to multiply everything by four. This will
		//  be washed out when I normalize the distribution.
		
		if (popA0 == popA1 && popA1 == popB0 && popB0 == popB1) {
			// First case, all populations are the same
			int population = popA0;
			
			distributions[ population ][ 2 * a0 + b0 ]++;
			distributions[ population ][ 2 * a0 + b1 ]++;
			distributions[ population ][ 2 * a1 + b0 ]++;
			distributions[ population ][ 2 * a1 + b1 ]++;
		}
		else if (popA0 == popA1) {
			// Two valid comparisons
			
			int population = popA0;
			if (popA0 == popB0) {
				distributions[ population ][ 2 * a0 + b0 ]++;
				distributions[ population ][ 2 * a1 + b0 ]++;
			}
			else if (popA0 == popB1) {
				distributions[ population ][ 2 * a0 + b1 ]++;
				distributions[ population ][ 2 * a1 + b1 ]++;
			}
		}			
		else if (popB0 == popB1) {
			// Two valid comparisons
			
			int population = popB0;
			if (popA0 == popB0) {
				distributions[ population ][ 2 * a0 + b0 ]++;
				distributions[ population ][ 2 * a0 + b1 ]++;
			}
			else if (popA1 == popB0) {
				distributions[ population ][ 2 * a1 + b0 ]++;
				distributions[ population ][ 2 * a1 + b1 ]++;
			}
		}
		else if (popA0 == popB0) {
			// One valid comparison
			int population = popA0;
			distributions[ population ][ 2 * a0 + b0 ]++;
		}
		else if (popA0 == popB1) {
			int population = popA0;
			distributions[ population ][ 2 * a0 + b1 ]++;
		}
		else if (popA1 == popB0) {
			int population = popA1;
			distributions[ population ][ 2 * a1 + b0 ]++;
		}
		else if (popA1 == popB1) {
			int population = popA1;
			distributions[ population ][ 2 * a1 + b1 ]++;
		}
	}

	public void countAlleleMatrix(byte[][][] hiddenVariables, int locusSpacing, int windowSize) {
		
		int numLags = windowSize / locusSpacing; // note integer division -- this will ignore remainder.
		int[] lags = new int[numLags];
		for (int i = 0; i < numLags; i++) {
			lags[i] = locusSpacing * (i+1);
		}
		
		double[][] populationMutualInformationSums;
		int[][] populationMutualInformationCounts;

		for (int replication = 0; replication < numReplications; replication++) {
			populationMutualInformationSums = new double[numPopulations][numLags];
			populationMutualInformationCounts = new int[numPopulations][numLags];

			byte[][] replicatedAlleles = new byte[numGenomes][numPositions];
			
			for (int g = 0; g < numGenomes; g++) {
				byte[][] hidden = hiddenVariables[g];
				replicatedAlleles[g] = replicate(hidden, genomePopulationWeights.get(g));
			}

			for (int locus1 = 0; locus1 < maxNumPositions; locus1 += locusSpacing) {
				for (int lagID = 0; lagID < lags.length; lagID++) {
					int lag = lags[lagID];
					int locus2 = locus1 + lag;
					
					if (locus2 >= maxNumPositions) { break; }
					
					int[][] distributions = new int[numPopulations][4];

					for (int g = 0; g < numGenomes; g++) {
						byte[][] hidden = hiddenVariables[g];
						byte[] alleles = replicatedAlleles[g];

						countAllelePair(hidden, alleles, distributions, locus1, locus2);
					}

					for (int population = 0; population < numPopulations; population++) {
						int totalPairs = sum(distributions[population]);
						if (totalPairs > 0) {
							populationMutualInformationSums[population][lagID] += 
								mutualInformation(distributions[population]);
							populationMutualInformationCounts[population][lagID] ++;
						}
					}
					
				}
			}
			
			for (int population = 0; population < numPopulations; population++) {
				for (int lagID = 0; lagID < lags.length; lagID++) {
					int lag = lags[lagID];
					System.out.format("%d\t%d\t%d\t%f\treplicated\n", numPopulations, population, lag,
						populationMutualInformationSums[population][lagID] / 
						populationMutualInformationCounts[population][lagID]);
				}
			}
		}
		
		populationMutualInformationSums = new double[numPopulations][numLags];
		populationMutualInformationCounts = new int[numPopulations][numLags];

		for (int locus1 = 0; locus1 < maxNumPositions; locus1 += locusSpacing) {
			for (int lagID = 0; lagID < lags.length; lagID++) {
				int lag = lags[lagID];
				int locus2 = locus1 + lag;
				
				if (locus2 >= maxNumPositions) { break; }

				int[][] distributions = new int[numPopulations][4];

				for (int g = 0; g < numGenomes; g++) {
					byte[][] hidden = hiddenVariables[g];
					byte[] alleles = data.get(g);

					countAllelePair(hidden, alleles, distributions, locus1, locus2);
				}

				for (int population = 0; population < numPopulations; population++) {
					int totalPairs = sum(distributions[population]);
					if (totalPairs > 0) {
						populationMutualInformationSums[population][lagID] += 
							mutualInformation(distributions[population]);
						populationMutualInformationCounts[population][lagID] ++;
					}
				}
				
			}
		}
		
		for (int population = 0; population < numPopulations; population++) {
			for (int lagID = 0; lagID < lags.length; lagID++) {
				int lag = lags[lagID];
				System.out.format("%d\t%d\t%d\t%f\treal\n", numPopulations, population, lag,
					populationMutualInformationSums[population][lagID] / 
					populationMutualInformationCounts[population][lagID]);
			}
		}
		
	}

	public static void main(String[] args) throws Exception {
		if (args.length == 5) {
			LaggedMutualInformationPPC ppc = new LaggedMutualInformationPPC(new File(args[0]), new File(args[1]), new File(args[2]));
			int locusSpacing = Integer.parseInt(args[3]);
			int windowSize = Integer.parseInt(args[4]);
			byte[][][] hiddenVariables = ppc.generateHiddenVariables();
			ppc.countAlleleMatrix(hiddenVariables, locusSpacing, windowSize);
		}
		else {
			System.out.println("Usage: ppc lag-mi [plink .bed file] [genome pop file] [snp pop file] [spacing] [window-size]");
		}
	}
}