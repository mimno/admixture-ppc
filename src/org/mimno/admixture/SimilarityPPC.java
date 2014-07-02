package org.mimno.admixture;

import java.util.*;
import java.util.zip.*;
import java.io.*;

public class SimilarityPPC extends Replicator {

	public SimilarityPPC(File dataFile, File genomePopFile, File posPopFile) throws IOException {
		super(dataFile, genomePopFile, posPopFile);
	}
	
	public void countWeightedSharedAlleles(byte[][] hidden1, byte[][] hidden2, byte[] alleles1, byte[] alleles2, double[] populationSums, int[] populationCounts) {
		
		double[] populationSharedAlleleCounts = new double[numPopulations];
		int[] populationSharedPairs = new int[numPopulations];
		
		for (int position = 0; position < maxNumPositions; position++) {
			byte a0 = getCopy(alleles1[position], 0);
			byte a1 = getCopy(alleles1[position], 1);
			byte b0 = getCopy(alleles2[position], 0);
			byte b1 = getCopy(alleles2[position], 1);
			
			int popA0 = hidden1[position][0];
			int popA1 = hidden1[position][1];
			int popB0 = hidden2[position][0];
			int popB1 = hidden2[position][1];

			if (popA0 == popA1 && popA1 == popB0 && popB0 == popB1) {
				// First case, all populations are the same
				int population = popA0;
				populationSharedAlleleCounts[population] += 
					0.25 * ((a0 == b0 ? 1 : 0) + (a0 == b1 ? 1 : 0) + (a1 == b0 ? 1 : 0) + (a1 == b1 ? 1 : 0));
				populationSharedPairs[population]++;
			}
			else if (popA0 == popA1) {
				// Two valid comparisons

				int population = popA0;
				if (popA0 == popB0) {
					populationSharedAlleleCounts[population] += 
						0.5 * ((a0 == b0 ? 1 : 0) + (a1 == b0 ? 1 : 0));
				}
				else if (popA0 == popB1) {
					populationSharedAlleleCounts[population] += 
						0.5 * ((a0 == b1 ? 1 : 0) + (a1 == b1 ? 1 : 0));
				}
				populationSharedPairs[population]++;
			}			
			else if (popB0 == popB1) {
				// Two valid comparisons

				int population = popB0;
				if (popA0 == popB0) {
					populationSharedAlleleCounts[population] += 
						0.5 * ((a0 == b0 ? 1 : 0) + (a0 == b1 ? 1 : 0));
				}
				else if (popA1 == popB0) {
					populationSharedAlleleCounts[population] += 
						0.5 * ((a1 == b0 ? 1 : 0) + (a1 == b1 ? 1 : 0));
				}
				populationSharedPairs[population]++;
			}
			else if (popA0 == popB0) {
				// One valid comparison
				int population = popA0;
				populationSharedAlleleCounts[population] += (a0 == b0 ? 1 : 0);
				populationSharedPairs[population]++;
			}
			else if (popA0 == popB1) {
				int population = popA0;
				populationSharedAlleleCounts[population] += (a0 == b1 ? 1 : 0);
				populationSharedPairs[population]++;
			}
			else if (popA1 == popB0) {
				int population = popA1;
				populationSharedAlleleCounts[population] += (a1 == b0 ? 1 : 0);
				populationSharedPairs[population]++;
			}
			else if (popA1 == popB1) {
				int population = popA1;
				populationSharedAlleleCounts[population] += (a1 == b1 ? 1 : 0);
				populationSharedPairs[population]++;
			}
		}			

		for (int population = 0; population < numPopulations; population++) {
			if (populationSharedPairs[population] > 500) {
				populationSums[population] += populationSharedAlleleCounts[population] / populationSharedPairs[population];
				populationCounts[population]++;
			}
		}

	}

	public void measureGenomeSimilarities(byte[][][] hiddenVariables) {
		
		double[] populationSums;
		int[] populationCounts;
		
		for (int replication = 0; replication < numReplications; replication++) {
			populationSums = new double[ numPopulations ];
			populationCounts = new int[ numPopulations ];
			
			byte[][] replicatedAlleles = new byte[numGenomes][numPositions];
			
			for (int g = 0; g < numGenomes; g++) {
				byte[][] hidden = hiddenVariables[g];
				replicatedAlleles[g] = replicate(hidden, genomePopulationWeights.get(g));
			}

			for (int g1 = 0; g1 < numGenomes-1; g1++) {
				byte[][] hidden1 = hiddenVariables[g1];
				for (int g2 = g1+1; g2 < numGenomes; g2++) {
					byte[][] hidden2 = hiddenVariables[g2];
					countWeightedSharedAlleles(hidden1, hidden2, replicatedAlleles[g1], replicatedAlleles[g2], populationSums, populationCounts);
				}
			}
			
			for (int population = 0; population < numPopulations; population++) {
				System.out.format("%d\t%d\t%f\treplicated\n", numPopulations, population, populationSums[population] / populationCounts[population]);
			}
		}
				
		populationSums = new double[ numPopulations ];
		populationCounts = new int[ numPopulations ];
		
		for (int g1 = 0; g1 < numGenomes-1; g1++) {
			byte[][] hidden1 = hiddenVariables[g1];
			for (int g2 = g1+1; g2 < numGenomes; g2++) {
				byte[][] hidden2 = hiddenVariables[g2];
				countWeightedSharedAlleles(hidden1, hidden2, data.get(g1), data.get(g2), populationSums, populationCounts);
			}
		}
		
		for (int population = 0; population < numPopulations; population++) {
			System.out.format("%d\t%d\t%f\treal\n", numPopulations, population, populationSums[population] / populationCounts[population]);
		}
	}

	public static void main(String[] args) throws Exception {
		SimilarityPPC ppc = new SimilarityPPC(new File(args[0]), new File(args[1]), new File(args[2]));
		byte[][][] hiddenVariables = ppc.generateHiddenVariables();
		ppc.measureGenomeSimilarities(hiddenVariables);
	}
}