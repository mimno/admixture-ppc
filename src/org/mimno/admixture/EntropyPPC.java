package org.mimno.admixture;

import java.util.*;
import java.util.zip.*;
import java.io.*;

public class EntropyPPC extends Replicator {

	public EntropyPPC(File dataFile, File genomePopFile, File posPopFile) throws IOException {
		super(dataFile, genomePopFile, posPopFile);
	}

	public double[] meanAlleleEntropies(byte[][][] hiddenVariables, boolean replicating) {
		double[] populationAlleleEntropies = new double[numPopulations];
		int[] populationAlleleCounts = new int[numPopulations];

		for (int genome = 0; genome < numGenomes; genome++) {
			byte[] alleles = data.get(genome);
			double[] populationWeights = genomePopulationWeights.get(genome);
			for (int position = 0; position < maxNumPositions; position++) {
				byte firstPopulation = hiddenVariables[genome][position][0];
	            byte secondPopulation = hiddenVariables[genome][position][1];

				byte firstAllele, secondAllele;
			
				if (replicating) {
	                 // copy 1
	                 firstAllele = random.nextDouble() < posPopMajorWeights[position][ firstPopulation ] ? (byte)1 : (byte)0;
	                 // copy 2
	                 secondAllele = random.nextDouble() < posPopMajorWeights[position][ secondPopulation ] ? (byte)1 : (byte)0;
	             }
	             else {
	                 // both copies at the same time
	                 firstAllele = alleles[position] > 0 ? (byte)1 : (byte)0;
	                 secondAllele = alleles[position] == 2 ? (byte)1 : (byte)0;
				}

				for (int population = 0; population < numPopulations; population++) {
					double prob = 
						populationWeights[population] * 
						(firstAllele * posPopMajorWeights[position][population] +
						 (1 - firstAllele) * posPopMinorWeights[position][population]);
					if (prob > 0.0) {
						populationAlleleEntropies[firstPopulation] -= prob * Math.log(prob);
					}
					assert ! Double.isNaN(populationAlleleEntropies[firstPopulation]) : population + " " + prob;
				}
				populationAlleleCounts[firstPopulation]++;
			
				for (int population = 0; population < numPopulations; population++) {
					double prob = 
						populationWeights[population] * 
						(secondAllele * posPopMajorWeights[position][population] +
						 (1 - secondAllele) * posPopMinorWeights[position][population]);
					if (prob > 0.0) {
						populationAlleleEntropies[secondPopulation] -= prob * Math.log(prob);
					}
					assert ! Double.isNaN(populationAlleleEntropies[secondPopulation]) : population + " " + prob;
				}
				populationAlleleCounts[secondPopulation]++;
			}
		}

		for (int population = 0; population < numPopulations; population++) {
			populationAlleleEntropies[population] /= populationAlleleCounts[population];
			assert ! Double.isNaN(populationAlleleEntropies[population]) : "divided by " + populationAlleleCounts[population];
		}

		return populationAlleleEntropies;
	}

	public void measurePosteriorEntropy(byte[][][] hiddenVariables) {

		for (int replication = 0; replication < numReplications; replication++) {
			double[] meanEntropies = meanAlleleEntropies(hiddenVariables, true);
			for (int population = 0; population < numPopulations; population++) {
				System.out.format("%d\t%d\t%f\t%s\n", numPopulations, population, meanEntropies[population], "replicated");
			}
		}

		double[] meanEntropies = meanAlleleEntropies(hiddenVariables, false);
		for (int population = 0; population < numPopulations; population++) {
			System.out.format("%d\t%d\t%f\t%s\n", numPopulations, population, meanEntropies[population], "real");
		}
	}
 
	public static void main(String[] args) throws Exception {
		EntropyPPC ppc = new EntropyPPC(new File(args[0]), new File(args[1]), new File(args[2]));
		byte[][][] hiddenVariables = ppc.generateHiddenVariables();
		ppc.measurePosteriorEntropy(hiddenVariables);
	}
}