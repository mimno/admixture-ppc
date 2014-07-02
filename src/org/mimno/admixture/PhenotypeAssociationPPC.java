package org.mimno.admixture;

import java.util.*;
import java.util.zip.*;
import java.io.*;

public class PhenotypeAssociationPPC extends Replicator {

	public PhenotypeAssociationPPC(File dataFile, File genomePopFile, File posPopFile) throws IOException {
		super(dataFile, genomePopFile, posPopFile);
	}

	public double phenotypeBayesFactor(byte[][][] hiddenVariables, boolean replicating, int[] phenotypes,
										  double negativePrior, double positivePrior) {
		double maxBayesFactor = Double.NEGATIVE_INFINITY;
		double sumBayesFactor = 0;

		double[] prior = new double[] { positivePrior, negativePrior };

		int[][][] populationAllelePhenotypeCounts;
		
		byte[][] positionAlleles = new byte[numGenomes][2];

		for (int position = 0; position < maxNumPositions; position++) {
			populationAllelePhenotypeCounts = new int[ numPopulations ][2][2];

			// Pass through once to record or sample alleles for all genomes at this position
			for (int genome = 0; genome < numGenomes; genome++) {
				byte firstPopulation = hiddenVariables[genome][position][0];
				byte secondPopulation = hiddenVariables[genome][position][1];

				if (replicating) {
					// copy 1
					positionAlleles[genome][0] = random.nextDouble() < posPopMajorWeights[position][ firstPopulation ] ? (byte)1 : (byte)0;
					// copy 2
					positionAlleles[genome][1] = random.nextDouble() < posPopMajorWeights[position][ secondPopulation ] ? (byte)1 : (byte)0;
				}
				else {
					// both copies at the same time
					byte[] alleles = data.get(genome);
					positionAlleles[genome][0] = alleles[position] > 0 ? (byte)1 : (byte)0;
					positionAlleles[genome][1] = alleles[position] == 2 ? (byte)1 : (byte)0;
				}

				// Now update aggregate counts
				populationAllelePhenotypeCounts[ firstPopulation ][ positionAlleles[genome][0] ][ phenotypes[genome] ]++;
				populationAllelePhenotypeCounts[ secondPopulation ][ positionAlleles[genome][1] ][ phenotypes[genome] ]++;
			}

			
			double[][][] cache = new double[ numPopulations ][2][2];
			for (int population = 0; population < numPopulations; population++) {
				int[][] counts = populationAllelePhenotypeCounts[ population ];
				for (int allele = 0; allele < 2; allele++) {
					for (int phenotype = 0; phenotype < 2; phenotype++) {
						cache[population][allele][phenotype] = 
							Math.log( (prior[ phenotype ] + counts[allele][ phenotype ]) /
									  (prior[0] + counts[allele][0] + prior[1] + counts[allele][1]) ) -
							Math.log( (prior[ phenotype ] + counts[0][ phenotype ] + counts[1][ phenotype ]) /
									  (prior[0] + prior[1] + counts[0][0] + counts[0][1] + counts[1][0] + counts[0][1]) );
					}
				}
			}

			int numAlleles = 0;
			double bayesFactor = 0.0;

			// Now pass through again to calculate likelihood ratios
			for (int genome = 0; genome < numGenomes; genome++) {
                byte firstPopulation = hiddenVariables[genome][position][0];
				byte secondPopulation = hiddenVariables[genome][position][1];
				
				// probability of first allele 
				byte allele = positionAlleles[genome][0];
				int[][] counts = populationAllelePhenotypeCounts[ firstPopulation ];
				bayesFactor += cache[ firstPopulation ][ allele ][ phenotypes[genome] ];
				numAlleles++;
				
				/*
				System.out.format("%d / (%d + %d), (%d + %d) / (%d + %d + %d + %d) = %f\n",
								  counts[ phenotypes[genome] ][allele], 
								  counts[ phenotypes[genome] ][0], counts[ phenotypes[genome] ][1],
								  counts[0][allele], counts[1][allele],
								  counts[0][0], counts[0][1], counts[1][0], counts[0][1],
								  logLikelihoodRatio);
				*/

				
				// probability of the second allele
				allele = positionAlleles[genome][1];
				counts = populationAllelePhenotypeCounts[ secondPopulation ];
				bayesFactor += cache[ secondPopulation ][ allele ][ phenotypes[genome] ];
				numAlleles++;
				
			}

			sumBayesFactor += bayesFactor / numAlleles;
			if (bayesFactor / numAlleles > maxBayesFactor) {
				maxBayesFactor = bayesFactor / numAlleles;
			}
			
		}
		
		//return sumBayesFactor / (numPositions * Math.log(10));
		return maxBayesFactor / Math.log(10);
	}

	public void measureGWASBayesFactor(byte[][][] hiddenVariables, int numPhenotypeSamples) {

		double[][] sumReplicatedBayesFactors = new double[numReplications][numPopulations];
		double[] sumRealBayesFactors = new double[numPopulations];

		for (int phenotypeSample = 0; phenotypeSample < numPhenotypeSamples; phenotypeSample++) {
		
			// Generate a random phenotype with probability 0.5 for target population, 0.1 for others.
			int[][] populationPhenotypes = new int[numPopulations][numGenomes];
			
			for (int population = 0; population < numPopulations; population++) {
				for (int genome = 0; genome < numGenomes; genome++) {
					double[] populationWeights = genomePopulationWeights.get(genome);
					double probability = 0.5 * populationWeights[population] + 0.1 * (1.0 - populationWeights[population]);
					populationPhenotypes[population][genome] = random.nextDouble() < probability ? 1 : 0;
				}
			}
			
			for (int replication = 0; replication < numReplications; replication++) {
				for (int population = 0; population < numPopulations; population++) {
					double bayesFactor = phenotypeBayesFactor(hiddenVariables, true, populationPhenotypes[population],
															  0.1, 0.1);
					sumReplicatedBayesFactors[replication][population] += bayesFactor;
				}
			}
			
			for (int population = 0; population < numPopulations; population++) {
				double bayesFactor = phenotypeBayesFactor(hiddenVariables, false, populationPhenotypes[population],
														  0.1, 0.1);
				sumRealBayesFactors[population] += bayesFactor;
			}
		}

		for (int replication = 0; replication < numReplications; replication++) {
			for (int population = 0; population < numPopulations; population++) {
				System.out.format("%d\t%d\t%f\t%s\n", numPopulations, population, sumReplicatedBayesFactors[replication][population] / numPhenotypeSamples, "replicated");
			}
		}

		for (int population = 0; population < numPopulations; population++) {
			System.out.format("%d\t%d\t%f\t%s\n", numPopulations, population, sumRealBayesFactors[population] / numPhenotypeSamples, "real");
		}
	}
	
	public static void main(String[] args) throws Exception {
		PhenotypeAssociationPPC ppc = new PhenotypeAssociationPPC(new File(args[0]), new File(args[1]), new File(args[2]));
		byte[][][] hiddenVariables = ppc.generateHiddenVariables();
		ppc.measureGWASBayesFactor(hiddenVariables, 1);
	}

}