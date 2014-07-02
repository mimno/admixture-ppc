package org.mimno.admixture;

import java.util.*;
import java.util.zip.*;
import java.io.*;

public class Replicator {

	int numReplications = 100;
	int positionCutoff = Integer.MAX_VALUE;
		
	int maxNumPositions = 0;
	int numPositions = 0;
	int numPopulations = -1;
	int numGenomes = 0;

	AlleleData data;
	ArrayList<double[]> genomePopulationWeights;

	double[][] posPopMajorWeights;
	double[][] posPopMinorWeights;

	Random random;

	public Replicator(File dataFile, File genomePopFile, File posPopFile) throws IOException {
		this.numPopulations = numPopulations;
		random = new Random();

		genomePopulationWeights = new ArrayList<double[]>();
		
		String line;
		BufferedReader in = new BufferedReader(new FileReader(genomePopFile));
		while ((line = in.readLine()) != null) {
			String[] fields = line.split("\\s+");
			
			if (numPopulations == -1) {
				numPopulations = fields.length;
			}
			
			double[] weights = new double[numPopulations];
			for (int population = 0; population < numPopulations; population++) {
				weights[population] = Double.parseDouble(fields[population]);
			}
			genomePopulationWeights.add(weights);
		}
		in.close();
		
		numGenomes = genomePopulationWeights.size();

		System.err.println("read " + numGenomes + " genomes");


		// Read SNP-population weights into a buffer to so we can count 
		//  the total number of SNPs and create a statically-defined array later.
		ArrayList<double[]> snpPopBuffer = new ArrayList<double[]>();

		in = new BufferedReader(new FileReader(posPopFile));
        while ((line = in.readLine()) != null) {
            String[] fields = line.split("\\s+");
			
			double[] weights = new double[numPopulations];
			for (int population = 0; population < numPopulations; population++) {
				weights[population] = Double.parseDouble(fields[population]);
			}
			snpPopBuffer.add(weights);
		}
		in.close();
		
		numPositions = snpPopBuffer.size();
		maxNumPositions = numPositions;
		
		// Fill out the summary statistics
		posPopMajorWeights = new double[numPositions][numPopulations];
		posPopMinorWeights = new double[numPositions][numPopulations];

		for (int position = 0; position < numPositions; position++) {
			double[] weights = snpPopBuffer.get(position);
			for (int population = 0; population < numPopulations; population++) {
				posPopMajorWeights[position][population] = weights[population];
				posPopMinorWeights[position][population] = 1.0 - weights[population];
			}
		}

		System.err.println("read " + numPositions + " snp probs");
		
		data = new AlleleData();
		data.readPlinkBinary(dataFile, numGenomes, numPositions);

		System.err.println("read data");
	}
	
	public byte[][][] generateHiddenVariables() {
		// Generate two hidden population indicators for every position in every genome.
		byte[][][] hiddenVariables = new byte[numGenomes][maxNumPositions][2];

		double[] samplingWeights = new double[ numPopulations ];
		for (int genome = 0; genome < numGenomes; genome++) {
			byte[] alleles = data.get(genome);
			double[] populationWeights = genomePopulationWeights.get(genome);

			for (int position = 0; position < maxNumPositions; position++) {
				double sum = 0.0;
				if (alleles[position] == 0) {
					// Sample two hidden populations conditioned on both observations being 0, which we'll call the minor allele.
					sum = 0;
					for (int population = 0; population < numPopulations; population++) {
						samplingWeights[population] = populationWeights[population] * posPopMinorWeights[position][population];
						sum += samplingWeights[population];
					}
					
					double sample = random.nextDouble() * sum;
					byte pop = 0;
					while (sample > samplingWeights[pop]) {
						sample -= samplingWeights[pop];
						pop++;
					}
					hiddenVariables[genome][position][0] = pop;

					sample = random.nextDouble() * sum;
                    pop = 0;
                    while (sample > samplingWeights[pop]) {
						sample -= samplingWeights[pop];
						pop++;
                    }
                    hiddenVariables[genome][position][1] = pop;
				}
				else if (alleles[position] == 1) {
					// Sample two hidden populations, first conditioned on the observation
					//  being the major allele 1, the second conditioned on the observation
					//  being the minor allele 0.
					sum = 0;
					for (int population = 0; population < numPopulations; population++) {
						samplingWeights[population] = populationWeights[population] * posPopMajorWeights[position][population];
						sum += samplingWeights[population];
					}
					
					double sample = random.nextDouble() * sum;
					byte pop = 0;
					while (sample > samplingWeights[pop]) {
						sample -= samplingWeights[pop];
						pop++;
					}
					hiddenVariables[genome][position][0] = pop;
 
					// Recalculate the sampling distribution
					sum = 0;
					for (int population = 0; population < numPopulations; population++) {
						samplingWeights[population] = populationWeights[population] * posPopMinorWeights[position][population];
						sum += samplingWeights[population];
					}

					sample = random.nextDouble() * sum;
                    pop = 0;
                    while (sample > samplingWeights[pop]) {
						sample -= samplingWeights[pop];
						pop++;
                    }
                    hiddenVariables[genome][position][1] = pop;
				}
				else {
					// Sample two hidden populations conditioned on both observations being 1, which we'll call the major allele.
					sum = 0;
					for (int population = 0; population < numPopulations; population++) {
						samplingWeights[population] = populationWeights[population] * posPopMajorWeights[position][population];
						sum += samplingWeights[population];
					}
					
					double sample = random.nextDouble() * sum;
					byte pop = 0;
					while (sample > samplingWeights[pop]) {
						sample -= samplingWeights[pop];
						pop++;
					}
					hiddenVariables[genome][position][0] = pop;

					sample = random.nextDouble() * sum;
                    pop = 0;
                    while (sample > samplingWeights[pop]) {
						sample -= samplingWeights[pop];
						pop++;
                    }
                    hiddenVariables[genome][position][1] = pop;
				}
			}
		}

		return hiddenVariables;
	}

	public byte[] replicate(byte[][] genomeHiddenVariables, double[] populationWeights) {
		byte[] replicant = new byte[maxNumPositions];
		
		for (int position = 0; position < maxNumPositions; position++) {
			// The first copy
			if (random.nextDouble() < posPopMajorWeights[position][ genomeHiddenVariables[position][0] ]) {
				replicant[position]++;
			}
			
			// The second copy
			if (random.nextDouble() < posPopMajorWeights[position][ genomeHiddenVariables[position][1] ]) {
				replicant[position]++;
			}
		}

		return replicant;
	}

	public void replaceData(byte[][][] hiddenVariables) {
		for (int g = 0; g < numGenomes; g++) {
            byte[][] hidden = hiddenVariables[g];
			byte[] replicatedAlleles = replicate(hidden, genomePopulationWeights.get(g));
			data.set(g, replicatedAlleles);
		}
		System.err.println("read replacement data");
	}
	
	/** 
	 * As in Plink, I'm making the assumption that 0 is 00, 1 is 10, and 2 is 11.
	 */
	public boolean matches(byte a, byte b, int copy) {
		if (copy == 0) {
			if (a >= 1 && b >= 1) { return true; }
			if (a == 0 && b == 0) { return true; }
			return false;
		}
		else {
			if (a == 2 && b == 2) { return true; }
			if (a <= 1 && b <= 1) { return true; }
			return false;			
		}
	}

	public byte getCopy(byte x, int copy) {
		if (copy == 0) {
			if (x >= 1) { return 1; }
			else { return 0; }
		}
		else {
			if (x == 2) { return 1; }
			else { return 0; }
		}
	}

	public static void main(String[] args) throws Exception {

		/*
		
		Replicator replicator = new Replicator(numPopulationsOption.value,
											   new File(inputFile.value), 
											   new File(outputPrefix.value + ".genome_pops"),
											   new File(outputPrefix.value + ".snp_pops"));

		if (labelsFile.value != null) {
			replicator.loadLabels(new File(labelsFile.value), labelsColumn.value);
		}

		byte[][][] hiddenVariables = replicator.generateHiddenVariables();
		//replicator.compareExact(hiddenVariables);

		if (replicateFromReplicationOption.value) {
			// set the "real" data to a replication
			replicator.replaceData(hiddenVariables);
			// and recalculate the hidden population assignment variables
			hiddenVariables = replicator.generateHiddenVariables();
		}

		if (genomeSimilaritiesOption.value) {
			replicator.measureGenomeSimilarities(hiddenVariables);
		}
		else if (fstOption.value) {
			replicator.measureFST(hiddenVariables);
		}
		else if (gwasOption.value) {
			replicator.measureGWASBayesFactor(hiddenVariables, 10);
		}
		else if (entropyOption.value) {
			replicator.measurePosteriorEntropy(hiddenVariables);
		}
		else if (adjacentPairsOption.value) {
			replicator.countAdjacentAllelePairs(hiddenVariables, snpSpacingOption.value);
		}
		else {
			replicator.countAlleleMatrix(hiddenVariables, snpSpacingOption.value, snpWindowOption.value);
		}
		*/
	}
}
	
