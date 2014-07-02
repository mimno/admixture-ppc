package org.mimno.admixture;

import java.util.*;
import java.util.zip.*;
import java.io.*;

public class FixationIndexPPC extends Replicator {
	HashMap<String,double[]> labelPopulationWeights = null;
	ArrayList<String> labels = null;

	public FixationIndexPPC(File dataFile, File genomePopFile, File posPopFile) throws IOException {
		super(dataFile, genomePopFile, posPopFile);
	}

	public void loadLabels(File labelsFile, int labelsColumn) throws IOException {
		labels = new ArrayList<String>();
		labelPopulationWeights = new HashMap<String,double[]>();

		BufferedReader in = new BufferedReader(new FileReader(labelsFile));
		String line = null;
		while ((line = in.readLine()) != null) {
			String[] fields = line.split("\\s+");
			String label = fields[ labelsColumn - 1 ];
			labels.add(label);
			if (! labelPopulationWeights.containsKey(label)) {
				labelPopulationWeights.put(label, new double[numPopulations]);
			}
		}
		in.close();
	}

	public double[] fixationIndex(byte[][][] hiddenVariables, boolean replicating, int[] individualCountries, int numCountries) {
		double[] populationFSTs = new double[ numPopulations ];
		int[] nonIdenticalSNPs = new int[ numPopulations ];

		int[][] populationCountryAlleles = new int[ numPopulations ][ numCountries ];
		int[][] populationCountryMajorAlleles = new int[ numPopulations ][ numCountries ];

		for (int position = 0; position < maxNumPositions; position++) {
			// Clear previous counts
			for (int population = 0; population < numPopulations; population++) {
				Arrays.fill(populationCountryAlleles[population], 0);
				Arrays.fill(populationCountryMajorAlleles[population], 0);
			}

			for (int genome = 0; genome < numGenomes; genome++) {
				byte[] alleles = data.get(genome);

				// First copy
				byte population = hiddenVariables[genome][position][0];
				int allele;
				if (replicating) {
					// Sample an allele
					allele = random.nextDouble() < posPopMajorWeights[position][ population ] ? 1 : 0;
				}
				else {
					// The real value
					allele = alleles[position] > 0 ? 1 : 0;
				}

				populationCountryAlleles[population][ individualCountries[genome] ]++;
				if (allele == 1) {
					populationCountryMajorAlleles[population][ individualCountries[genome] ]++;
				}

				// Second copy
				population = hiddenVariables[genome][position][1];
				if (replicating) {
					// Sample an allele
					allele = random.nextDouble() < posPopMajorWeights[position][ population ] ? 1 : 0;
				}
				else {
					// The real value
					allele = alleles[position] == 2 ? 1 : 0;
				}

				populationCountryAlleles[population][ individualCountries[genome] ]++;
				if (allele == 1) {
					populationCountryMajorAlleles[population][ individualCountries[genome] ]++;
				}
			}

			// Do population calculations

			for (int population = 0; population < numPopulations; population++) {
				double majorProbability = 0.0;
				int totalAlleles = 0;
                                
				// First get the average
				for (int country = 0; country < numCountries; country++) {
					majorProbability += populationCountryMajorAlleles[population][country];
					totalAlleles += populationCountryAlleles[population][country];
				}
                                
				majorProbability /= totalAlleles;

				// Now compute the variance
				double sumSquares = 0.0;
				int effectiveNumberOfCountries = 0;
				for (int country = 0; country < numCountries; country++) {
					if (populationCountryAlleles[population][country] > 0) {
						double diff = (double) populationCountryMajorAlleles[population][country] / populationCountryAlleles[population][country] - majorProbability;
						sumSquares += diff * diff;
						effectiveNumberOfCountries++;
					}
				}
                                
				if (majorProbability > 0.0 && majorProbability < 1.0 && effectiveNumberOfCountries > 1) {
					populationFSTs[population] += sumSquares / ((effectiveNumberOfCountries - 1) * majorProbability * (1.0 - majorProbability));
					nonIdenticalSNPs[population]++;
				}
			}
		}

		for (int population = 0; population < numPopulations; population++) {
			populationFSTs[population] /= nonIdenticalSNPs[population];
		}

		return populationFSTs;
	}

	public void measureFST(byte[][][] hiddenVariables) {

		int[] individualCountries = new int[ numGenomes ];
		HashMap<String,Integer> countryAlphabet = new HashMap<String,Integer>();
		for (int id = 0; id < numGenomes; id++) {
			String country = labels.get(id);
			int countryID;
			if (! countryAlphabet.containsKey(country)) {
				countryID = countryAlphabet.size();
				countryAlphabet.put(country, countryID);
			}
			else {
				countryID = countryAlphabet.get(country);
			}
			individualCountries[id] = countryID;
		}

		int numCountries = countryAlphabet.size();
		
		for (int replication = 0; replication < numReplications; replication++) {
			double[] fst = fixationIndex(hiddenVariables, true, individualCountries, numCountries);

			for (int population = 0; population < numPopulations; population++) {
				System.out.format("%d\t%d\t%f\t%s\n", numPopulations, population, fst[population], "replicated");
			}
		}

		double[] fst = fixationIndex(hiddenVariables, false, individualCountries, numCountries);
		
		for (int population = 0; population < numPopulations; population++) {
			System.out.format("%d\t%d\t%f\t%s\n", numPopulations, population, fst[population], "real");
		}
		
	}
	
	public static void main(String[] args) throws Exception {
		if (args.length == 4) {
			FixationIndexPPC ppc = new FixationIndexPPC(new File(args[0]), new File(args[1]), new File(args[2]));
			ppc.loadLabels(new File(args[3]), 1);
			byte[][][] hiddenVariables = ppc.generateHiddenVariables();
			ppc.measureFST(hiddenVariables);
		} 
		else if (args.length == 5) {
			FixationIndexPPC ppc = new FixationIndexPPC(new File(args[0]), new File(args[1]), new File(args[2]));
			ppc.loadLabels(new File(args[3]), Integer.parseInt(args[4]));
			byte[][][] hiddenVariables = ppc.generateHiddenVariables();
			ppc.measureFST(hiddenVariables);
		}
		else {
			System.out.println("Usage: ppc fst [plink .bed file] [genome pop file] [snp pop file] [labels file] [label column (default 1)]");
		}
	}
}