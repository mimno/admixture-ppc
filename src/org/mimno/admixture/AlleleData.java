package org.mimno.admixture;

import java.io.*;
import java.nio.*;
import java.util.*;
import java.util.zip.*;

public class AlleleData {

	int MISSING_VALUE = -9;

	int numGenomes = 0;
	int numPositions = 0;

	public static final int SPACE_CODEPOINT = Character.codePointAt(" ", 0);
	public static final int ZERO_CODEPOINT = Character.codePointAt("0", 0);
	public static final int ONE_CODEPOINT = Character.codePointAt("1", 0);
	public static final int TWO_CODEPOINT = Character.codePointAt("2", 0);
	
	ArrayList<byte[]> genomePositionAlleles;

	public AlleleData() {
		genomePositionAlleles = new ArrayList<byte[]>();
	}
	
	public void readText(File dataFile) throws IOException {
		BufferedReader in = null;
		if (dataFile.getName().endsWith(".gz")) {
			//Process process = Runtime.getRuntime().exec("gzcat " + dataFile);
			//in = new BufferedReader(new InputStreamReader(process.getInputStream()));
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(dataFile))));
		}
		else {
			in = new BufferedReader(new FileReader(dataFile));
		}
		String line = null;
		
		long startTime = System.currentTimeMillis();
		int rowNumber = 1;
		while ((line = in.readLine()) != null) {
			if (numPositions == 0) { 
				numPositions = line.split("\\s+").length;
			}

			byte[] alleles = new byte[numPositions];

			int numCodePoints = Character.codePointCount(line, 0, line.length());
			int position = 0;
			for (int i = 0; i < numCodePoints; i++) {
				int codePoint = line.codePointAt(i);
				if (codePoint == ZERO_CODEPOINT) {
					alleles[position] = 0;
					position++;
				}
				else if (codePoint == ONE_CODEPOINT) {
					alleles[position] = 1;
					position++;
				}
				else if (codePoint == TWO_CODEPOINT) {
					alleles[position] = 2;
					position++;
				}
			}

			genomePositionAlleles.add(alleles);
			//if (rowNumber % 10 == 0) { System.out.format("%d %d\n", rowNumber, (System.currentTimeMillis() - startTime)); }
			rowNumber++;
		}

		numGenomes = genomePositionAlleles.size();
		System.err.format("Loaded %d individuals, %d SNPs\n", numGenomes, numPositions);
	}

	public void writeBinary(File outfile) throws IOException {
		DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile)));

		out.writeInt(numGenomes);
		out.writeInt(numPositions);
		
		for (int genome = 0; genome < numGenomes; genome++) {
			byte[] alleles = genomePositionAlleles.get(genome);
			for (int position = 0; position < numPositions; position++) {
				out.writeByte(alleles[position]);
			}
		}

		out.close();
	}

	public void readBinary(File infile) throws IOException {
		DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(infile)));
		
		numGenomes = in.readInt();
		numPositions = in.readInt();
		
		for (int genome = 0; genome < numGenomes; genome++) {
            byte[] alleles = new byte[numPositions];
            for (int position = 0; position < numPositions; position++) {
                alleles[position] = in.readByte();
            }
			genomePositionAlleles.add(alleles);
        }
		
		in.close();
	}
	
	public void readPlinkBinary(File infile, int numGenomes, int numPositions) throws IOException {
		this.numGenomes = numGenomes;
		this.numPositions = numPositions;

		int snp, genome;
		
		for (genome = 0; genome < numGenomes; genome++) {
			genomePositionAlleles.add(new byte[numPositions]);
		}

		FileInputStream in = new FileInputStream(infile);

		// Two bytes specifying format
		in.read();
		in.read();
		
		byte[][] data = new byte[numGenomes][numPositions];
		
		// Indicator specifying SNP x Individual (1) or Individual x SNP (0)
		int snpFirst = in.read();

		int c;
		
		snp = 0;
		genome = 0;
		
		
		while ((c = in.read()) != -1) {
			for (int i = 0; i < 4 && genome < numGenomes; i++) {
				genomePositionAlleles.get(genome)[snp] = (byte) (((c & 2) >> 1) + (c & 1));
				c = c >> 2;
				genome++;
			}
			
			if (genome == numGenomes) {
				genome = 0;
				snp++;
			}
		}
		
		in.close();
	}

	public byte[] get(int genome) {
		return genomePositionAlleles.get(genome);
	}

	public void set(int genome, byte[] alleles) {
		genomePositionAlleles.set(genome, alleles);
	}

	public static void main(String[] args) throws Exception {
		AlleleData data = new AlleleData();
		data.readText(new File(args[0]));
		data.writeBinary(new File(args[1]));
	}

}