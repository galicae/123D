package sscc;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import loader.Converter;

import org.apache.commons.io.FileUtils;

import resources.AminoAcids;

public class SsccReader {
	private String ssccString;

	public SsccReader(String dir) {
		this.ssccString = dir;
	}

	/**
	 * 
	 * @return
	 * @throws IOException
	 */
	public SsccFile findAndRead() throws IOException {
		File ssccFile = new File("sscc/" + ssccString + ".sscc");
		ArrayList<Integer> sequence = new ArrayList<Integer>();
		StringBuilder structure = new StringBuilder();
		ArrayList<Integer> localContacts = new ArrayList<Integer>();
		ArrayList<Integer> globalContacts = new ArrayList<Integer>();
		for (String line : FileUtils.readLines(ssccFile)) {
			String[] lineArray = line.split("\\s+");
			int seq = AminoAcids.getIntRepresentation(lineArray[0]);
			char struc = lineArray[1].charAt(0);
			int lcont = Integer.parseInt(lineArray[2]);
			int gcont = Integer.parseInt(lineArray[3]);

			sequence.add(seq);
			structure.append(struc);
			localContacts.add(lcont);
			globalContacts.add(gcont);
		}
		Converter c = new Converter();
		int[] intStructure = c.convertSecStruc(structure.toString());

		int[] seqArray = new int[sequence.size()];
		for (int i = 0; i < seqArray.length; i++) {
			seqArray[i] = sequence.get(i);
		}

		int[] lcontArray = new int[localContacts.size()];
		for (int i = 0; i < lcontArray.length; i++) {
			lcontArray[i] = localContacts.get(i);
		}

		int[] gcontArray = new int[globalContacts.size()];
		for (int i = 0; i < gcontArray.length; i++) {
			gcontArray[i] = globalContacts.get(i);
		}
		
		SsccFile resultFile = new SsccFile(seqArray, lcontArray, gcontArray, intStructure,
				ssccString);
		return resultFile;
	}

}
