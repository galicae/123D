package PSS;

import static resources.AminoAcids.*;

import java.io.File;
import java.io.IOException;
import java.util.Collection;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.TrueFileFilter;

import resources.AminoAcids;

public class PssReader {
	private String ssccString;
	private static final int ALPHA = 0, BETA = 1, OTHER = 2;

	public PssReader(String dir) {
		this.ssccString = dir;
	}

	public PSSResult calcSecStructPref() throws IOException {
		File ssccDir = new File(ssccString);
		Collection<File> ssccFiles = FileUtils.listFiles(ssccDir,
				TrueFileFilter.INSTANCE, TrueFileFilter.INSTANCE);

		int[][] prefArray = new int[20][3];
		int[] occurAmino = new int[20];
		int[] occurStruct = new int[3];
		int totalAmino = 0;

		for (File f : ssccFiles) {
			String[] file = FileUtils.readFileToString(f).split("\n");
			totalAmino += file.length;
			for (String s : file) {
				String[] line = s.split(" ");
				int currAmAcid = AminoAcids.getIntRepresentation(line[0]);
				occurAmino[currAmAcid]++;
				if (line[1].startsWith("a")) {
					prefArray[currAmAcid][ALPHA]++;
					occurStruct[ALPHA]++;
				} else if (line[1].startsWith("b")) {
					prefArray[currAmAcid][BETA]++;
					occurStruct[BETA]++;
				} else {
					prefArray[currAmAcid][OTHER]++;
					occurStruct[OTHER]++;
				}
			}
		}
		return new PSSResult(prefArray, occurAmino, occurStruct, totalAmino);
	}

}
