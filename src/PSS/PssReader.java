package PSS;

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

	/**
	 * A wrapper function, to make things look a bit better on the outside
	 * 
	 * @return the Secondary Structure Preference scores in form of a double
	 *         array
	 * @throws IOException
	 */
	public double[][] calcPSSScore() throws IOException {
		PSSMeasurement m = measureFrequencies();
		return calcSeqStrucPref(m);
	}

	/**
	 * this function reads the sscc files and collects all data needed for the
	 * calculation of the secondary structure preference score, as defined in
	 * the paper. This measurements are saved as a PSSMeasurement object.
	 * 
	 * @return a {@link PSSMeasurement} object, containing all information
	 *         needed for calculation of the SSP score. This information is
	 *         N(i,s) [number of amino acids i in structure s], N(i) [number of
	 *         amino acids i], N(s) [number of amino acids in structure s] and N
	 *         [total number of amino acids]
	 * @throws IOException
	 */
	private PSSMeasurement measureFrequencies() throws IOException {
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

		return new PSSMeasurement(prefArray, occurAmino, occurStruct,
				totalAmino);
	}

	/**
	 * Performs the log odds calculation needed to provide the Secondary
	 * Structure Preference scores for all amino acids.
	 * 
	 * @param m
	 *            a {@link PSSMeasurement} object
	 * @return the score table, containing the PSS score for every amino acid
	 *         and every secondary structure (here alpha helix, beta strand and
	 *         coil)
	 */
	private double[][] calcSeqStrucPref(PSSMeasurement m) {
		double[][] secStrucPreference = new double[20][3];
		double[][] expectedSecStruc = new double[20][3];
		double n_i = 0;
		double n_s = 0;
		double up = 0;
		double n_total = (m.getTotalAmino() * 1.0);
		for (int i = 0; i < 20; i++) {
			for (int j = 0; j < 3; j++) {
				n_i = m.getOccurAmino()[i];
				n_s = m.getOccurStruct()[j];
				up = n_i * n_s;
				expectedSecStruc[i][j] = up;
			}
		}

		double odds = 0;
		for (int i = 0; i < 20; i++) {
			for (int j = 0; j < 3; j++) {
				odds = m.getPrefArray()[i][j] * n_total
						/ expectedSecStruc[i][j];
				// System.out.println(odds);
				secStrucPreference[i][j] = -Math.log(odds);
			}
		}
		return secStrucPreference;
	}

}
