package gotoh;

import static resources.AminoAcids.REVERSE;

import java.util.LinkedList;

import resources.Matrix;
import sscc.SsccFile;

public class FreeshiftAligner extends Aligner {
	private GotohProfile profile;
	private int[] seq1, seq2, struct2, cont2;
	public double[][] score;
	private double[][] ins, del;
	private double checkScore = 0;
	double[] max = { 0, 0, 0 };
	private LinkedList<int[]> tracebackList;
	private String seq1ID, seq2ID;

	public FreeshiftAligner(GotohProfile profile, int[] seq1, SsccFile sscc,
			String seq1ID) {
		this.profile = profile;
		this.seq1 = seq1;
		this.seq2 = sscc.getSequence();
		this.struct2 = sscc.getStructure();
		this.seq1ID = seq1ID;
		this.seq2ID = sscc.getID();
		this.score = new double[seq1.length + 1][seq2.length + 1];
		this.ins = new double[seq1.length + 1][seq2.length + 1];
		this.del = new double[seq1.length + 1][seq2.length + 1];
		tracebackList = new LinkedList<int[]>();
	}

	public void initialize() {
		score[0][0] = 0;
		for (int i = 1; i < seq1.length + 1; i++) {
			score[i][0] = 0;
			del[i][0] = 0;
		}
		for (int i = 1; i < seq2.length + 1; i++) {
			score[0][i] = 0;
			ins[0][i] = 0;
		}
		ins[0][0] = del[0][0] = 0;

	}

	/**
	 * this function aligns the two sequences globally
	 */
	public void align() {
		for (int x = 1; x <= seq1.length; x++) {
			for (int y = 1; y <= seq2.length; y++) {
				double w1 = profile.getGextend() + profile.getGopen();
				ins[x][y] = Math.max(score[x - 1][y] + w1, ins[x - 1][y]
						+ profile.getGextend());
				del[x][y] = Math.max(score[x][y - 1] + w1, del[x][y - 1]
						+ profile.getGextend());
				double temp = score[x - 1][y - 1] + match(x, y);
				score[x][y] = Math.max(temp, Math.max(ins[x][y], del[x][y]));
			}
		}
	}

	private double match(int x, int y) {
		double seqScore = profile.getMatrixScore(seq1[x - 1], seq2[y - 1]);
		double prefScore = Matrix.pss[seq1[x - 1]][struct2[y]];
		double contScore = Matrix.ccp[struct2[y]][cont2[y]][seq1[x - 1]];
		double result = profile.getContWeight() * contScore
				+ profile.getPrefWeight() * prefScore + profile.getSeqWeight()
				* seqScore;
		return result;
	}

	public void trace(int x, int y) {
		if (y == 0 || x == 0) {

		} else {
			int[] res = { x, y };
			checkScore += profile.getMatrixScore(seq1[x], seq2[y]);
			tracebackList.push(res);
			if (score[x][y] == ins[x][y]) {
				int k = 1;
				while (score[x - k][y] + profile.getGextend() * k
						+ profile.getGopen() != score[x][y]) {
					int[] resn = { x - k, y };
					checkScore += profile.getMatrixScore(seq1[x - k], seq2[y])
							+ profile.getGextend() * k + profile.getGopen();
					tracebackList.push(resn);
					if (x - k == 0)
						break;
					k++;
				}
				trace(x - k, y);
			} else if (score[x][y] == del[x][y]) {
				int k = 1;
				while (score[x][y - k] + profile.getGextend() * k
						+ profile.getGopen() != score[x][y]) {
					int[] resn = { x, y - k };
					checkScore += profile.getMatrixScore(seq1[x], seq2[y - k])
							+ profile.getGextend() * k + profile.getGopen();
					tracebackList.push(resn);
					if (y - k == 0)
						break;
					k++;
				}
				trace(x, y - k);
			} else {
				trace(x - 1, y - 1);
			}
		}
	}

	public String[] interpretTraceback() {
		String[] result = new String[2];
		result[0] = result[1] = "";
		int[] temp = new int[2];
		int[] prev = new int[2];
		prev = tracebackList.pop();
		// this element has y=0, so I have to align every x before
		// prev[0],prev[1] with gaps, or x=0, so the other way round
		if (prev[0] > prev[1]) {
			for (int i = 1; i < prev[0]; i++) {
				result[0] += Character.toString(REVERSE[(char) seq1[i - 1]]);
				result[1] += "-";
			}
		} else {
			for (int i = 1; i < prev[1]; i++) {
				result[1] += Character.toString(REVERSE[(char) seq2[i - 1]]);
				result[0] += "-";
			}
		}
		temp = tracebackList.pop();
		result[0] += Character.toString(REVERSE[(char) seq1[temp[0] - 1]]);
		result[1] += Character.toString(REVERSE[(char) seq2[temp[1] - 1]]);
		while (!tracebackList.isEmpty()) {
			temp = tracebackList.pop();
			if (temp[0] == prev[0]) {
				result[0] += "-";
				result[1] += Character
						.toString(REVERSE[(char) seq2[temp[1] - 1]]);
			} else if (temp[1] == prev[1]) {
				result[0] += Character
						.toString(REVERSE[(char) seq1[temp[0] - 1]]);
				result[1] += "-";
			} else {
				result[0] += Character
						.toString(REVERSE[(char) seq1[temp[0] - 1]]);
				result[1] += Character
						.toString(REVERSE[(char) seq2[temp[1] - 1]]);
			}
			prev = temp;
		}
		// now we are at the end of the freeshift; it remains to recover the
		// portion of the alignment that is parallel with the y axis
		if (prev[0] < prev[1]) {
			for (int i = prev[0]; i <= seq1.length; i++) {
				result[0] += Character.toString(REVERSE[(char) seq1[i - 1]]);
				result[1] += "-";
			}
		} else {
			for (int i = prev[1]; i <= seq2.length; i++) {
				result[1] += Character.toString(REVERSE[(char) seq2[i - 1]]);
				result[0] += "-";
			}
		}
		return result;
	}

	@Override
	public void printAlignment() {
		for (int y = 0; y <= seq2.length; y++) {
			for (int x = 0; x <= seq1.length; x++) {
				System.out.print(score[x][y] + "\t");
			}
			System.out.println();
		}
	}

	@Override
	public GotohAnswer alignPair() {
		initialize();
		align();

		for (int i = 0; i < seq1.length; i++) { // check last row for max
			if (score[i][seq2.length] >= max[2]) {
				max[0] = i;
				max[1] = seq2.length;
				max[2] = score[i][seq2.length];
			}
		}

		for (int i = 0; i < seq2.length; i++) { // check last column for max
			if (score[seq1.length][i] >= max[2]) {
				max[0] = seq1.length;
				max[1] = i;
				max[2] = score[seq1.length][i];
			}
		}

		trace((int) max[0] - 1, (int) max[1] - 1);
		String[] sresult = new String[2];
		sresult = interpretTraceback();

		GotohAnswer result = new GotohAnswer(seq1ID, seq2ID, sresult[0],
				sresult[1], max[2], profile);
		// System.out.println(seq1ID + " " + seq2ID + " " + max[2]);
		return result;
	}

	public double getCheckScore() {
		return this.checkScore;
	}
}
