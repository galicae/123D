package gotoh;

import static resources.AminoAcids.REVERSE;
import static resources.AminoAcids.matrixNumbering;

import java.util.LinkedList;

import resources.Matrix;
import sscc.SsccFile;

public class LocalAligner extends Aligner {
	public LocalAligner(GotohProfile profile, int[] seq1, SsccFile sscc,
			String seq1ID) {
		super.profile = profile;
		super.seq1 = seq1;
		super.seq2 = sscc.getSequence();
		super.struct2 = sscc.getStructure();
		super.seq1ID = seq1ID;
		super.seq2ID = sscc.getID();
		super.localConts = sscc.getLocalContacts();
		super.globalConts = sscc.getGlobalContacts();
		super.score = new double[seq1.length + 1][seq2.length + 1];
		super.ins = new double[seq1.length + 1][seq2.length + 1];
		super.del = new double[seq1.length + 1][seq2.length + 1];
		tracebackList = new LinkedList<int[]>();
	}

	public void initializeF() {
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
	 * super function aligns the two sequences globally
	 */
	public void alignF() {
		for (int x = 1; x <= seq1.length; x++) {
			for (int y = 1; y <= seq2.length; y++) {
				double w1 = profile.getGextend()
						* profile.getGapExtend()[struct2[y - 1]]
						+ profile.getGopen()
						* profile.getGapInsert()[struct2[y - 1]];
				ins[x][y] = Math.max(score[x - 1][y] + w1, ins[x - 1][y]
						+ profile.getGextend()
						* profile.getGapExtend()[struct2[y - 1]]);
				del[x][y] = Math.max(score[x][y - 1] + w1, del[x][y - 1]
						+ profile.getGextend()
						* profile.getGapExtend()[struct2[y - 1]]);
				double temp = score[x - 1][y - 1] + match(x, y);
				score[x][y] = Math.max(Math.max(temp, 0),
						Math.max(ins[x][y], del[x][y]));
				// save the max; we need that for traceback (?)
				if (score[x][y] >= max[2]) {
					max[0] = x - 1;
					max[1] = y - 1;
					max[2] = score[x][y];
				}
			}
		}
	}

	private double match(int x, int y) {
		double seqScore = profile.getMatrixScore(matrixNumbering(seq1[x - 1]),
				matrixNumbering(seq2[y - 1]));
		double prefScore = Matrix.pss[struct2[y - 1]][seq1[x - 1]];
		double lcontScore = Matrix.ccp[struct2[y - 1]][seq1[x - 1]][localConts[y - 1]];
		double gcontScore = Matrix.ccp[struct2[y - 1]][seq1[x - 1]][globalConts[y - 1]];
		double result = profile.getLocalContWeight()[struct2[y - 1]]
				* lcontScore + profile.getGlobalContWeight()[struct2[y - 1]]
				* gcontScore + profile.getPrefWeight()[struct2[y - 1]]
				* prefScore + profile.getSeqWeight()[struct2[y - 1]] * seqScore;
		return result;
	}

	private static final double epsilon = 0.0001d;

	private static boolean isInEpsilon(double a, double b) {
		return (a > (b - epsilon)) && (a < (b + epsilon));
	}

	public void trace(int x, int y) {
		int[] tres = { x, y };
		tracebackList.push(tres);
		while (score[x][y] != 0) {
			if (isInEpsilon(score[x][y], ins[x][y])) {
				// find the x-k,y position where the gap was opened

				boolean found = false;
				while (!found) {
					if (ins[x - 1][y] + profile.getGextend() == ins[x][y]) {
						int[] res = { x - 1, y };
						tracebackList.push(res);
						x--;
						continue;
					} else {
						int[] res = { x - 1, y };
						x--;
						tracebackList.push(res);
						found = true;
					}
				}
			} else if (isInEpsilon(score[x][y], del[x][y])) {
				// find the x,y-k position where the gap was opened
				boolean found = false;
				while (!found) {
					if (del[x][y - 1] + profile.getGextend() == del[x][y]) {
						int[] res = { x, y - 1 };
						tracebackList.push(res);
						y--;
						continue;
					} else {
						int[] res = { x, y - 1 };
						y--;
						tracebackList.push(res);
						found = true;
					}
				}
			} else {// that means we came from score[x-1][y-1]
				int[] res = { x - 1, y - 1 };
				// checkScore += profile.getMatrixScore(seq1[x - 1], seq2[y -
				// 1]);
				tracebackList.push(res);
				x--;
				y--;
				continue;
			}
		}
	}

	public String[] interpretTraceback() {
		String[] result = new String[2];
		result[0] = result[1] = "";
		int[] temp = new int[2];
		int[] prev = new int[2];
		if (tracebackList.isEmpty()) {
			prev[0] = (int) max[0];
			prev[1] = (int) max[1];
		} else {
			prev = tracebackList.pop();
		}
		// super element has y=0, so I have to align every x before
		// prev[0],prev[1] with gaps, or x=0, so the other way round

		for (int i = 0; i < prev[0]; i++) {
			result[0] += Character.toString(REVERSE[(char) seq1[i]]);
			result[1] += "-";
		}

		for (int i = 0; i < prev[1]; i++) {
			result[1] += Character.toString(REVERSE[(char) seq2[i]]);
			result[0] += "-";
		}

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

		if (prev[0] != seq1.length && prev[1] != seq2.length) {
			result[0] += Character.toString(REVERSE[(char) seq1[prev[0]]]);
			result[1] += Character.toString(REVERSE[(char) seq2[prev[1]]]);
		}

		// now we are at the end of the freeshift; it remains to recover the
		// portion of the alignment that is parallel with the y axis

		for (int i = prev[0]; i <= seq1.length; i++) {
			result[0] += Character.toString(REVERSE[(char) seq1[i - 1]]);
			result[1] += "-";
		}

		for (int i = prev[1]; i <= seq2.length; i++) {
			result[1] += Character.toString(REVERSE[(char) seq2[i - 1]]);
			result[0] += "-";
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
		initializeF();
		alignF();
		trace((int) max[0], (int) max[1]);
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
