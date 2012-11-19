package gotoh;

import static resources.AminoAcids.REVERSE;
import static resources.AminoAcids.matrixNumbering;

import java.util.LinkedList;

import resources.Matrix;
import sscc.SsccFile;

public class LocalAligner extends Aligner {
	private GotohProfile profile;
	private int[] seq1, seq2, struct2, localConts, globalConts;
	public double[][] score;
	private double[][] ins, del;
	private double checkScore = 0;
	double[] max = { 0, 0, 0 };
	private LinkedList<int[]> tracebackList;
	private String seq1ID, seq2ID;

	public LocalAligner(GotohProfile profile, int[] seq1, SsccFile sscc,
			String seq1ID) {
		this.profile = profile;
		this.seq1 = seq1;
		this.seq2 = sscc.getSequence();
		this.struct2 = sscc.getStructure();
		this.seq1ID = seq1ID;
		this.seq2ID = sscc.getID();
		this.localConts = sscc.getLocalContacts();
		this.globalConts = sscc.getGlobalContacts();
		this.score = new double[seq1.length + 1][seq2.length + 1];
		this.ins = new double[seq1.length + 1][seq2.length + 1];
		this.del = new double[seq1.length + 1][seq2.length + 1];
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
	 * this function aligns the two sequences globally
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

	public void trace(int x, int y) {
		if (score[x][y] == 0) {

		} else {
			int[] res = { x, y };
			tracebackList.push(res);
			if (score[x][y] == ins[x][y]) {
				checkScore -= (profile.getGopen() * profile.getGapInsert()[struct2[seq2[y]]]);
				trace(x - 1, y);
			} else if (score[x][y] == del[x][y]) {
				trace(x, y - 1);
				checkScore -= (profile.getGopen() * profile.getGapInsert()[struct2[seq2[y]]]);
			} else {
				trace(x - 1, y - 1);
				checkScore -= profile.getMatrixScore(seq1[x - 1], seq2[y - 1]);
			}
		}
	}

	public String[] interpretTraceback() {
		String[] result = new String[2];
		result[0] = result[1] = "";
		int[] temp = new int[2];
		int[] prev = new int[2];
		temp = tracebackList.pop();
		// the traceback should only contain a limited number of elements. For
		// all sequence positions before prev[0] align seq1 with "-" and then
		// seq2 with "-"
		for (int i = 1; i < temp[0]; i++) {
			result[0] += Character.toString(REVERSE[(char) seq1[i - 1]]);
			result[1] += "-";
		}
		for (int i = 1; i < temp[1]; i++) {
			result[1] += Character.toString(REVERSE[(char) seq2[i - 1]]);
			result[0] += "-";
		}
		prev = temp;
		result[0] += Character.toString(REVERSE[(char) seq1[prev[0] - 1]]);
		result[1] += Character.toString(REVERSE[(char) seq2[prev[1] - 1]]);
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
		for (int i = prev[0]; i < seq1.length; i++) {
			result[0] += Character.toString(REVERSE[(char) seq1[i]]);
			result[1] += "-";
		}
		for (int i = prev[1]; i < seq2.length; i++) {
			result[1] += Character.toString(REVERSE[(char) seq2[i]]);
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
