package gotoh;

import static resources.AminoAcids.REVERSE;
import static resources.AminoAcids.matrixNumbering;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;

import resources.Matrix;
import sscc.SsccFile;

public class FreeshiftDouble extends Aligner {
	private GotohProfile profile;
	private int[] seq1, seq2, struct2, localConts, globalConts;
	public double[][] scoreF, scoreB;
	private double[][] insF, delF, insB, delB;
	private double checkScore = 0;
	double[] max = { 0, 0, 0 };
	private LinkedList<int[]> tracebackList;
	private String seq1ID, seq2ID;

	public FreeshiftDouble(GotohProfile profile, int[] seq1, SsccFile sscc,
			String seq1ID) {
		this.profile = profile;
		this.seq1 = seq1;
		this.seq2 = sscc.getSequence();
		this.struct2 = sscc.getStructure();
		this.seq1ID = seq1ID;
		this.seq2ID = sscc.getID();
		this.localConts = sscc.getLocalContacts();
		this.globalConts = sscc.getGlobalContacts();
		this.scoreF = new double[seq1.length + 1][seq2.length + 1];
		this.scoreB = new double[seq1.length + 1][seq2.length + 1];
		this.insF = new double[seq1.length + 1][seq2.length + 1];
		this.delF = new double[seq1.length + 1][seq2.length + 1];
		this.insB = new double[seq1.length + 1][seq2.length + 1];
		this.delB = new double[seq1.length + 1][seq2.length + 1];
		tracebackList = new LinkedList<int[]>();
	}

	public void initializeF() {
		scoreF[0][0] = 0;
		for (int i = 1; i < seq1.length + 1; i++) {
			scoreF[i][0] = 0;
			delF[i][0] = Double.NEGATIVE_INFINITY;
		}
		for (int i = 1; i < seq2.length + 1; i++) {
			scoreF[0][i] = 0;
			insF[0][i] = Double.NEGATIVE_INFINITY;
		}
		insF[0][0] = delF[0][0] = Double.NEGATIVE_INFINITY;

	}

	public void initializeB() {
		scoreB[0][0] = 0;
		for (int i = 1; i < seq1.length + 1; i++) {
			scoreB[i][0] = 0;
			delB[i][0] = Double.NEGATIVE_INFINITY;
		}
		for (int i = 1; i < seq2.length + 1; i++) {
			scoreB[0][i] = 0;
			insB[0][i] = Double.NEGATIVE_INFINITY;
		}
		insB[0][0] = delB[0][0] = Double.NEGATIVE_INFINITY;

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
				insF[x][y] = Math.max(scoreF[x - 1][y] + w1, insF[x - 1][y]
						+ profile.getGextend()
						* profile.getGapExtend()[struct2[y - 1]]);
				delF[x][y] = Math.max(scoreF[x][y - 1] + w1, delF[x][y - 1]
						+ profile.getGextend()
						* profile.getGapExtend()[struct2[y - 1]]);
				double temp = scoreF[x - 1][y - 1] + match(x, y);
				scoreF[x][y] = Math.max(temp, Math.max(insF[x][y], delF[x][y]));
			}
		}
	}

	public void alignB() {
		int seq1Length = this.seq1.length;
		int[] seq1 = new int[seq1Length];
		for (int i = 0; i < seq1Length; i++) {
			seq1[i] = this.seq1[seq1Length - i - 1];
		}

		int seq2Length = this.seq2.length;
		int[] seq2 = new int[seq2Length];
		for (int i = 0; i < seq2Length; i++) {
			seq2[i] = this.seq2[seq2Length - i - 1];
		}

		for (int x = 1; x <= seq1Length; x++) {
			for (int y = 1; y <= seq2.length; y++) {
				double w1 = profile.getGextend()
						* profile.getGapExtend()[struct2[y - 1]]
						+ profile.getGopen()
						* profile.getGapInsert()[struct2[y - 1]];
				insB[x][y] = Math.max(scoreB[x - 1][y] + w1, insB[x - 1][y]
						+ profile.getGextend()
						* profile.getGapExtend()[struct2[y - 1]]);
				delB[x][y] = Math.max(scoreB[x][y - 1] + w1, delB[x][y - 1]
						+ profile.getGextend()
						* profile.getGapExtend()[struct2[y - 1]]);
				double temp = scoreB[x - 1][y - 1] + match(x, y);
				scoreB[x][y] = Math.max(temp, Math.max(insB[x][y], delB[x][y]));
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
		while (x != 0 && y != 0) {
			if (scoreF[x][y] == insF[x][y]) {
				// find the x-k,y position where the gap was opened
				boolean found = false;
				int k = 1;
				while (!found) {
					double diff = scoreF[x - k][y] + profile.getGextend()
							* profile.getGapExtend()[struct2[y]] * (k - 1)
							+ profile.getGopen()
							* profile.getGapInsert()[struct2[y]];
					if (isInEpsilon(diff, scoreF[x][y]) && k < y) {
						int[] res = { x - k, y };
						tracebackList.push(res);
						k++;
					} else {
						int[] res = { x - k, y };
						tracebackList.push(res);
						found = true;
						x -= k;
					}
					continue;
				}
			} else if (scoreF[x][y] == delF[x][y]) {
				// find the x-k,y position where the gap was opened
				boolean found = false;
				int k = 1;
				while (!found) {
					double diff = scoreF[x][y - k] + profile.getGextend()
							* profile.getGapExtend()[struct2[y - k]] * (k - 1)
							+ profile.getGopen()
							* profile.getGapInsert()[struct2[y - k]];
					if (isInEpsilon(diff, scoreF[x][y]) && k < x) {
						int[] res = { x, y - k };
						tracebackList.push(res);
						k++;
					} else {
						int[] res = { x, y - k };
						tracebackList.push(res);
						found = true;
						y -= k;
					}
					continue;
				}
			} else {// that means we came from score[x-1][y-1]
				int[] res = { x - 1, y - 1 };
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
		if (!tracebackList.isEmpty()) {
			temp = tracebackList.pop();
			result[0] += Character.toString(REVERSE[(char) seq1[temp[0] - 1]]);
			result[1] += Character.toString(REVERSE[(char) seq2[temp[1] - 1]]);
		} else {
			result[0] += Character.toString(REVERSE[(char) seq1[prev[0] - 1]]);
			result[1] += Character.toString(REVERSE[(char) seq2[prev[1] - 1]]);
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
				System.out.print(scoreF[x][y] + "\t");
			}
			System.out.println();
		}
	}

	@Override
	public GotohAnswer alignPair() {
		initializeF();
		alignF();

		initializeB();
		alignB();

		for (int i = 0; i < seq1.length; i++) { // check last row for max
			// System.out.println(score[i][seq2.length]);
			if (scoreF[i][seq2.length] >= max[2]) {
				max[0] = i;
				max[1] = seq2.length;
				max[2] = scoreF[i][seq2.length];
			}
		}
		// System.out.println("##############");

		for (int i = 0; i < seq2.length; i++) { // check last column for max
			// System.out.println(score[seq1.length][i]);
			if (scoreF[seq1.length][i] >= max[2]) {
				max[0] = seq1.length;
				max[1] = i;
				max[2] = scoreF[seq1.length][i];
			}
		}

		// System.out.println(max[2]);
		try {
			FileWriter fstream = new FileWriter("outF.html");
			BufferedWriter out = new BufferedWriter(fstream);
			out.write(printMatriceHtml(scoreF, seq1, seq2));
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			FileWriter fstream = new FileWriter("outB.html");
			BufferedWriter out = new BufferedWriter(fstream);
			out.write(printMatriceHtml(scoreB, reverseSequence(seq1), reverseSequence(seq2)));
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		trace((int) max[0] - 1, (int) max[1] - 1);
		String[] sresult = new String[2];
		sresult = interpretTraceback();

		GotohAnswer result = new GotohAnswer(seq1ID, seq2ID, sresult[0],
				sresult[1], max[2], profile);

		return result;
	}

	public double getCheckScore() {
		return this.checkScore;
	}

	public String printMatriceHtml(double[][] matrix, int[] seq1, int[] seq2) {
		StringBuilder b = new StringBuilder();
		b.append("<!DOCTYPE html PUBLIC \"...\">");
		b.append("<html>");
		b.append("\t<head>");
		b.append("&nbsp;" + "&nbsp;" + "&nbsp;" + "&nbsp;" + "</head>");
		b.append("&nbsp;" + "&nbsp;" + "&nbsp;" + "&nbsp;" + "<body>"
				+ "<code>");
		b.append("<table border=\"1\"> <tr> <td>  </td>");
		b.append("<td>  </td>");
		for (int i = 0; i < matrix.length - 1; i++)
			b.append("<td><b>" + REVERSE[(char) seq1[i]] + "</b></td>");
		b.append("</tr><br>");
		int seq2Pointer = 0;
		for (int y = 0; y < matrix[0].length; y++) {
			b.append("<tr>");
			for (int x = 0; x < matrix.length; x++) {
				if (y == 0 && x == 0) {
					b.append("<td></td>");
				} else if (x == 0) {
					b.append("<td><b>" + REVERSE[seq2[seq2Pointer]]
							+ "</b></td>");
					seq2Pointer++;
				}
				double printout = Math.round(scoreF[x][y] * 10) / 10.0;
				b.append("<td>" + printout + "</td>");
			}
		}
		b.append("</code><br>");
		return b.toString();
	}

	private static final double epsilon = 0.0001d;

	private static boolean isInEpsilon(double a, double b) {
		return (a > (b - epsilon)) && (a < (b + epsilon));
	}

	private int[] reverseSequence(int[] revSeq) {
		int seq1Length = revSeq.length;
		int[] resSeq = new int[seq1Length];
		for(int i = 0; i < seq1Length; i++) {
			resSeq[i] = revSeq[seq1Length - i - 1];
		}
		return resSeq;
	}
}
