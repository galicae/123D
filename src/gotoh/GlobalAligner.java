package gotoh;

import static resources.AminoAcids.*;

import java.util.LinkedList;

import resources.Matrix;
import sscc.SsccFile;

public class GlobalAligner extends Aligner {
	public GlobalAligner(GotohProfile profile, int[] seq1, SsccFile sscc,
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
			del[i][0] = Double.NEGATIVE_INFINITY;
		}
		for (int i = 1; i < seq2.length + 1; i++) {
			score[0][i] = 0;
			ins[0][i] = Double.NEGATIVE_INFINITY;
		}
		ins[0][0] = del[0][0] = Double.NEGATIVE_INFINITY;

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
				score[x][y] = Math.max(temp, Math.max(ins[x][y], del[x][y]));
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
		while (x > 0 && y > 0) {
			if (score[x][y] == ins[x][y]) {
				// find the x-k,y position where the gap was opened
				boolean found = false;
				int k = 1;
				while (!found) {
					double diff = score[x - k][y] + profile.getGextend()
							* profile.getGapExtend()[struct2[y]] * (k - 1)
							+ profile.getGopen()
							* profile.getGapInsert()[struct2[y]];
					if (isInEpsilon(diff, score[x][y]) && k < y) {
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
			} else if (score[x][y] == del[x][y]) {
				// find the x-k,y position where the gap was opened
				boolean found = false;
				int k = 1;
				while (!found) {
					double diff = score[x][y - k] + profile.getGextend()
							* profile.getGapExtend()[struct2[y - k]] * (k - 1)
							+ profile.getGopen()
							* profile.getGapInsert()[struct2[y - k]];
					if (isInEpsilon(diff, score[x][y]) && k < x) {
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
		result[0] = "";
		result[1] = "";
		int[] temp = new int[2];
		int[] prev = new int[2];
		prev = tracebackList.pop();
		// global mode
		while (prev[0] == 0 && prev[1] == 0)
			prev = tracebackList.pop();
		if (prev[0] > prev[1] && prev[1] == 0) {
			for (int i = 1; i < prev[0]; i++) {
				result[0] += Character.toString(REVERSE[(char) seq1[i - 1]]);
				result[1] += "-";
			}
		} else if (prev[1] > prev[0] && prev[0] == 0) {
			for (int i = 1; i < prev[1]; i++) {
				result[1] += Character.toString(REVERSE[(char) seq2[i - 1]]);
				result[0] += "-";
			}
		}
		prev = tracebackList.pop();
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
		for (int i = prev[0] + 1; i <= seq1.length; i++) {
			result[0] += Character.toString(REVERSE[(char) seq1[i - 1]]);
			result[1] += "-";
		}
		for (int i = prev[1] + 1; i <= seq2.length; i++) {
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

	public GotohAnswer alignPair() {
		initializeF();
		alignF();
		trace(seq1.length - 1, seq2.length - 1);
		String[] sresult = new String[2];
		sresult = interpretTraceback();
		if (sresult[0].startsWith("-") && sresult[1].startsWith("-")) {
			sresult[0] = sresult[0].substring(1);
			sresult[1] = sresult[1].substring(1);
		}

		GotohAnswer result = new GotohAnswer(seq1ID, seq2ID, sresult[0],
				sresult[1], score[seq1.length][seq2.length], profile);

		// System.out.println(seq1ID + " " + seq2ID + " " + result.getScore());
		return result;
	}

	public String printMatriceHtml(double[][] matrix) {
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
				double printout = Math.round(score[x][y] * 10) / 10.0;
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

	private double getCheck() {
		return this.checkScore;
	}
}
