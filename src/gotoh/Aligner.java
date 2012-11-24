package gotoh;

import static resources.AminoAcids.REVERSE;

import java.util.LinkedList;

public class Aligner {

	public Aligner() {

	}

	GotohProfile profile;
	int[] seq1, seq2, struct2, localConts, globalConts;
	public double[][] score;
	double[][] ins, del;
	double checkScore = 0;
	double[] max = { 0, 0, 0 };
	LinkedList<int[]> tracebackList;
	String seq1ID, seq2ID;

	public Aligner(GotohProfile profile, int[] seq1, int[] seq2,
			String seq1ID, String seq2ID) {
		this.profile = profile;
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.seq1ID = seq1ID;
		this.seq2ID = seq2ID;
		this.score = new double[seq1.length + 1][seq2.length + 1];
		this.ins = new double[seq1.length + 1][seq2.length + 1];
		this.del = new double[seq1.length + 1][seq2.length + 1];
		tracebackList = new LinkedList<int[]>();
	}

	/**
	 * depending on the mode the initialization of the score matrix is
	 * performed. For a local or free shift alignment the first row and column
	 * are all set to 0; for a global alignment one needs to apply the gap
	 * function.
	 */
	public void initialize() {

	}

	/**
	 * this function aligns the two sequences
	 */
	public void align() {

	}

	public void trace(int x, int y) {

	}

	public double getCheckScore() {
		return 0;
	}

	@Deprecated
	public void printAlignment() {

	}

	public String[] interpretTraceback() {
		return null;
	}

	public GotohAnswer alignPair() {
		return null;
	}

	public String printHtmlPreamble() {
		StringBuilder b = new StringBuilder();
		b.append("<!DOCTYPE html PUBLIC \"...\">");
		b.append("<html>");
		b.append("\t<head>");
		b.append("&nbsp;" + "&nbsp;" + "&nbsp;" + "&nbsp;" + "</head>");
		b.append("&nbsp;" + "&nbsp;" + "&nbsp;" + "&nbsp;" + "<body>"
				+ "<code>");
		return b.toString();
	}

	public String printMatriceHtml(double[][] matrix, int[] seq1, int[] seq2) {
		StringBuilder b = new StringBuilder();
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
				if (matrix[x][y] != Double.NEGATIVE_INFINITY) {
					double printout = Math.round(matrix[x][y] * 10) / 10.0;
					b.append("<td>" + printout + "</td>");
				} else {
					b.append("<td>" + "-inf" + "</td>");
				}
			}
		}
		b.append("</code><br>");
		return b.toString();
	}

	public String printMatriceTxt(double[][] matrix, int[] seq1, int[] seq2) {
		StringBuilder b = new StringBuilder();
		b.append("\t\t");
		for (int i = 0; i < matrix.length - 1; i++)
			b.append(REVERSE[(char) seq1[i]] + "\t");
		int seq2Pointer = 0;
		for (int y = 0; y < matrix[0].length; y++) {
			b.append("\n");
			for (int x = 0; x < matrix.length; x++) {
				if (y == 0 && x == 0) {
					b.append("\t");
				} else if (x == 0) {
					b.append(REVERSE[seq2[seq2Pointer]] + "\t");
					seq2Pointer++;
				}
				double printout = Math.round(matrix[x][y] * 10) / 10.0;
				b.append(printout + "\t");
			}
		}
		b.append("\n");
		return b.toString();
	}

	public int[] getSeq1() {
		return seq1;
	}

	public int[] getSeq2() {
		return seq2;
	}

	public void printMatrices() {
		if (profile.getPrintmatrices().equals("html")) {
			System.out.println(printHtmlPreamble());
			System.out.println("<h2>Score matrix</h2>");
			System.out.println(printMatriceHtml(score, getSeq1(), getSeq2()));
			System.out.println("<h2>Insertion matrix</h2>");
			System.out.println(printMatriceHtml(ins, getSeq1(), getSeq2()));
			System.out.println("<h2>Deletion matrix</h2>");
			System.out.println(printMatriceHtml(del, getSeq1(), getSeq2()));
		} else if (profile.getPrintmatrices().equals("txt")) {
			System.out.println("Score matrix");
			System.out.println(printMatriceTxt(score, getSeq1(), getSeq2()));
			System.out.println("Insertion matrix");
			System.out.println(printMatriceTxt(score, getSeq1(), getSeq2()));
			System.out.println("Deletion matrix");
			System.out.println(printMatriceTxt(score, getSeq1(), getSeq2()));
		}
	}

}
