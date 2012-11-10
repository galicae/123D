package gotoh;


import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Locale;

public class GotohAnswer {
	private String seq1ID, seq2ID, seq1, seq2;
	private double score;
	private GotohProfile prof;

	public GotohAnswer() {

	}

	public GotohAnswer(String seq1id, String seq2Id, String seq1, String seq2,
			double score, GotohProfile prof) {
		seq1ID = seq1id;
		this.seq2ID = seq2Id;
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.score = Math.round(score * 1000.0) / 1000.0;
		this.prof = prof;
	}

	public String getSeq1ID() {
		return seq1ID;
	}

	public void setSeq1ID(String seq1id) {
		seq1ID = seq1id;
	}

	public String getSeq2Id() {
		return seq2ID;
	}

	public void setSeq2Id(String seq2ID) {
		this.seq2ID = seq2ID;
	}

	public String getSeq1() {
		return seq1;
	}

	public void setSeq1(String seq1) {
		this.seq1 = seq1;
	}

	public String getSeq2() {
		return seq2;
	}

	public void setSeq2(String seq2) {
		this.seq2 = seq2;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}

	public void printAlignment() {
		NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
		DecimalFormat df = (DecimalFormat)nf;
		df.setMinimumFractionDigits(4);
		if (!prof.getPrintmatrices().equals("html")) {
			if (prof.getPrintali()) {
				System.out.println(">" + seq1ID + " " + seq2ID + " " + df.format(score));
				System.out.println(seq1ID + ": " + seq1);
				System.out.println(seq2ID + ": " + seq2);
			} else {
				System.out.println(seq1ID + " " + seq2ID + " " + df.format(score));
			}
		}
		else {
			if (prof.getPrintali()) {
				System.out.println(">" + seq1ID + " " + seq2ID + " " + df.format(score) + "<br>");
				System.out.println(seq1ID + ": " + seq1 + "<br>");
				System.out.println(seq2ID + ": " + seq2 + "<br>");
			} else {
				System.out.println(seq1ID + " " + seq2ID + " " + df.format(score) + "<br>");
			}
		}
	}

}
