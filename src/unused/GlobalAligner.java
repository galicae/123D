package unused;

import static resources.AminoAcids.REVERSE;

import gotoh.Aligner;
import gotoh.GotohAnswer;
import gotoh.GotohProfile;

import java.util.LinkedList;

public class GlobalAligner extends Aligner {
	private GotohProfile profile;
	private String seq1ID, seq2ID;
	private int[] seq1, seq2;
	public double[][] score;
	private double checkScore = 0;
	private double[][] ins, del;
	private LinkedList<int[]> tracebackList;

	public GlobalAligner(GotohProfile profile, int[] seq1, int[] seq2,
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

	public void initialize() {
		score[0][0] = 0;
		for (int i = 1; i < seq1.length + 1; i++) {
			score[i][0] = profile.getGextend() * i + profile.getGopen();
			del[i][0] = Double.NEGATIVE_INFINITY;
		}
		for (int i = 1; i < seq2.length + 1; i++) {
			score[0][i] = profile.getGextend() * i + profile.getGopen();
			ins[0][i] = Double.NEGATIVE_INFINITY;
		}
		ins[0][0] = del[0][0] = Double.NEGATIVE_INFINITY;
	}

	/**
	 * this function aligns the two sequences globally
	 */
	public void align() {
		for (int x = 1; x <= seq1.length; x++) {
			for (int y = 1; y <= seq2.length; y++) {
				double w1 = profile.getGopen() + profile.getGextend();
				ins[x][y] = Math.max(score[x - 1][y] + w1, ins[x - 1][y]
						+ profile.getGextend());
				del[x][y] = Math.max(score[x][y - 1] + w1, del[x][y - 1]
						+ profile.getGextend());
				double temp = score[x - 1][y - 1]
						+ profile.getMatrixScore(seq1[x - 1], seq2[y - 1]);
				score[x][y] = Math.max(temp, Math.max(ins[x][y], del[x][y]));
			}
		}
	}

	public void trace(int x, int y) {
		while (x > 0 && y > 0) {
			if (score[x][y] == ins[x][y]) {
				// find the x-k,y position where the gap was opened
				boolean found = false;
				int k = 1;
				while (!found) {
					double diff = score[x - k][y] + profile.getGextend()
							* (k - 1) + profile.getGopen();
					if (isInEpsilon(diff, score[x][y]) && k < y) {
						int[] res = { x - k, y };
						checkScore += profile.getMatrixScore(seq1[x - k], seq2[y])
								+ profile.getGextend() * (k - 1)
								+ profile.getGopen();
						tracebackList.push(res);
						k++;
					} else {
						int[] res = { x - k, y };
						checkScore += profile.getMatrixScore(seq1[x - k], seq2[y])
								+ profile.getGextend() * (k - 1)
								+ profile.getGopen();
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
							* (k - 1) + profile.getGopen();
					if (isInEpsilon(diff, score[x][y]) && k < x) {
						int[] res = { x, y - k };
						checkScore += profile.getMatrixScore(seq1[x], seq2[y - k])
								+ profile.getGextend() * (k - 1)
								+ profile.getGopen();
						tracebackList.push(res);
						k++;
					} else {
						int[] res = { x, y - k };
						checkScore += profile.getMatrixScore(seq1[x], seq2[y - k])
								+ profile.getGextend() * (k - 1)
								+ profile.getGopen();
						tracebackList.push(res);
						found = true;
						y -= k;
					}
					continue;
				}
			} else {// that means we came from score[x-1][y-1]
				int[] res = { x - 1, y - 1 };
				checkScore += profile.getMatrixScore(seq1[x - 1], seq2[y - 1]);
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
		initialize();
		align();
		trace(seq1.length - 1, seq2.length - 1);
		String[] sresult = new String[2];
		// while (!tracebackList.isEmpty()) {
		// int[] p = new int[2];
		// p = tracebackList.pop();
		// System.out.println(p[0] + ", " + p[1]);
		// }
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

	private static final double epsilon = 0.0001d;

	private static boolean isInEpsilon(double a, double b) {
		return (a > (b - epsilon)) && (a < (b + epsilon));
	}

	private double getCheck() {
		return this.checkScore;
	}
}
