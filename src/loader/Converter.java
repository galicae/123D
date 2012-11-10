package loader;

import static resources.AminoAcids.*;

public class Converter {
	private static final int ALPHA = 0;
	private static final int BETA = 1;
	private static final int OTHER = 2;

	public int[] convertSeq(String seq) {
		int[] result = new int[seq.length()];
		int i;
		char c;
		for (i = 0; i < seq.length(); i++) {
			c = seq.charAt(i);
			switch (c) {
			case 65:
				result[i] = A;
				continue;
			case 82:
				result[i] = R;
				continue;
			case 78:
				result[i] = N;
				continue;
			case 68:
				result[i] = D;
				continue;
			case 67:
				result[i] = C;
				continue;
			case 81:
				result[i] = Q;
				continue;
			case 69:
				result[i] = E;
				continue;
			case 71:
				result[i] = G;
				continue;
			case 72:
				result[i] = H;
				continue;
			case 73:
				result[i] = I;
				continue;
			case 76:
				result[i] = L;
				continue;
			case 75:
				result[i] = K;
				continue;
			case 77:
				result[i] = M;
				continue;
			case 70:
				result[i] = F;
				continue;
			case 80:
				result[i] = P;
				continue;
			case 83:
				result[i] = S;
				continue;
			case 84:
				result[i] = T;
				continue;
			case 87:
				result[i] = W;
				continue;
			case 89:
				result[i] = Y;
				continue;
			case 86:
				result[i] = V;
				continue;
				// default: result[i] = Integer.MIN_VALUE;
			}
		}

		return result;
	}

	public int[] convertSecStruc(String struct) {
		int[] result = new int[struct.length()];
		int i;
		char c;
		for (i = 0; i < struct.length(); i++) {
			c = struct.charAt(i);
			switch (c) {
			case 97:
				result[i] = ALPHA;
			case 98:
				result[i] = BETA;
			default:
				result[i] = OTHER;
			}
		}

		return result;
	}
}
