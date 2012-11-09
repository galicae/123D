package resources;

public final class AminoAcids {

	// THRee letter code
	public static int ALA = 0;
	public static int CYS = 1;
	public static int ASP = 2;
	public static int GLU = 3;
	public static int PHE = 4;
	public static int GLY = 5;
	public static int HIS = 6;
	public static int ILE = 7;
	public static int LYS = 8;
	public static int LEU = 9;
	public static int MET = 10;
	public static int ASN = 11;
	public static int PRO = 12;
	public static int GLN = 13;
	public static int ARG = 14;
	public static int SER = 15;
	public static int THR = 16;
	public static int VAL = 17;
	public static int TRP = 18;
	public static int TYR = 19;

	// one letter code
	public static int A = 0;
	public static int C = 1;
	public static int D = 2;
	public static int E = 3;
	public static int F = 4;
	public static int G = 5;
	public static int H = 6;
	public static int I = 7;
	public static int K = 8;
	public static int L = 9;
	public static int M = 10;
	public static int N = 11;
	public static int P = 12;
	public static int Q = 13;
	public static int R = 14;
	public static int S = 15;
	public static int T = 16;
	public static int V = 17;
	public static int W = 18;
	public static int Y = 19;

	/**
	 * kinda bad function to make string representation of aminoacid to int
	 * 
	 * @param s
	 *            string representation of aminoacid
	 * @return -1 if input is somehow wrong, the int code of the aminoacid in
	 *         the normal scenario
	 */
	public static int getIntRepresentation(String s) {
		if (s.equals("A")) {
			return A;
		} else if (s.equals("C")) {
			return C;
		} else if (s.equals("D")) {
			return D;
		} else if (s.equals("E")) {
			return E;
		} else if (s.equals("F")) {
			return F;
		} else if (s.equals("G")) {
			return G;
		} else if (s.equals("H")) {
			return H;
		} else if (s.equals("I")) {
			return I;
		} else if (s.equals("K")) {
			return K;
		} else if (s.equals("L")) {
			return L;
		} else if (s.equals("M")) {
			return M;
		} else if (s.equals("N")) {
			return N;
		} else if (s.equals("P")) {
			return P;
		} else if (s.equals("Q")) {
			return Q;
		} else if (s.equals("R")) {
			return R;
		} else if (s.equals("S")) {
			return S;
		} else if (s.equals("T")) {
			return T;
		} else if (s.equals("V")) {
			return V;
		} else if (s.equals("W")) {
			return W;
		} else if (s.equals("Y")) {
			return Y;
		} else if (s.equals("ALA")) {
			return A;
		} else if (s.equals("CYS")) {
			return C;
		} else if (s.equals("ASP")) {
			return D;
		} else if (s.equals("GLU")) {
			return E;
		} else if (s.equals("PHE")) {
			return F;
		} else if (s.equals("GLY")) {
			return G;
		} else if (s.equals("HIS")) {
			return H;
		} else if (s.equals("ILE")) {
			return I;
		} else if (s.equals("LYS")) {
			return K;
		} else if (s.equals("LEU")) {
			return L;
		} else if (s.equals("MET")) {
			return M;
		} else if (s.equals("ASN")) {
			return N;
		} else if (s.equals("PRO")) {
			return P;
		} else if (s.equals("GLN")) {
			return Q;
		} else if (s.equals("ARG")) {
			return R;
		} else if (s.equals("SER")) {
			return S;
		} else if (s.equals("THR")) {
			return T;
		} else if (s.equals("VAL")) {
			return V;
		} else if (s.equals("TRP")) {
			return W;
		} else if (s.equals("TYR")) {
			return Y;
		} else if (Character.isUpperCase(s.codePointAt(0))) {
			return C;
		} else
			return -1;
	}

	public static String[] reverse = { "A", "C", "D", "E", "F", "G", "H", "I",
			"K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" };

}
