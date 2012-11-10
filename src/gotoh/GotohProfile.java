package gotoh;

public class GotohProfile {
	private String matrixName, pairfile, seqlibfile;
	private double gopen, gextend;
	private String mode, printmatrices, wmode;
	private boolean check, printali;
	public double[][] matrix;
	private double[] gapInsert = { 2.00, 6.63, 6.05 };
	private double[] gapExtend = { 5.04, 1.66, 2.37 };
	private double[] seqWeight = { 10.57, 7.74, 11.70 };
	private double[] prefWeight = { 4.15, 0.06, 0.06 };
	private double[] localContWeight = { 1.81, 0.05, 0.03 };
	private double[] globalContWeight = { 1.24, 6.44, 1.25 };

	// defaults!!
	public GotohProfile() {
		matrixName = "dayhoff";
		matrix = new double[21][21];
		this.matrix = (double[][]) resources.Matrix.dayhoff.clone();
		gopen = -15;
		gextend = -3;
		mode = "freeshift";
		printali = false;
		printmatrices = "";
		check = false;
	}

	public String getMatrixName() {
		return matrixName;
	}

	public void setMatrix(String matrixName) {
		this.matrixName = matrixName;
		if (matrixName.equals("dayhoff"))
			this.matrix = (double[][]) resources.Matrix.dayhoff.clone();
		else if (matrixName.equals("threader"))
			this.matrix = (double[][]) resources.Matrix.threader.clone();
		else if (matrixName.equals("blakeCohen"))
			this.matrix = (double[][]) resources.Matrix.blakeCohen.clone();
		else if (matrixName.equals("blosum50"))
			this.matrix = (double[][]) resources.Matrix.blosum50.clone();
		else
			this.matrix = (double[][]) resources.Matrix.pam250.clone();
	}

	public String getPairs() {
		return pairfile;
	}

	public void setPairs(String pairs) {
		this.pairfile = pairs;
	}

	public String getSeqlib() {
		return seqlibfile;
	}

	public void setSeqlib(String seqlib) {
		this.seqlibfile = seqlib;
	}

	public double getGopen() {
		return gopen;
	}

	public void setGopen(double gopen) {
		this.gopen = gopen;
	}

	public double getGextend() {
		return gextend;
	}

	public void setGextend(double gextend) {
		this.gextend = gextend;
	}

	public String getMode() {
		return mode;
	}

	public void setMode(String mode) {
		this.mode = mode;
	}

	public String getPrintmatrices() {
		return printmatrices;
	}

	public void setPrintmatrices(String printmatrices) {
		this.printmatrices = printmatrices;
	}

	public boolean isCheck() {
		return check;
	}

	public void setCheck(boolean check) {
		this.check = check;
	}

	public boolean isPrintali() {
		return printali;
	}

	public void setPrintali(boolean printali) {
		this.printali = printali;
	}

	public boolean getPrintali() {
		return this.printali;
	}

	public double getMatrixScore(int x, int y) {
		if (x < y)
			return matrix[y][x];
		else
			return matrix[x][y];
	}

	public double[][] getMatrix() {
		return matrix;
	}

	public String getSeqlibfile() {
		return seqlibfile;
	}

	public void setSeqlibfile(String seqlibfile) {
		this.seqlibfile = seqlibfile;
	}

	public double[] getPrefWeight() {
		return prefWeight;
	}

	public void setPrefWeight(double[] prefWeight) {
		this.prefWeight = prefWeight;
	}

	public String getWmode() {
		return wmode;
	}

	public void setWmode(String wmode) {
		this.wmode = wmode;
	}

	public double[] getSeqWeight() {
		return seqWeight;
	}

	public void setSeqWeight(double[] seqWeight) {
		this.seqWeight = seqWeight;
	}

	public double[] getLocalContWeight() {
		return localContWeight;
	}

	public void setLocalContWeight(double[] localContWeight) {
		this.localContWeight = localContWeight;
	}

	public double[] getGlobalContWeight() {
		return globalContWeight;
	}

	public void setGlobalContWeight(double[] globalContWeight) {
		this.globalContWeight = globalContWeight;
	}

	public double[] getGapInsert() {
		return gapInsert;
	}

	public void setGapInsert(double[] gapInsert) {
		this.gapInsert = gapInsert;
	}

	public double[] getGapExtend() {
		return gapExtend;
	}

	public void setGapExtend(double[] gapExtend) {
		this.gapExtend = gapExtend;
	}

	
}
