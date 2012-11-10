package gotoh;

public class GotohProfile {
	private String matrixName, pairfile, seqlibfile;
	private double gopen, gextend, seqWeight, prefWeight, contWeight;
	private String mode, printmatrices, wmode;
	private boolean check, printali;
	public double[][] matrix;

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

	public double getSeqWeight() {
		return seqWeight;
	}

	public void setSeqWeight(double seqWeight) {
		this.seqWeight = seqWeight;
	}

	public double getPrefWeight() {
		return prefWeight;
	}

	public void setPrefWeight(double prefWeight) {
		this.prefWeight = prefWeight;
	}

	public double getContWeight() {
		return contWeight;
	}

	public void setContWeight(double contWeight) {
		this.contWeight = contWeight;
	}

	public String getWmode() {
		return wmode;
	}

	public void setWmode(String wmode) {
		this.wmode = wmode;
	}
	
	
}
