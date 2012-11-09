package PSS;

public class PSSMeasurement {

	private int[][] prefArray;
	private int[] occurAmino;
	private int[] occurStruct;
	private int totalAmino;

	public PSSMeasurement(int[][] prefArray, int[] occurAmino, int[] occurStruct,
			int totalAmino) {
		super();
		this.prefArray = prefArray;
		this.occurAmino = occurAmino;
		this.occurStruct = occurStruct;
		this.totalAmino = totalAmino;
	}

	public int[][] getPrefArray() {
		return prefArray;
	}

	public void setPrefArray(int[][] prefArray) {
		this.prefArray = prefArray;
	}

	public int[] getOccurAmino() {
		return occurAmino;
	}

	public void setOccurAmino(int[] occurAmino) {
		this.occurAmino = occurAmino;
	}

	public int[] getOccurStruct() {
		return occurStruct;
	}

	public void setOccurStruct(int[] occurStruct) {
		this.occurStruct = occurStruct;
	}

	public int getTotalAmino() {
		return totalAmino;
	}

	public void setTotalAmino(int totalAmino) {
		this.totalAmino = totalAmino;
	}

}
