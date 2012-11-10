package gotoh;

public abstract class Aligner {
	
	public Aligner() {
		
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

	@Deprecated
	public void printAlignment() {
		
	}
	
	public String[] interpretTraceback() {
		return null;
	}
	
	public GotohAnswer alignPair() {
		return null;
	}
	
}
