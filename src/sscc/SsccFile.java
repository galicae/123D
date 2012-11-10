package sscc;

public class SsccFile {
	private int[] sequence, localContacts, globalContacts, structure;
	private String ID;

	public SsccFile(int[] sequence, int[] localContacts, int[] globalContacts,
			int[] structure, String ID) {
		super();
		this.sequence = sequence;
		this.localContacts = localContacts;
		this.globalContacts = globalContacts;
		this.structure = structure;
		this.ID = ID;
	}

	public int[] getSequence() {
		return sequence;
	}

	public void setSequence(int[] sequence) {
		this.sequence = sequence;
	}

	public int[] getLocalContacts() {
		return localContacts;
	}

	public void setLocalContacts(int[] contacts) {
		this.localContacts = contacts;
	}

	public int[] getStructure() {
		return structure;
	}

	public void setStructure(int[] structure) {
		this.structure = structure;
	}

	public String getID() {
		return ID;
	}

	public void setID(String iD) {
		ID = iD;
	}

	public int[] getGlobalContacts() {
		return globalContacts;
	}

	public void setGlobalContacts(int[] globalContacts) {
		this.globalContacts = globalContacts;
	}

}
