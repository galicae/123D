package sscc;
public class SsccFile {
	private int[] sequence, contacts, structure;
	private String ID;

	public SsccFile(int[] sequence, int[] contacts, int[] structure, String ID) {
		super();
		this.sequence = sequence;
		this.contacts = contacts;
		this.structure = structure;
		this.ID = ID;
	}

	public int[] getSequence() {
		return sequence;
	}

	public void setSequence(int[] sequence) {
		this.sequence = sequence;
	}

	public int[] getContacts() {
		return contacts;
	}

	public void setContacts(int[] contacts) {
		this.contacts = contacts;
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
	
	
	

}
